#!/usr/bin/env python3
"""Create validated 2,001-bp REF and ALT sequences for COPD candidates.

The frozen Section 2 candidate table is GRCh38.  Every usable record places
the first base of the normalized VCF allele at zero-based sequence index 1000.
REF and ALT sequences share the same 1,000-bp left flank; their right flanks
are independently cropped so both model inputs remain exactly 2,001 bp.
All input records, including unresolved tags, remain in the metadata audit.
"""

from __future__ import annotations

import hashlib
import json
import platform
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import pysam


ROOT = Path(__file__).resolve().parents[4]
SECTION = ROOT / "diseases" / "COPD" / "04_modeling"
INPUT = (
    ROOT
    / "diseases"
    / "COPD"
    / "02_gwas"
    / "results"
    / "COPD-S2-R006E_candidate_variants_grch38.tsv.gz"
)
FASTA = ROOT / "models" / "TREDNET_v2" / "fasta" / "hg38.fa"
DATA = SECTION / "data"
RESULTS = SECTION / "results"
LOGS = SECTION / "logs"
REF_FASTA = DATA / "COPD_candidate_variants_ref_2001bp.fa"
ALT_FASTA = DATA / "COPD_candidate_variants_alt_2001bp.fa"
METADATA = RESULTS / "COPD-S4-R002_candidate_sequence_audit.tsv.gz"
MANIFEST = RESULTS / "COPD-S4-R002_candidate_sequence_manifest.json"

LEFT_FLANK = 1000
INPUT_LENGTH = 2001


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    RESULTS.mkdir(parents=True, exist_ok=True)
    LOGS.mkdir(parents=True, exist_ok=True)

    candidates = pd.read_csv(INPUT, sep="\t", low_memory=False)
    genome = pysam.FastaFile(str(FASTA))
    references = set(genome.references)
    audit_rows: list[dict[str, object]] = []

    with REF_FASTA.open("w") as ref_handle, ALT_FASTA.open("w") as alt_handle:
        for row in candidates.to_dict(orient="records"):
            record = dict(row)
            candidate_id = str(row["candidate_record_id"])
            status = "ok"
            observed_ref = ""
            ref_sequence = ""
            alt_sequence = ""

            bed_eligible = (
                row["bed_eligible"]
                if isinstance(row["bed_eligible"], bool)
                else str(row["bed_eligible"]).lower() == "true"
            )
            if not bed_eligible:
                status = "missing_coordinate"
            elif pd.isna(row["ref"]) or pd.isna(row["alt"]):
                status = "unresolved_alleles"
            else:
                chrom = str(int(float(row["chromosome_grch38"])))
                fasta_chrom = f"chr{chrom}" if f"chr{chrom}" in references else chrom
                ref = str(row["ref"]).upper()
                alt = str(row["alt"]).upper()
                pos1 = int(float(row["position_grch38"]))
                pos0 = pos1 - 1

                if fasta_chrom not in references:
                    status = "missing_fasta_chromosome"
                elif set(ref + alt) - set("ACGT"):
                    status = "ambiguous_allele"
                elif pos0 - LEFT_FLANK < 0:
                    status = "chromosome_boundary"
                else:
                    observed_ref = genome.fetch(
                        fasta_chrom, pos0, pos0 + len(ref)
                    ).upper()
                    if observed_ref != ref:
                        status = "reference_mismatch"
                    else:
                        left = genome.fetch(
                            fasta_chrom, pos0 - LEFT_FLANK, pos0
                        ).upper()
                        right = genome.fetch(
                            fasta_chrom,
                            pos0 + len(ref),
                            pos0 + len(ref) + INPUT_LENGTH,
                        ).upper()
                        ref_sequence = (left + ref + right)[:INPUT_LENGTH]
                        alt_sequence = (left + alt + right)[:INPUT_LENGTH]
                        if len(ref_sequence) != INPUT_LENGTH or len(alt_sequence) != INPUT_LENGTH:
                            status = "chromosome_boundary"
                        elif set(ref_sequence + alt_sequence) - set("ACGTN"):
                            status = "unsupported_sequence_symbol"

            if status == "ok":
                ref_handle.write(f">{candidate_id}\n{ref_sequence}\n")
                alt_handle.write(f">{candidate_id}\n{alt_sequence}\n")

            record.update(
                {
                    "sequence_status": status,
                    "observed_reference_allele": observed_ref,
                    "sequence_length": INPUT_LENGTH if status == "ok" else pd.NA,
                    "variant_first_base_index_zero_based": LEFT_FLANK if status == "ok" else pd.NA,
                }
            )
            audit_rows.append(record)

    genome.close()
    audit = pd.DataFrame(audit_rows)
    audit.to_csv(
        METADATA,
        sep="\t",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    counts = Counter(audit["sequence_status"])

    manifest = {
        "result_id": "COPD-S4-R002",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "genome_build": "GRCh38",
        "input_length_bp": INPUT_LENGTH,
        "variant_first_base_index_zero_based": LEFT_FLANK,
        "indel_policy": (
            "shared 1000-bp left flank; substitute normalized VCF REF or ALT; "
            "crop independently on the right to 2001 bp"
        ),
        "input_records": len(candidates),
        "status_counts": dict(sorted(counts.items())),
        "software": {
            "python": platform.python_version(),
            "pandas": pd.__version__,
            "pysam": pysam.__version__,
        },
        "files": {},
    }
    for label, path in {
        "candidate_input": INPUT,
        "reference_fasta": FASTA,
        "ref_sequences": REF_FASTA,
        "alt_sequences": ALT_FASTA,
        "sequence_audit": METADATA,
        "script": Path(__file__).resolve(),
    }.items():
        manifest["files"][label] = {
            "path": str(path),
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        }
    MANIFEST.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")

    print(f"Input candidate records: {len(candidates)}")
    for status, count in sorted(counts.items()):
        print(f"{status}: {count}")
    print(f"REF FASTA: {REF_FASTA}")
    print(f"ALT FASTA: {ALT_FASTA}")
    print(f"Audit: {METADATA}")


if __name__ == "__main__":
    main()
