#!/usr/bin/env python3
"""Download selected donor-matched severe-emphysema ENCODE peak files.

Selection was made after the complete released human-lung metadata inventory in
00_query_encode_lung.py. All nine experiments use donor ENCDO520EJG, whose
ENCODE health-status field states "lungs with severe emphysema". For each of
three lobes, the preferred pseudoreplicated GRCh38 narrowPeak file is retained
for ATAC-seq, H3K27ac, and H3K27me3.
"""

from __future__ import annotations

import gzip
import hashlib
from pathlib import Path

import pandas as pd
import requests


SECTION = Path(__file__).resolve().parents[1]
OUT = SECTION / "data" / "encode_peaks"
OUT.mkdir(parents=True, exist_ok=True)
BASE = "https://www.encodeproject.org"
HEADERS = {"Accept": "application/json"}
ACCESS_DATE = "2026-10-01"

FILES = [
    ("upper_right", "ATAC-seq", "accessibility", "ENCSR990NNX", "ENCFF906HOT"),
    ("upper_right", "Histone ChIP-seq", "H3K27ac", "ENCSR954JMZ", "ENCFF149QUM"),
    ("upper_right", "Histone ChIP-seq", "H3K27me3", "ENCSR337LGY", "ENCFF685CNW"),
    ("lower_right", "ATAC-seq", "accessibility", "ENCSR923VTG", "ENCFF189LZA"),
    ("lower_right", "Histone ChIP-seq", "H3K27ac", "ENCSR837DVF", "ENCFF299FWI"),
    ("lower_right", "Histone ChIP-seq", "H3K27me3", "ENCSR930JRJ", "ENCFF591PDC"),
    ("lower_left", "ATAC-seq", "accessibility", "ENCSR555ZDH", "ENCFF899TQV"),
    ("lower_left", "Histone ChIP-seq", "H3K27ac", "ENCSR875OQJ", "ENCFF014OZD"),
    ("lower_left", "Histone ChIP-seq", "H3K27me3", "ENCSR842LST", "ENCFF725NIC"),
]


def get_json(path: str) -> dict:
    response = requests.get(f"{BASE}{path}", headers=HEADERS, timeout=120)
    response.raise_for_status()
    return response.json()


def download(url: str, output: Path) -> None:
    if output.exists() and output.stat().st_size > 0:
        try:
            with gzip.open(output, "rt") as handle:
                next(handle)
            return
        except (OSError, EOFError, StopIteration):
            output.unlink()
    temporary = output.with_suffix(output.suffix + ".part")
    with requests.get(url, stream=True, timeout=300) as response:
        response.raise_for_status()
        with temporary.open("wb") as handle:
            for block in response.iter_content(chunk_size=1024 * 1024):
                if block:
                    handle.write(block)
    with gzip.open(temporary, "rt") as handle:
        next(handle)
    temporary.replace(output)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    rows: list[dict[str, object]] = []
    donor_accessions: set[str] = set()
    for lobe, assay, mark, experiment, file_accession in FILES:
        metadata = get_json(f"/files/{file_accession}/?format=json")
        if metadata.get("dataset", "").strip("/").split("/")[-1] != experiment:
            raise RuntimeError(f"{file_accession}: unexpected experiment {metadata.get('dataset')}")
        expected = {
            "status": "released",
            "assembly": "GRCh38",
            "file_format": "bed",
            "file_type": "bed narrowPeak",
            "output_type": "pseudoreplicated peaks",
        }
        for key, value in expected.items():
            if metadata.get(key) != value:
                raise RuntimeError(
                    f"{file_accession}: expected {key}={value!r}, found {metadata.get(key)!r}"
                )

        biosamples = metadata.get("biosample_ontology", {})
        experiment_metadata = get_json(f"/experiments/{experiment}/?format=json")
        sample_ids = {
            rep.get("library", {}).get("biosample", {}).get("accession")
            for rep in experiment_metadata.get("replicates", [])
        } - {None}
        donors: set[str] = set()
        health: set[str] = set()
        donor_ages: set[str] = set()
        donor_age_units: set[str] = set()
        donor_sexes: set[str] = set()
        donor_ethnicities: set[str] = set()
        donor_races: set[str] = set()
        for sample_id in sample_ids:
            sample = get_json(f"/biosamples/{sample_id}/?format=json")
            donor = sample.get("donor", {})
            if donor.get("accession"):
                donors.add(donor["accession"])
                donor_accessions.add(donor["accession"])
            if donor.get("health_status"):
                health.add(donor["health_status"])
            if donor.get("age") is not None:
                donor_ages.add(str(donor["age"]))
            if donor.get("age_units"):
                donor_age_units.add(donor["age_units"])
            if donor.get("sex"):
                donor_sexes.add(donor["sex"])
            donor_ethnicities.update(donor.get("ethnicity") or [])
            donor_races.update(donor.get("race") or [])

        if donors != {"ENCDO520EJG"}:
            raise RuntimeError(f"{experiment}: unexpected donors {sorted(donors)}")
        if not any("severe emphysema" in value.lower() for value in health):
            raise RuntimeError(f"{experiment}: donor metadata lacks severe emphysema")

        href = metadata["href"]
        url = BASE + href
        output = OUT / f"{file_accession}.bed.gz"
        download(url, output)
        rows.append(
            {
                "lobe": lobe,
                "assay": assay,
                "mark": mark,
                "experiment": experiment,
                "file_accession": file_accession,
                "biosample_term": biosamples.get("term_name", ""),
                "biosample_accessions": ";".join(sorted(sample_ids)),
                "donor_accessions": ";".join(sorted(donors)),
                "donor_health_status": ";".join(sorted(health)),
                "donor_age": ";".join(sorted(donor_ages)),
                "donor_age_units": ";".join(sorted(donor_age_units)),
                "donor_sex": ";".join(sorted(donor_sexes)),
                "donor_ethnicity": ";".join(sorted(donor_ethnicities)),
                "donor_race": ";".join(sorted(donor_races)),
                "assembly": metadata.get("assembly"),
                "output_type": metadata.get("output_type"),
                "biological_replicates": ";".join(map(str, metadata.get("biological_replicates", []))),
                "technical_replicates": ";".join(metadata.get("technical_replicates", [])),
                "bytes": output.stat().st_size,
                "sha256": sha256(output),
                "path": str(output.resolve()),
                "source_url": url,
                "retrieved_on": ACCESS_DATE,
            }
        )

    if donor_accessions != {"ENCDO520EJG"}:
        raise RuntimeError(f"Peak set is not donor matched: {sorted(donor_accessions)}")
    frame = pd.DataFrame(rows)
    frame.to_csv(OUT / "manifest.tsv", sep="\t", index=False)
    print(frame[["lobe", "mark", "experiment", "file_accession", "bytes"]].to_string(index=False))
    print(f"files\t{len(frame)}")
    print(f"bytes\t{int(frame['bytes'].sum())}")
    print("donor\tENCDO520EJG")


if __name__ == "__main__":
    main()
