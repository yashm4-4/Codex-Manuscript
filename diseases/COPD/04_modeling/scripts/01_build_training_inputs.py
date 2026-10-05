#!/usr/bin/env python3
"""Build reproducible TREDNet enhancer and silencer input partitions."""

from __future__ import annotations

import subprocess
import sys
from collections import Counter
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[4]
MODEL = ROOT / "models" / "TREDNET_v2"
S4 = ROOT / "diseases" / "COPD" / "04_modeling"
PEAKS = S4 / "data" / "training_peaks"
OUT = S4 / "trednet" / "input_training_data"
LOGS = S4 / "logs"
RESULTS = S4 / "results"
OUT.mkdir(parents=True, exist_ok=True)
LOGS.mkdir(parents=True, exist_ok=True)
RESULTS.mkdir(parents=True, exist_ok=True)

SEED = 20261001
RATIO = 2
TRAIN = {f"chr{x}" for x in list(range(1, 7)) + list(range(10, 23)) + ["X", "Y"]}
VALIDATION = {"chr7"}
TEST = {"chr8", "chr9"}

MODELS = {
    "enhancer": {
        "eid": "COPD_SevereEmphysema_Lung_Enhancer_DHS_x2",
        "histone": PEAKS / "COPD_ENCDO520EJG_H3K27ac_all_lobes.narrowPeak.gz",
    },
    "silencer": {
        "eid": "COPD_SevereEmphysema_Lung_Silencer_DHS_x2",
        "histone": PEAKS / "COPD_ENCDO520EJG_H3K27me3_all_lobes.narrowPeak.gz",
    },
}


def partition(chrom: str) -> str:
    if chrom in TRAIN:
        return "train"
    if chrom in VALIDATION:
        return "validation"
    if chrom in TEST:
        return "test"
    return "excluded"


def count_partitions(path: Path) -> Counter[str]:
    counts: Counter[str] = Counter()
    with path.open() as handle:
        for line in handle:
            if line.strip():
                counts[partition(line.split("\t", 1)[0])] += 1
    return counts


def main() -> None:
    accessibility = PEAKS / "COPD_ENCDO520EJG_ATAC_all_lobes.narrowPeak.gz"
    common = [
        "--dnase", str(accessibility),
        "--gencode-gtf", str(ROOT / "data" / "gencode" / "gencode.v50.annotation.gtf.gz"),
        "--blacklist", str(ROOT / "data" / "blacklist" / "hg38-blacklist.v2.bed.gz"),
        "--control-pool", str(MODEL / "input_training_data" / "control" / "allDHS.merge.nonPromoterExon"),
        "--outdir", str(OUT),
        "--ratio", str(RATIO),
        "--seed", str(SEED),
    ]

    rows: list[dict[str, object]] = []
    for model_type, config in MODELS.items():
        command = [
            sys.executable,
            str(MODEL / "make_input_training_data.py"),
            "--eid", config["eid"],
            "--h3k27ac", str(config["histone"]),
            *common,
        ]
        completed = subprocess.run(
            command,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=True,
        )
        (LOGS / f"01_build_{model_type}_inputs.log").write_text(
            "$ " + " ".join(command) + "\n" + completed.stdout
        )

        for label in ("positive", "control"):
            path = OUT / f"{config['eid']}_{label}_1kb.bed"
            counts = count_partitions(path)
            for split in ("train", "validation", "test", "excluded"):
                rows.append(
                    {
                        "model_type": model_type,
                        "eid": config["eid"],
                        "donor": "ENCDO520EJG",
                        "disease_context": "severe emphysema",
                        "class": label,
                        "partition": split,
                        "chromosomes": (
                            "chr1-6,chr10-22,chrX,chrY" if split == "train" else
                            "chr7" if split == "validation" else
                            "chr8-9" if split == "test" else "noncanonical"
                        ),
                        "n_regions": counts[split],
                        "control_ratio": RATIO,
                        "sampling_seed": SEED,
                        "source_accessibility": str(accessibility),
                        "source_histone": str(config["histone"]),
                        "source_control_pool": str(MODEL / "input_training_data" / "control" / "allDHS.merge.nonPromoterExon"),
                    }
                )

    frame = pd.DataFrame(rows)
    frame.to_csv(RESULTS / "COPD-S4-R001_training_partitions.tsv", sep="\t", index=False)
    print(frame.to_string(index=False))


if __name__ == "__main__":
    main()
