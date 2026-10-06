#!/bin/bash
# Root submits this array only after phase-I and execution-integrity gates PASS.
#SBATCH --job-name=COPD_V2_internal
#SBATCH --partition=gpu
#SBATCH --gres=gpu:a100:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --array=0-17%4
#SBATCH --output=/vf/users/Dcode/yash/Codex-Manuscript/disease-regulatory-genomics-workflow/diseases/COPD/07_gap_closure/internal-training-1.0/logs/train_%A_%a.log

set -euo pipefail
# Slurm executes a spool copy: never resolve stage relative to BASH_SOURCE.
TREDNET_STAGE=/vf/users/Dcode/yash/Codex-Manuscript/disease-regulatory-genomics-workflow/diseases/COPD/07_gap_closure/internal-training-1.0
TREDNET_ARRAY_INDEX=${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}
if [[ ! "$TREDNET_ARRAY_INDEX" =~ ^([0-9]|1[0-7])$ ]]; then
    echo "Frozen array index must be an integer from 0 through 17" >&2
    exit 2
fi
TREDNET_CONFIGURATIONS=(V2-A V2-B V2-C)
TREDNET_MODELS=(enhancer h3k27me3)
TREDNET_SEEDS=(104729 130363 155921)
TREDNET_CONFIGURATION=${TREDNET_CONFIGURATIONS[$((TREDNET_ARRAY_INDEX / 6))]}
TREDNET_MODEL=${TREDNET_MODELS[$(((TREDNET_ARRAY_INDEX % 6) / 3))]}
TREDNET_SEED=${TREDNET_SEEDS[$((TREDNET_ARRAY_INDEX % 3))]}
TREDNET_RUN_ID=${TREDNET_CONFIGURATION}_${TREDNET_MODEL}_seed${TREDNET_SEED}
export PYTHONHASHSEED="$TREDNET_SEED"
exec bash "$TREDNET_STAGE/scripts/runtime.sh" "$TREDNET_STAGE/scripts/run_logged.py" \
    --stage "$TREDNET_STAGE" --label "${TREDNET_RUN_ID}_attempt001" -- \
    bash "$TREDNET_STAGE/scripts/runtime.sh" "$TREDNET_STAGE/scripts/train_phase_two.py" \
    --stage "$TREDNET_STAGE" --configuration "$TREDNET_CONFIGURATION" \
    --model "$TREDNET_MODEL" --seed "$TREDNET_SEED" --attempt 1
