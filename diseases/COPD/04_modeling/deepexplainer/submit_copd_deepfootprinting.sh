#!/usr/bin/env bash
set -euo pipefail

# Submit deterministic TREDNet DeepExplainer/TF-MoDISco jobs for COPD models.
# With no arguments, both models are considered; an incomplete model is skipped.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TREDNET_DIR="$(cd "${SCRIPT_DIR}/../trednet" && pwd)"
TEMPLATE="${SCRIPT_DIR}/DeepExplainer_TREDNet_seeded.template.py"
MODEL_I_DIR="${TREDNET_DIR}/model_phase_I"
MODEL_II_ROOT="${TREDNET_DIR}/models_output"
INPUT_DIR="${TREDNET_DIR}/input_training_data"
FASTA_FILE="${TREDNET_DIR}/fasta/hg38.fa"
MEME_FILE="${SCRIPT_DIR}/backup/JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.meme"
DEEPSHAP_ENV="/vf/users/Dcode/gaetano/conda/envs/deepshap_env"
DEEPSHAP_PYTHON="${DEEPSHAP_ENV}/bin/python"
EXPERIMENT="COPD_Section4"
NUM_SEQS_TO_EXPLAIN=1000
BACKGROUND_CONTROLS=100

ENHANCER_EID="COPD_SevereEmphysema_Lung_Enhancer_DHS_x2"
SILENCER_EID="COPD_SevereEmphysema_Lung_Silencer_DHS_x2"

if [[ "$#" -gt 0 ]]; then
    EIDS=("$@")
else
    EIDS=("${ENHANCER_EID}" "${SILENCER_EID}")
fi

COMMON_REQUIRED=(
    "${TEMPLATE}"
    "${MODEL_I_DIR}/phase_one_model.h5"
    "${MODEL_I_DIR}/phase_one_weights.h5"
    "${MODEL_I_DIR}/features_list.txt"
    "${FASTA_FILE}"
    "${MEME_FILE}"
    "${DEEPSHAP_PYTHON}"
)
for required_file in "${COMMON_REQUIRED[@]}"; do
    if [[ ! -f "${required_file}" ]]; then
        echo "ERROR: missing required file: ${required_file}" >&2
        exit 1
    fi
done

mkdir -p "${SCRIPT_DIR}/output/${EXPERIMENT}"

submitted=0
skipped=0
for eid in "${EIDS[@]}"; do
    case "${eid}" in
        "${ENHANCER_EID}"|"${SILENCER_EID}") ;;
        *)
            echo "ERROR: unsupported COPD EID: ${eid}" >&2
            exit 1
            ;;
    esac

    model_dir="${MODEL_II_ROOT}/${eid}"
    positive_bed="${INPUT_DIR}/${eid}_positive_1kb.bed"
    control_bed="${INPUT_DIR}/${eid}_control_1kb.bed"
    model_complete="${model_dir}/auc.txt"
    required_files=(
        "${model_dir}/phase_two_model.keras"
        "${model_dir}/${eid}_phase_two_weights.weights.h5"
        "${positive_bed}"
        "${control_bed}"
    )
    missing=0
    for required_file in "${required_files[@]}"; do
        if [[ ! -f "${required_file}" ]]; then
            echo "WARNING: ${eid} is missing ${required_file}; skipping." >&2
            missing=1
        fi
    done
    if [[ "${missing}" -ne 0 || ! -s "${model_complete}" ]]; then
        echo "SKIP: ${eid} does not yet have a nonempty training-completion marker (${model_complete})."
        skipped=$((skipped + 1))
        continue
    fi

    timestamp="$(date -u +%Y-%m-%dT%H-%M-%SZ)"
    output_dir="${SCRIPT_DIR}/output/${EXPERIMENT}/${eid}_${timestamp}"
    mkdir -p "${output_dir}"
    configured_script="${output_dir}/DeepExplainer_TREDNet.py"
    cp "${TEMPLATE}" "${configured_script}"

    sed -i \
        -e "s|TREDNET_FASTA_FILE|${FASTA_FILE}|g" \
        -e "s|TREDNet_MODEL_I_PATH|${MODEL_I_DIR}|g" \
        -e "s|TREDNet_MODEL_II_PATH|${model_dir}|g" \
        -e "s|INPUT_SEQUENCES_POSITIVE|${positive_bed}|g" \
        -e "s|INPUT_SEQUENCES_CONTROL|${control_bed}|g" \
        -e "s|BIOS_ID|${eid}|g" \
        -e "s|OUTPUT_DIR|${output_dir}|g" \
        -e "s|PLOTTING|False|g" \
        -e "s|NUM_SEQS_TO_EXPLAIN|${NUM_SEQS_TO_EXPLAIN}|g" \
        -e "s|EXPLAIN_FLAG|True|g" \
        -e "s|MEME_PATH|${MEME_FILE}|g" \
        "${configured_script}"

    if grep -Eq 'TREDNET_FASTA_FILE|TREDNet_MODEL_[I]+_PATH|INPUT_SEQUENCES_|BIOS_ID|OUTPUT_DIR|NUM_SEQS_TO_EXPLAIN|MEME_PATH' "${configured_script}"; then
        echo "ERROR: unresolved configuration placeholder in ${configured_script}" >&2
        exit 1
    fi

    job_id="$(
        sbatch --parsable \
            --time=8:00:00 \
            --partition=gpu \
            --gres=gpu:a100:1 \
            --cpus-per-task=4 \
            --qos=global \
            --mem=100g \
            --error="${output_dir}/${eid}_log.err" \
            --output="${output_dir}/${eid}_log.out" \
            --job-name="COPD_${eid#COPD_}_shap" \
            --wrap="export PATH='${DEEPSHAP_ENV}/bin':\$PATH PYTHONHASHSEED=20261001 TF_DETERMINISTIC_OPS=1 MPLBACKEND=Agg; exec '${DEEPSHAP_PYTHON}' '${configured_script}'"
    )"

    printf 'field\tvalue\njob_id\t%s\neid\t%s\nsubmitted_utc\t%s\nseed\t20261001\npositive_sequences_requested\t%s\nbackground_controls\t%s\nmodel_completion_marker\t%s\nconfigured_script\t%s\noutput_directory\t%s\n' \
        "${job_id}" "${eid}" "${timestamp}" "${NUM_SEQS_TO_EXPLAIN}" \
        "${BACKGROUND_CONTROLS}" "${model_complete}" "${configured_script}" "${output_dir}" \
        > "${output_dir}/submission_manifest.tsv"
    echo "SUBMITTED: ${eid} as Slurm job ${job_id}; output ${output_dir}"
    submitted=$((submitted + 1))
done

echo "Submission summary: ${submitted} submitted, ${skipped} skipped."
