#!/bin/bash

############################### USER CONFIGURATION ###############################
EID="MotorNeuron_Enhancer_DHS_x2"
EXPERIMENT="ALS_Section4"

TREDNET_PATH="/data/Dcode/gaetano/projects/Writing_Scientific_Manuscript_Codex/models/TREDNET_v2"

# Paths to TREDNet models (⚠️ Do NOT add a trailing slash `/` at the end)
TREDNET_MODEL_I_PATH="${TREDNET_PATH}/model_phase_I"  # Pretrained TREDNet I (hg38)
TREDNET_MODEL_II_PATH="${TREDNET_PATH}/models_output/${EID}"  # Fine-tuned TREDNet II

# Input BED files
INPUT_SEQS_POS="${TREDNET_PATH}/input_training_data/${EID}_positive_1kb.bed"
INPUT_SEQS_CTRL="${TREDNET_PATH}/input_training_data/${EID}_control_1kb.bed"

# Enable Plotting of Sequence to Explain (slower if True)
PLOTTING=True # True or False capitalized the first letter

# If False it will run only the Summary plot to get Feature Importance, if True it will run the DeepExplainer and TFmodisco
EXPLAIN_FLAG=True # True or False capitalized the first letter

# Number of Sequences to explain (if 0/negative/greater than # of sequences, it will explain all sequences)
NUM_SEQS_TO_EXPLAIN=1000

# GPU settings
NUM_GPUS=1
TYP_GPU="a100" # v100x: faster to get or a100 for best performance

# MEME FIle
MEME_PATH="/data/Dcode/gaetano/projects/Writing_Scientific_Manuscript_Codex/models/DeepExplainer_TREDNet/backup/JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.meme"
###################################################################################


############################### INPUT VALIDATION #################################
# Fail fast if any file required by DeepExplainer_TREDNet.py is missing,
# instead of failing inside the GPU job after waiting in the queue.
REQUIRED_FILES=(
    "${TREDNET_MODEL_I_PATH}/phase_one_model.h5"
    "${TREDNET_MODEL_I_PATH}/phase_one_weights.h5"
    "${TREDNET_MODEL_I_PATH}/features_list.txt"
    "${TREDNET_MODEL_II_PATH}/phase_two_model.keras"
    "${TREDNET_MODEL_II_PATH}/${EID}_phase_two_weights.weights.h5"
    "${INPUT_SEQS_POS}"
    "${INPUT_SEQS_CTRL}"
    "${MEME_PATH}"
)
MISSING=0
for f in "${REQUIRED_FILES[@]}"; do
    if [ ! -f "$f" ]; then
        echo "❌ Missing required file: $f"
        MISSING=1
    fi
done
if [ "$MISSING" -eq 1 ]; then
    echo "🚫 Job NOT submitted for ${EID}. Fix the missing file(s) above and rerun."
    exit 1
fi
###################################################################################


############################### AUTO SETUP #######################################
# Resolve script directory
SCRIPT_PATH="$(readlink -f "$0")"
SCRIPT_DIR="$(dirname "$SCRIPT_PATH")"

# Create timestamped output directory
TIMESTAMP=$(date +%Y-%m-%d_%H-%M-%S)
OUTPUT_DIR="${SCRIPT_DIR}/output/${EXPERIMENT}/${EID}_${TIMESTAMP}"
mkdir -p "${OUTPUT_DIR}"

# Backup configuration
cp "${SCRIPT_DIR}/backup/args.sh" "${OUTPUT_DIR}/args.sh"
sed -i "s|OUTPUT_DIR|${OUTPUT_DIR}|g" "${OUTPUT_DIR}/args.sh"

# Copy and configure Python script
cp "${SCRIPT_DIR}/backup/DeepExplainer_TREDNet.py" "${OUTPUT_DIR}/DeepExplainer_TREDNet.py"

sed -i "s|TREDNet_MODEL_I_PATH|${TREDNET_MODEL_I_PATH}|g" "${OUTPUT_DIR}/DeepExplainer_TREDNet.py"
sed -i "s|TREDNet_MODEL_II_PATH|${TREDNET_MODEL_II_PATH}|g" "${OUTPUT_DIR}/DeepExplainer_TREDNet.py"
sed -i "s|INPUT_SEQUENCES_POSITIVE|${INPUT_SEQS_POS}|g" "${OUTPUT_DIR}/DeepExplainer_TREDNet.py"
sed -i "s|INPUT_SEQUENCES_CONTROL|${INPUT_SEQS_CTRL}|g" "${OUTPUT_DIR}/DeepExplainer_TREDNet.py"
sed -i "s|BIOS_ID|${EID}|g" "${OUTPUT_DIR}/DeepExplainer_TREDNet.py"
sed -i "s|OUTPUT_DIR|${OUTPUT_DIR}|g" "${OUTPUT_DIR}/DeepExplainer_TREDNet.py"
sed -i "s|PLOTTING|${PLOTTING}|g" "${OUTPUT_DIR}/DeepExplainer_TREDNet.py"
sed -i "s|NUM_SEQS_TO_EXPLAIN|${NUM_SEQS_TO_EXPLAIN}|g" "${OUTPUT_DIR}/DeepExplainer_TREDNet.py"
sed -i "s|EXPLAIN_FLAG|${EXPLAIN_FLAG}|g" "${OUTPUT_DIR}/DeepExplainer_TREDNet.py"
sed -i "s|MEME_PATH|${MEME_PATH}|g" "${OUTPUT_DIR}/DeepExplainer_TREDNet.py"
###################################################################################


############################### LOG OUTPUT #######################################
echo "✅ Job prepared for ${EID} at ${TIMESTAMP} on ${NUM_GPUS} Biowulf GPUs: ${TYP_GPU}"
echo "📂 Output directory: ${OUTPUT_DIR}"
###################################################################################

sbatch --time=8:00:00 \
       --partition=gpu \
       --gres=gpu:${TYP_GPU}:${NUM_GPUS} \
       --cpus-per-task=4 \
       --qos=gpunlm2025.2 \
       --mem=100g \
       --error=${OUTPUT_DIR}/${EID}_log.err \
       --output=${OUTPUT_DIR}/${EID}_log.out \
       --job-name=${EID}_shap \
       ${OUTPUT_DIR}/args.sh
       
