# TREDNet workspace (created 2026-10-01)

Run everything FROM THIS DIRECTORY. Code and the frozen phase-I model are
read from the shared folder; everything you generate is written here.

    SHARED=/vf/users/Dcode/gaetano/projects/Writing_Scientific_Manuscript_Codex/models/TREDNET_v2
    PY=$SHARED/.venv/bin/python

## 1. Build training inputs for a new cell line / tissue
Peaks (DNase/ATAC + H3K27ac narrowPeak, hg38) go in your section's data/
directory first. Then:

    $PY $SHARED/make_input_training_data.py \
      --eid <CellLine>_Enhancer_DHS_x2 \
      --dnase <your_dnase.narrowPeak.gz> \
      --h3k27ac <your_h3k27ac.narrowPeak.gz> \
      --gencode-gtf ../../../../data/gencode/gencode.v50.annotation.gtf.gz \
      --blacklist ../../../../data/blacklist/hg38-blacklist.v2.bed.gz \
      --outdir input_training_data

## 2. Train phase II (GPU node; adapt $SHARED/run_trednet_motorneuron.sh,
##    but cd HERE instead of the shared folder)

    export TREDNET_EID=<CellLine>_Enhancer_DHS_x2
    $PY $SHARED/TREDNet_v2.py        # writes models_output/$TREDNET_EID/

## 3. Score sequences
Your new model:      --models-dir ./models_output
Existing shared models (e.g. MotorNeuron_Enhancer_DHS_x2):
                     --models-dir $SHARED/models_output

    $PY $SHARED/TREDNet_v2_inference.py \
      --eid <EID> --models-dir <see above> \
      --seqs ../data/<variants.fa> --out ../results/<scores.tsv>
