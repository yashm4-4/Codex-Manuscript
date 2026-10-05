import os

# The upstream training entry point does not set a seed.  Keep this local copy
# so the shared, read-only TREDNet installation remains untouched while the
# COPD run has an explicit and recorded random state.
SEED = int(os.getenv("TREDNET_SEED", "20261001"))
os.environ.setdefault("PYTHONHASHSEED", str(SEED))
os.environ.setdefault("TF_DETERMINISTIC_OPS", "1")
os.environ.setdefault("TF_CUDNN_DETERMINISTIC", "1")

import sys
import numpy as np
import random
import time
import glob
import gc
import hashlib
import platform
from Bio import SeqIO
from pybedtools import BedTool
from sklearn import metrics
import h5py

# Disable XLA to prevent large intermediate allocations in Colab
os.environ['TF_XLA_FLAGS'] = '--tf_xla_enable_xla_devices=false'

import tensorflow as tf
import keras
import datetime

from keras.models import load_model
from keras.callbacks import ModelCheckpoint, EarlyStopping
from keras.optimizers import Adadelta

random.seed(SEED)
np.random.seed(SEED)
keras.utils.set_random_seed(SEED)
try:
    tf.config.experimental.enable_op_determinism()
except Exception as exc:
    print(f"WARNING: TensorFlow deterministic operations could not be enabled: {exc}")


train_chromosomes = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr10", "chr11", "chr12", "chr13",
                     "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]
validation_chromosomes = ["chr7"]
test_chromosomes = ["chr8", "chr9"]

INPUT_LENGTH = 2001
EPOCH = int(os.getenv("TREDNET_EPOCHS", "50"))
BATCH_SIZE = int(os.getenv("TREDNET_BATCH_SIZE", "64"))
PRED_BATCH_SIZE = int(os.getenv("TREDNET_PRED_BATCH_SIZE", "32"))
GPUS = 4


nucleotides = ['A', 'C', 'G', 'T']

MODEL_DIR = "./model_phase_I"
MODEL_DIR_phase_II = "./model_phase_II"
FASTA_FILE = "./fasta/hg38.fa"

# Configurable experiment ID and input directory (set via environment variables)
EID = os.getenv("TREDNET_EID", "K562_Enhancer_DHS_x2")
source_input = os.getenv("TREDNET_INPUT_DIR", "./input_training_data/")

SAVE_DIR = f"./models_output/{EID}"


def _sha256(path, chunk_size=1024 * 1024):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        while True:
            block = handle.read(chunk_size)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def write_runtime_manifest():
    """Record the exact run configuration beside the model outputs."""
    os.makedirs(SAVE_DIR, exist_ok=True)
    phase_one_weights = _resolve_existing_file(
        os.path.join(MODEL_DIR, "phase_one_weights"),
        [".weights.h5", ".h5", ".hdf5"],
    )
    positive = os.path.join(source_input, f"{EID}_positive_1kb.bed")
    control = os.path.join(source_input, f"{EID}_control_1kb.bed")
    rows = [
        ("timestamp", datetime.datetime.now(datetime.timezone.utc).isoformat()),
        ("eid", EID),
        ("seed", SEED),
        ("epochs_max", EPOCH),
        ("batch_size", BATCH_SIZE * GPUS),
        ("prediction_batch_size", PRED_BATCH_SIZE),
        ("sequence_stream_batch_size", int(os.getenv("TREDNET_STREAM_BATCH_SIZE", "4096"))),
        ("dataset_builder", "memory_bounded_streaming_v2"),
        ("train_chromosomes", ",".join(train_chromosomes)),
        ("validation_chromosomes", ",".join(validation_chromosomes)),
        ("test_chromosomes", ",".join(test_chromosomes)),
        ("python", sys.version.replace("\n", " ")),
        ("platform", platform.platform()),
        ("tensorflow", tf.__version__),
        ("keras", keras.__version__),
        ("numpy", np.__version__),
        ("devices", ";".join(str(device) for device in tf.config.list_physical_devices())),
        ("phase_one_weights", phase_one_weights),
        ("phase_one_weights_sha256", _sha256(phase_one_weights)),
        ("positive_input", positive),
        ("positive_input_sha256", _sha256(positive)),
        ("control_input", control),
        ("control_input_sha256", _sha256(control)),
        ("training_script", os.path.abspath(__file__)),
        ("training_script_sha256", _sha256(__file__)),
        ("deterministic_ops", os.getenv("TF_DETERMINISTIC_OPS", "")),
    ]
    manifest = os.path.join(SAVE_DIR, "training_manifest.tsv")
    with open(manifest, "w") as handle:
        handle.write("key\tvalue\n")
        for key, value in rows:
            handle.write(f"{key}\t{value}\n")
    print(f"Runtime manifest: {manifest}")
    for key, value in rows:
        if key not in {"python", "platform", "devices"}:
            print(f"CONFIG {key}={value}")


def _resolve_existing_file(base_path, extensions):
    for ext in extensions:
        path = base_path + ext
        if os.path.exists(path):
            return path
    return base_path + extensions[0]


###############################################################################################################################################
def create_dataset():

    print("running create_dataset")
    positive_bed_file = os.path.join(source_input, f"{EID}_positive_1kb.bed")
    print(positive_bed_file)
    negative_bed_file = os.path.join(source_input, f"{EID}_control_1kb.bed")
    print(negative_bed_file)
    print(os.path.join(SAVE_DIR, f"{EID}_phase_two_dataset.hdf5"))
    dataset_save_file = os.path.join(SAVE_DIR, f"{EID}_phase_two_dataset.hdf5")

    create_dataset_for_phase_two(positive_bed_file, negative_bed_file, dataset_save_file)
	
###############################################################################################################################################
def get_chrom2seq(capitalize=True):

    chrom2seq = {}
    for seq in SeqIO.parse(FASTA_FILE, "fasta"):
        chrom_name = seq.description.split()[0]
        seq_value = seq.seq.upper() if capitalize else seq.seq
        chrom2seq[chrom_name] = seq_value

        # Make lookup robust to chr-prefix differences between FASTA and BED files.
        if chrom_name.startswith("chr"):
            chrom2seq[chrom_name[3:]] = seq_value
        else:
            chrom2seq["chr" + chrom_name] = seq_value

    return chrom2seq

###############################################################################################################################################
def seq2one_hot(seq):
    alphabet = np.array(['A','C','G','T'])
    seq = seq.upper()
    char_to_idx = {ch: i for i, ch in enumerate(alphabet)}
    idx = np.array([char_to_idx.get(ch, -1) for ch in seq])
    one_hot = np.zeros((len(seq), 4), dtype=np.float32)
    valid = idx >= 0
    one_hot[np.arange(len(seq))[valid], idx[valid]] = 1
    return one_hot


###############################################################################################################################################
def _build_phase_one_model():
    """Reconstruct the phase-one Sequential CNN in native Keras 3."""
    max_norm = keras.constraints.MaxNorm(max_value=0.9, axis=0)
    glorot_uniform = keras.initializers.GlorotUniform()
    model = keras.Sequential([
        keras.layers.Conv1D(320, 8, activation="relu", padding="valid",
                            kernel_constraint=max_norm, kernel_initializer=glorot_uniform,
                            input_shape=(INPUT_LENGTH, 4), name="conv1d_1"),
        keras.layers.Conv1D(320, 8, activation="relu", padding="valid",
                            kernel_constraint=max_norm, kernel_initializer=glorot_uniform,
                            name="conv1d_2"),
        keras.layers.Dropout(0.2, name="dropout_1"),
        keras.layers.MaxPooling1D(pool_size=6, name="max_pooling1d_1"),
        keras.layers.Conv1D(480, 8, activation="relu", padding="valid",
                            kernel_constraint=max_norm, kernel_initializer=glorot_uniform,
                            name="conv1d_3"),
        keras.layers.Conv1D(480, 8, activation="relu", padding="valid",
                            kernel_constraint=max_norm, kernel_initializer=glorot_uniform,
                            name="conv1d_4"),
        keras.layers.Dropout(0.2, name="dropout_2"),
        keras.layers.MaxPooling1D(pool_size=6, name="max_pooling1d_2"),
        keras.layers.Conv1D(640, 8, activation="relu", padding="valid",
                            kernel_constraint=max_norm, kernel_initializer=glorot_uniform,
                            name="conv1d_5"),
        keras.layers.Conv1D(640, 8, activation="relu", padding="valid",
                            kernel_constraint=max_norm, kernel_initializer=glorot_uniform,
                            name="conv1d_6"),
        keras.layers.Dropout(0.5, name="dropout_3"),
        keras.layers.Flatten(name="flatten_1"),
        keras.layers.Dense(4560, activation="relu", name="dense_1"),
        keras.layers.Dense(4560, activation="linear", name="dense_2"),
        keras.layers.Activation("sigmoid", name="activation_1"),
    ])
    return model


def _build_phase_two_model():
    """Reconstruct the phase-two Sequential CNN in native Keras 3."""
    model = keras.Sequential([
        keras.layers.Conv1D(64, 4, activation="relu", padding="valid",
                            input_shape=(4560, 1), name="conv1d_1"),
        keras.layers.BatchNormalization(name="batch_normalization_1"),
        keras.layers.MaxPooling1D(pool_size=2, name="max_pooling1d_1"),
        keras.layers.Dropout(0.4, name="dropout_1"),
        keras.layers.Conv1D(128, 2, activation="relu", padding="valid", name="conv1d_2"),
        keras.layers.Dropout(0.4, name="dropout_2"),
        keras.layers.Flatten(name="flatten_1"),
        keras.layers.Dense(100, activation="relu", name="dense_1"),
        keras.layers.Dense(50, activation="relu", name="dense_2"),
        keras.layers.Dense(1, activation="sigmoid", name="dense_3"),
    ])
    return model


def get_phase_one_model():

    print("running get phase 1 model")
    json_file = os.path.join(MODEL_DIR, "model.json")
    phase_one_weights_file = _resolve_existing_file(
        os.path.join(MODEL_DIR, "phase_one_weights"),
        [".weights.h5", ".h5", ".hdf5"],
    )

    model = _build_phase_one_model()
    if os.path.exists(phase_one_weights_file):
        model.load_weights(phase_one_weights_file)
    else:
        print(f"WARNING: weights file not found at {phase_one_weights_file}, model will use random weights")

    return model


def predict_in_chunks(model, data, chunk_size=8):
    """Predict as one Keras stream, using ``chunk_size`` as its batch size.

    The upstream helper starts a new ``model.predict`` call for every batch.
    That repeatedly rebuilds a ``tf.data`` pipeline and makes large tissue
    datasets impractically slow.  A single call retains the same numerical
    batching while allowing Keras to stream all batches efficiently.
    """
    return model.predict(data, batch_size=chunk_size, verbose=1)


def _one_hot_batches(bed_file, allowed_chromosomes, chrom2seq, batch_size):
    """Yield bounded one-hot batches from one chromosome partition.

    This preserves the upstream 1-kb-window expansion exactly: subtract 501
    bases from the BED start and add 500 bases to the BED end, producing a
    2,001-bp sequence for every valid input window.
    """
    batch = []
    skipped_missing_chrom = 0
    skipped_length = 0
    with open(bed_file) as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            chrom = fields[0]
            if chrom not in allowed_chromosomes:
                continue
            if chrom not in chrom2seq:
                skipped_missing_chrom += 1
                continue
            start = int(fields[1]) - 501
            stop = int(fields[2]) + 500
            if start < 0:
                skipped_length += 1
                continue
            sequence = chrom2seq[chrom][start:stop]
            if len(sequence) != INPUT_LENGTH:
                skipped_length += 1
                continue
            batch.append(seq2one_hot(sequence))
            if len(batch) == batch_size:
                yield np.stack(batch), skipped_missing_chrom, skipped_length
                batch = []
                skipped_missing_chrom = 0
                skipped_length = 0
    if batch or skipped_missing_chrom or skipped_length:
        matrix = (np.stack(batch) if batch else
                  np.empty((0, INPUT_LENGTH, 4), dtype=np.float32))
        yield matrix, skipped_missing_chrom, skipped_length


def create_dataset_for_phase_two(positive_bed_file, negative_bed_file, dataset_save_file):

    print("running memory-bounded create_dataset_for_phase_two")
    os.makedirs(SAVE_DIR, exist_ok=True)
    chrom2seq = get_chrom2seq()
    model = get_phase_one_model()
    stream_batch_size = int(os.getenv("TREDNET_STREAM_BATCH_SIZE", "4096"))
    partitions = {
        "train": set(train_chromosomes),
        "val": set(validation_chromosomes),
        "test": set(test_chromosomes),
    }
    class_inputs = ((positive_bed_file, 1.0, "positive"),
                    (negative_bed_file, 0.0, "control"))

    print(f"Streaming phase-I embeddings to {dataset_save_file}")
    print(f"Sequence stream batch size: {stream_batch_size}")
    with h5py.File(dataset_save_file, "w") as of:
        of.attrs["dataset_builder"] = "memory_bounded_streaming_v2"
        of.attrs["sequence_stream_batch_size"] = stream_batch_size
        of.attrs["phase_one_prediction_batch_size"] = PRED_BATCH_SIZE
        for split, allowed_chromosomes in partitions.items():
            output = of.create_dataset(
                name=f"{split}_data",
                shape=(0, 4560),
                maxshape=(None, 4560),
                chunks=(256, 4560),
                dtype=np.float32,
                compression="gzip",
                shuffle=True,
            )
            class_counts = {}
            for bed_file, label, class_name in class_inputs:
                n_written = 0
                missing_chrom = 0
                invalid_length = 0
                batches = _one_hot_batches(
                    bed_file,
                    allowed_chromosomes,
                    chrom2seq,
                    stream_batch_size,
                )
                for batch, missing, invalid in batches:
                    missing_chrom += missing
                    invalid_length += invalid
                    if len(batch) == 0:
                        continue
                    embeddings = model.predict(
                        batch,
                        batch_size=PRED_BATCH_SIZE,
                        verbose=0,
                    )
                    previous = output.shape[0]
                    output.resize(previous + len(embeddings), axis=0)
                    output[previous:] = embeddings
                    n_written += len(embeddings)
                    del batch, embeddings
                class_counts[class_name] = n_written
                print(
                    f"{split} {class_name}: {n_written} embeddings; "
                    f"missing_chrom={missing_chrom}; invalid_length={invalid_length}"
                )
                gc.collect()

            labels = np.concatenate((
                np.ones(class_counts["positive"], dtype=np.float32),
                np.zeros(class_counts["control"], dtype=np.float32),
            ))
            if len(labels) != output.shape[0]:
                raise RuntimeError(
                    f"{split}: label/data mismatch {len(labels)} != {output.shape[0]}"
                )
            of.create_dataset(
                name=f"{split}_labels",
                data=labels,
                compression="gzip",
                shuffle=True,
            )
            of.attrs[f"{split}_positive"] = class_counts["positive"]
            of.attrs[f"{split}_control"] = class_counts["control"]
            of.flush()
            print(f"Completed {split}: {output.shape[0]} total embeddings")

    del model, chrom2seq
    gc.collect()
    print(f"Finished streaming dataset: {dataset_save_file}")

###############################################################################################################################################
def train_model():

    print("running train model()")
    import argparse
    import os
    import sys

    data_file = os.path.join(SAVE_DIR, f"{EID}_phase_two_dataset.hdf5")
    model = _build_phase_two_model()
    print("Loading the dataset from:", data_file)
    data = load_dataset(data_file)
    print("Launching the training of model")
    print("Model files and performance evaluation results will be written in:")
    print("        " + SAVE_DIR)
    run_model(data, model, SAVE_DIR)

###############################################################################################################################################
def load_dataset(data_file):

    print("running load dataset()")
    data = {}

    with h5py.File(data_file, "r") as inf:
        for _key in inf:
            data[_key] = inf[_key][()]

    data["train_data"] = data["train_data"][..., np.newaxis]
    data["test_data"] = data["test_data"][..., np.newaxis]
    data["val_data"] = data["val_data"][..., np.newaxis]

    return data

###############################################################################################################################################
def run_model(data, model, save_dir):
    
    print("running run_model()")
    weights_file = os.path.join(SAVE_DIR, f"{EID}_phase_two_weights.weights.h5")
    model_file = os.path.join(SAVE_DIR, "phase_two_model.keras")
	
    os.makedirs(save_dir, exist_ok=True)
    model.save(model_file)

    # Adadelta is recommended to be used with default values
    opt = Adadelta()

    # parallel_model = ModelMGPU(model, gpus=GPUS)
    parallel_model = model
    parallel_model.compile(loss='binary_crossentropy', optimizer=opt, metrics=["accuracy"])

    X_train = data["train_data"]
    Y_train = data["train_labels"]
    X_validation = data["val_data"]
    Y_validation = data["val_labels"]
    X_test = data["test_data"]
    Y_test = data["test_labels"]
    
    _callbacks = []
    checkpoint = ModelCheckpoint(
        filepath=weights_file,
        save_best_only=True,
        save_weights_only=True,
    )
    _callbacks.append(checkpoint)
    earlystopping = EarlyStopping(
        monitor="val_loss",
        patience=15,
        restore_best_weights=True,
    )
    _callbacks.append(earlystopping)
    _callbacks.append(
        keras.callbacks.CSVLogger(os.path.join(SAVE_DIR, "training_history.csv"))
    )

    parallel_model.fit(X_train,
                       Y_train,
                       batch_size=BATCH_SIZE * GPUS,
                       epochs=EPOCH,
                       validation_data=(X_validation, Y_validation),
                       shuffle=True,
                       callbacks=_callbacks,
                       verbose=2)

    # Evaluate and persist the best validation-loss checkpoint, rather than
    # whichever epoch happened to run last.
    parallel_model.load_weights(weights_file)
    parallel_model.save(model_file)
    Y_pred = parallel_model.predict(X_test).ravel()

    auc = metrics.roc_auc_score(Y_test, Y_pred)

    with open(os.path.join(save_dir, "auc.txt"), "w") as of:
        of.write("AUC: %f\n" % auc)

    [fprs, tprs, thrs] = metrics.roc_curve(Y_test, Y_pred)

    sort_ix = np.argsort(np.abs(fprs - 0.1))
    fpr10_thr = thrs[sort_ix[0]]

    sort_ix = np.argsort(np.abs(fprs - 0.05))
    fpr5_thr = thrs[sort_ix[0]]

    sort_ix = np.argsort(np.abs(fprs - 0.03))
    fpr3_thr = thrs[sort_ix[0]]

    sort_ix = np.argsort(np.abs(fprs - 0.01))
    fpr1_thr = thrs[sort_ix[0]]

    with open(os.path.join(save_dir, "fpr_threshold_scores.txt"), "w") as of:
        of.write("10 \t %f\n" % fpr10_thr)
        of.write("5 \t %f\n" % fpr5_thr)
        of.write("3 \t %f\n" % fpr3_thr)
        of.write("1 \t %f\n" % fpr1_thr)

    with open(os.path.join(save_dir, "roc_values.txt"), "w") as of:
        of.write("FPR\tTPR\tTHR\n")
        for fpr, tpr, thr in zip(fprs, tprs, thrs):
            of.write("%f\t%f\t%f\n" % (fpr, tpr, thr))

    [pr, rc, thresholds] = metrics.precision_recall_curve(Y_test, Y_pred)
    auprc = metrics.auc(rc, pr)

    with open(os.path.join(SAVE_DIR, "prc.txt"), "w") as of:
        of.write("PRC: %f\n" % auprc)
        
    with open(os.path.join(SAVE_DIR, "prc_values.txt"), "w") as of:
        of.write("PR\tRC\tTHR\n")
        for pr_, rc_, thr in zip(pr, rc, thresholds):
            of.write("%f\t%f\t%f\n" % (pr_, rc_, thr))

###############################################################################################################################################
if __name__ == "__main__":

    write_runtime_manifest()
    create_dataset()
    train_model()
