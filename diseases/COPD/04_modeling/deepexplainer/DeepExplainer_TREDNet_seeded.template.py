#!/usr/bin/env python
# coding: utf-8

# # DeepExplainer-TREDNet

# ## Inputs

# In[1]:

############################## CONFIGURATION ##############################
FASTA_FILE = "TREDNET_FASTA_FILE"

PHASE_ONE_MODEL_FILE = 'TREDNet_MODEL_I_PATH/phase_one_model.h5' ## Not need to change if using hg38
PHASE_ONE_WEIGHTS_FILE = 'TREDNet_MODEL_I_PATH/phase_one_weights.h5' ## Not need to change if using hg38

bios_id = "BIOS_ID" ###### ONLY USED TO SELECT TREDNET P2 AND POSITIVE/CONTROL BED FILES ########

PHASE_TWO_MODEL_FILE = f'TREDNet_MODEL_II_PATH/phase_two_model.keras'
PHASE_TWO_WEIGHTS_FILE = f'TREDNet_MODEL_II_PATH/{bios_id}_phase_two_weights.weights.h5'

POSITIVE_BED_1kb = 'INPUT_SEQUENCES_POSITIVE'
CONTROL_BED_1kb = 'INPUT_SEQUENCES_CONTROL'

# feature list (output TREDNet phase I features)
FEATURES_PATH = "TREDNet_MODEL_I_PATH/features_list.txt" 
save_dir = "OUTPUT_DIR"

MEME_FILE = "MEME_PATH"
SEED = 20261001
############################ END CONFIGURATION ##############################

# ## UTILS

# In[ ]:


import os
import json
import h5py as h5
from tqdm import tqdm
from datetime import datetime
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt

from Bio import SeqIO
from pybedtools import BedTool

import pandas as pd
import numpy as np
import random
import tensorflow as tf

from tensorflow.keras.models import Model
from keras.models import load_model, Sequential
from keras.layers import Conv1D, BatchNormalization, MaxPooling1D, Dropout, Flatten, Dense

import shap
from deeplift.visualization import viz_sequence
from modisco.visualization import viz_sequence
from deeplift.dinuc_shuffle import dinuc_shuffle

os.environ.setdefault("PYTHONHASHSEED", str(SEED))
os.environ.setdefault("TF_DETERMINISTIC_OPS", "1")
random.seed(SEED)
np.random.seed(SEED)
tf.random.set_seed(SEED)
try:
    tf.config.experimental.enable_op_determinism()
except Exception as exc:
    print(f"WARNING: deterministic TensorFlow operations unavailable: {exc}")

# #### FASTA

# In[3]:


def get_chrom2seq(capitalize=True):

    chrom2seq = {}
    for seq in SeqIO.parse(FASTA_FILE, "fasta"):
        chrom2seq[seq.description.split()[0]] = seq.seq.upper() if capitalize else seq.seq

    return chrom2seq

chrom2seq = get_chrom2seq()


# #### from bed to seq

# In[4]:


## load sequence from bed file 
def get_sequences_from_bed(bed_file, chrom2seq):
    bed_list = list(BedTool(bed_file))
    for r in bed_list:
        r.start -= 501
        r.stop += 500

    data_list = [] 
    for r in bed_list:
        _seq = chrom2seq[r.chrom][r.start:r.stop]
        data_list.append(str(_seq))

    return data_list

def get_sequences_position_from_bed(bed_file, position):
    chr, s_e = position.split(":")
    s, e = s_e.split("-")
    s, e = int(s), int(e)

    bed_list = list(BedTool(bed_file))
    for idx, r in enumerate(bed_list):
        if r.chrom == chr:
            if r.start == s-1 and r.stop == e:
                return idx
    return -1


# #### One-Hot Encoder

# In[5]:


def one_hot_encode_along_channel_axis(sequence):
    to_return = np.zeros((len(sequence), 4), dtype=np.int8)
    seq_to_one_hot_fill_in_array(zeros_array=to_return, sequence=sequence, one_hot_axis=1)
    return to_return


def seq_to_one_hot_fill_in_array(zeros_array, sequence, one_hot_axis):
    assert one_hot_axis in (0, 1)
    if one_hot_axis == 0:
        assert zeros_array.shape[1] == len(sequence)
    elif one_hot_axis == 1:
        assert zeros_array.shape[0] == len(sequence)
    # will mutate zeros_array
    for i, char in enumerate(sequence):
        if char in ("A", "a"):
            char_idx = 0
        elif char in ("C", "c"):
            char_idx = 1
        elif char in ("G", "g"):
            char_idx = 2
        elif char in ("T", "t"):
            char_idx = 3
        elif char in ("N", "n"):
            continue  # leave that pos as all 0's
        else:
            raise RuntimeError("Unsupported character: " + str(char))
        if one_hot_axis == 0:
            zeros_array[char_idx, i] = 1
        elif one_hot_axis == 1:
            zeros_array[i, char_idx] = 1

def shuffle_several_times(s):
    s = np.squeeze(s)
    return dinuc_shuffle(s, num_shufs=100)

# ## TREDNet Model

# ### Phase I and II

# In[6]:


def get_phase_one_model():

    print("running get phase 1 model")
    phase_one_model_file = PHASE_ONE_MODEL_FILE
    phase_one_weights_file = PHASE_ONE_WEIGHTS_FILE

    model = load_model(phase_one_model_file)
    model.load_weights(phase_one_weights_file)

    return model

def get_phase_two_model():
    phase_two_weights_file = PHASE_TWO_WEIGHTS_FILE
    # Rebuild instead of deserializing the Keras-3 .keras file: this workflow
    # runs TensorFlow/Keras 2.15 for SHAP compatibility.
    model = Sequential([
        Conv1D(64, 4, activation="relu", padding="valid", input_shape=(4560, 1), name="conv1d_1"),
        BatchNormalization(name="batch_normalization_1"),
        MaxPooling1D(pool_size=2, name="max_pooling1d_1"),
        Dropout(0.4, name="dropout_1"),
        Conv1D(128, 2, activation="relu", padding="valid", name="conv1d_2"),
        Dropout(0.4, name="dropout_2"),
        Flatten(name="flatten_1"),
        Dense(100, activation="relu", name="dense_1"),
        Dense(50, activation="relu", name="dense_2"),
        Dense(1, activation="sigmoid", name="dense_3"),
    ])
    model.load_weights(phase_two_weights_file)
    return model


# ### Merge Phase I and II

# In[7]:

print("----------------------------------")
## get phase I model 
print('Phase I model')
phase_one_model = get_phase_one_model()

## get phase II model
print('Phase II model')
phase_two_model = get_phase_two_model()

## merge phase I and phase II
merged_output = phase_two_model(phase_one_model.output)
merged_model = Model(inputs=phase_one_model.input, outputs=merged_output)
print('Merged model DONE')


# ### Get Data

# #### Get Positive and control bg same model

# In[ ]:


control_seqs = get_sequences_from_bed(CONTROL_BED_1kb,chrom2seq)
positive_seqs = get_sequences_from_bed(POSITIVE_BED_1kb,chrom2seq)

# Remove sequences that are only made up of 'N's (case-insensitive)
control_seqs = [seq for seq in control_seqs if set(seq.upper()) != {'N'}]
control_seqs = [seq for seq in control_seqs if len(seq) == 2001]

# DeepExplainer uses at most 100 background controls. Sampling only those 100
# is distributionally identical to sampling a larger control subset and then
# sub-sampling it, while avoiding a multi-gigabyte unused one-hot matrix.
control_sample_n = min(len(control_seqs), 100)
selected_control_seqs = random.sample(control_seqs, control_sample_n)
ohe_seqs_control = np.array([
    one_hot_encode_along_channel_axis(seq) for seq in selected_control_seqs
])

positive_seqs = [seq for seq in positive_seqs if set(seq.upper()) != {'N'}]
positive_seqs = [seq for seq in positive_seqs if len(seq) == 2001]

if len(positive_seqs) > NUM_SEQS_TO_EXPLAIN and NUM_SEQS_TO_EXPLAIN > 0:
    ohe_seqs_positive = np.array([one_hot_encode_along_channel_axis(seq) for seq in random.sample(positive_seqs, NUM_SEQS_TO_EXPLAIN)])     
else:
    ohe_seqs_positive = np.array([one_hot_encode_along_channel_axis(seq) for seq in positive_seqs])

print('OHE DONE')
print("----------------------------------")

os.makedirs(save_dir, exist_ok=True)
with open(os.path.join(save_dir, "selection_manifest.json"), "w") as handle:
    json.dump(
        {
            "seed": SEED,
            "positive_sequences_available": len(positive_seqs),
            "positive_sequences_explained": len(ohe_seqs_positive),
            "control_sequences_available": len(control_seqs),
            "background_controls": len(ohe_seqs_control),
            "positive_bed": POSITIVE_BED_1kb,
            "control_bed": CONTROL_BED_1kb,
            "fasta": FASTA_FILE,
            "phase_one_weights": PHASE_ONE_WEIGHTS_FILE,
            "phase_two_weights": PHASE_TWO_WEIGHTS_FILE,
            "meme_file": MEME_FILE,
        },
        handle,
        indent=2,
        sort_keys=True,
    )
    handle.write("\n")

# ## Importance Score

# #### Compute Shap Score, same background as TREDNet P.2

# In[9]:

# Compute SHAP values

if len(ohe_seqs_control) < 100:
    background = ohe_seqs_control
else:
    background = ohe_seqs_control[np.random.choice(len(ohe_seqs_control), 100, replace=False)] ## 100 sequences for background per each sequence

explainer_control_bg = shap.DeepExplainer((merged_model.input, merged_model.output[:, 0]), background)


seqs_to_explain = ohe_seqs_positive
batch_size = 64  # You can tune this
if EXPLAIN_FLAG:
    raw_shap_all = []
    shap_explanations = []
    print("----------------------------------")
    print('Number of sequences to explain', len(seqs_to_explain))

    for i in tqdm(range(0, len(seqs_to_explain), batch_size), desc="Explaining sequences"):
        batch = seqs_to_explain[i:i + batch_size]

        raw_batch_shap = explainer_control_bg.shap_values(batch, check_additivity=False).squeeze()
        if raw_batch_shap.ndim == 2:
            raw_batch_shap = raw_batch_shap[np.newaxis, :, :]
        raw_shap_all.append(raw_batch_shap)

        shap_batch_sum = np.sum(raw_batch_shap, axis=-1)
        shap_explanations.append(shap_batch_sum)

    # Combine all batches
    raw_shap_explanations = np.concatenate(raw_shap_all, axis=0)      # shape: (N, L, 4)
    shap_explanation = np.concatenate(shap_explanations, axis=0)      # shape: (N, L)

    print('SHAP DONE')

# ### Save Results

# In[10]:


    # Make the directory if it doesn't exist
    os.makedirs(save_dir, exist_ok=True)

    actual_shap = shap_explanation[:,:,None]*seqs_to_explain

    # Save files
    np.save(os.path.join(save_dir, "raw_shap_explanations.npy"), raw_shap_explanations)
    np.save(os.path.join(save_dir, "actual_shap.npy"), actual_shap)
    np.save(os.path.join(save_dir, "ohe_seqs_to_explain.npy"), seqs_to_explain)

    print(f"Saved SHAP results in: {save_dir}")
    print("----------------------------------")

# ### Visualize the first 3 sequences

# In[ ]:


# Assuming actual_shap and seqs_to_explain are already defined

if PLOTTING & EXPLAIN_FLAG:
    for idx, seq in enumerate(seqs_to_explain[:3]):
        print('Plotting sequence:', idx)
        
        # Create a figure and axis for each plot
        fig = plt.figure(figsize=(20, 5))  # Adjust the size as needed
        ax = fig.add_subplot(111)  # Create a single subplot
        
        # Plot the SHAP values for the current sequence
        viz_sequence.plot_weights_given_ax(ax, actual_shap[idx], subticks_frequency=100)
        
        # Save the plot
        save_path = os.path.join(save_dir, f"shap_sequence_{idx}.png")
        plt.savefig(save_path, bbox_inches='tight', dpi=300)
        plt.show()
        plt.close(fig)  # Close the figure to free memory
        


# ### Visualize first hypothetical score

# In[12]:

if PLOTTING & EXPLAIN_FLAG:
    viz_sequence.plot_weights(raw_shap_explanations[0], subticks_frequency=100)
    print("----------------------------------")

# ## Explain Features TREDNet input P2

# In[13]:


# Step 1: Get inputs to phase_two_model from phase_one_model
input_phase_two = phase_one_model.predict(ohe_seqs_positive)
control_input_phase_two = phase_one_model.predict(background)

new_input_phase_two = input_phase_two[..., np.newaxis]
new_control_input_phase_two = control_input_phase_two[..., np.newaxis]

# Step 2: Create SHAP explainer (try DeepExplainer first; fallback to KernelExplainer if needed)
try:
    explainer = shap.DeepExplainer((phase_two_model.input, phase_two_model.output[:, 0]), new_control_input_phase_two)
except Exception as e:
    print(f"DeepExplainer failed: {e}")

# Step 4: Explain predictions for the first 10 samples
shap_values = explainer.shap_values(new_input_phase_two,check_additivity=False)


# #### Features TREDNet phase I output

# In[14]:


df_features_TREDNet = pd.read_csv(FEATURES_PATH, sep="\t")
df_features_TREDNet['SRC_TYPE_TARGET_CL_NAME'] = df_features_TREDNet['SRC_TYPE'] + "_" + df_features_TREDNet['TARGET']+"_"+df_features_TREDNet['CL_NAME']
df_features_TREDNet


# ### Summary plot

# In[15]:
shap_values_squeezed = np.squeeze(shap_values)  # removes all singleton dimensions
np.save(os.path.join(save_dir, "shap_values_features_TREDNet.npy"), shap_values_squeezed)

gray_red_cmap = LinearSegmentedColormap.from_list("gray_red", ["white", "red"])

# Step 3: Summary plot with top 20 features
# Create summary plot and save to file
plt.figure()  # Start a new figure
shap.summary_plot(
    shap_values_squeezed,
    input_phase_two,
    feature_names=df_features_TREDNet['SRC_TYPE_TARGET_CL_NAME'].values,
    max_display=20,
    cmap=gray_red_cmap,
    show=False  # prevent immediate display
)

# Save to PNG (or PDF/SVG etc.)
plt.savefig(os.path.join(save_dir, "shap_summary_plot_TREDNet.png"), bbox_inches='tight', dpi=300)
plt.close()

# ## TF-MODISCO LITE

# ### OHE Sequences

# In[ ]:

if EXPLAIN_FLAG:
    # Load the data
    data = np.load(f'{save_dir}/ohe_seqs_to_explain.npy')  # shape: (a, b, c)

    # Swap the 2nd and 3rd dimensions
    data_swapped = np.transpose(data, (0, 2, 1))  # shape: (a, c, b)

    # Save as .npz
    np.savez(f'{save_dir}/ohe_tfmodisco.npz', arr_0=data_swapped)


    # ### Hypothetical Scores

    # In[ ]:


    # Load the data
    data = np.load(f'{save_dir}/raw_shap_explanations.npy')  # shape: (a, b, c)

    # Swap the 2nd and 3rd dimensions
    data_swapped = np.transpose(data, (0, 2, 1))  # shape: (a, c, b)

    # Save as .npz
    np.savez(f'{save_dir}/hyp_score_tfmodisco.npz', arr_0=data_swapped)


    # ### Run TF-modisco lite

    # In[ ]:

    print("----------------------------------")
    print('Running TF-modisco')
    os.makedirs(f"{save_dir}/tfmodisco_output", exist_ok=True)
    os.system(
        f"modisco motifs "
        f"-s {save_dir}/ohe_tfmodisco.npz "
        f"-a {save_dir}/hyp_score_tfmodisco.npz "
        f"-n 2000 -o {save_dir}/tfmodisco_output/modisco_results.h5 -w 200 -v"
    )


    # ### Visualize OUTPUT
    os.system(
        f"modisco report "
        f"-i {save_dir}/tfmodisco_output/modisco_results.h5 "
        f"-o {save_dir}/tfmodisco_output/ "
        f"-s ./ "
        f"-m {MEME_FILE}"
    )

    # In[ ]:


    fh = h5.File(f"{save_dir}/tfmodisco_output/modisco_results.h5", "r")
    if "pos_patterns" in fh.keys():
        for akey in fh['pos_patterns'].keys():
            contrib_scores = fh['pos_patterns'][akey]['contrib_scores'][:]

            # Create a figure and axis
            fig = plt.figure(figsize=(20, 5))  # You can adjust the size based on your preference
            ax = fig.add_subplot(111)  # Adds a subplot to the figure

            # Plot the weights on the specified axis
            viz_sequence.plot_weights_given_ax(ax, contrib_scores)  # Pass the axis here

            # Save the figure
            fig.savefig(f"{save_dir}/tfmodisco_output/tf_modisco_pattern_{akey}_contrib_scores.png", dpi=300, bbox_inches='tight')
            plt.show()
            plt.close(fig)

    if "neg_patterns" in fh.keys():
        for akey in fh['neg_patterns'].keys():
            contrib_scores = fh['neg_patterns'][akey]['contrib_scores'][:]

            # Create a figure and axis
            fig = plt.figure(figsize=(20, 5))  # You can adjust the size based on your preference
            ax = fig.add_subplot(111)  # Adds a subplot to the figure

            # Plot the weights on the specified axis
            viz_sequence.plot_weights_given_ax(ax, contrib_scores)  # Pass the axis here

            # Save the figure
            fig.savefig(f"{save_dir}/tfmodisco_output/tf_modisco_pattern_{akey}_contrib_scores_neg.png", dpi=300, bbox_inches='tight')
            plt.show()
            plt.close(fig)

    fh.close()
print("----------------------------------")
print("----------------------------------")
print("             ALL DONE             ")
print("----------------------------------")
print("----------------------------------")
# %%
