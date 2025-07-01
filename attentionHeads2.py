# Visualize Attention Heads at Binding Sites in ESM-2
# Expects a CSV with: sequence,label,binding_sites (comma-separated residue indices)

import torch
from transformers import AutoModelForMaskedLM, AutoTokenizer
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import os
# Load ESM-2 model and tokenizer
model_name = "facebook/esm2_t33_650M_UR50D"
model = AutoModelForMaskedLM.from_pretrained(model_name, output_attentions=True)
tokenizer = AutoTokenizer.from_pretrained(model_name)
model.eval()

# Load sequences and binding site annotations
csv_path = "BLI_bounds.txt" 
csv_data = pd.read_csv(csv_path,sep='\t')
csv_data_cleaned = csv_data.dropna()
sequences = csv_data_cleaned["sequence"].tolist()
binding_sites_list = [list(map(int, s.split(','))) for s in csv_data_cleaned["active_sites"]]

labels = csv_data_cleaned["bound"].tolist()


selected_layer = 28
selected_head = 1
os.chdir("../Layer28Attn1")


###
binding_heatmaps_by_label = {0: [], 1: []}

for idx, (seq, binding_sites, label) in enumerate(zip(sequences, binding_sites_list, labels)):
    if not binding_sites:
        print(f"Skipping sequence {idx} due to empty binding_sites.")
        continue
    print(f"\nSequence {idx}: {seq[:30]}... (length {len(seq)})")
    inputs = tokenizer(seq, return_tensors="pt")
    with torch.no_grad():
        outputs = model(**inputs, output_attentions=True)
        attentions = outputs.attentions
    seq_len = inputs['input_ids'].shape[1]
    # Extract binding-to-binding attention from selected layer and head
    head_attn = attentions[selected_layer][0, selected_head]  # (seq_len, seq_len)
    binding_to_binding = head_attn[binding_sites, :][:, binding_sites].cpu().numpy()  # (num_sites, num_sites)
    # Pad to max size for averaging
    max_sites = 20
    padded = np.full((max_sites, max_sites), np.nan)
    padded[:binding_to_binding.shape[0], :binding_to_binding.shape[1]] = binding_to_binding
    binding_heatmaps_by_label[label].append(padded)

# Average heatmaps
heatmaps_0 = np.array(binding_heatmaps_by_label[0])
heatmaps_1 = np.array(binding_heatmaps_by_label[1])

avg_heatmap_0 = np.nanmean(heatmaps_0, axis=0)
avg_heatmap_1 = np.nanmean(heatmaps_1, axis=0)
diff_heatmap = avg_heatmap_1 - avg_heatmap_0

# Plot
plt.figure(figsize=(6, 5))
sns.heatmap(avg_heatmap_0, cmap="viridis")
plt.title(f"Binding-to-Binding Attention (Unbound)\nLayer {selected_layer}, Head {selected_head}")
plt.tight_layout()
plt.savefig("binding_to_binding_heatmap_label0.png")
plt.close()

plt.figure(figsize=(6, 5))
sns.heatmap(avg_heatmap_1, cmap="viridis")
plt.title(f"Binding-to-Binding Attention (Bound)\nLayer {selected_layer}, Head {selected_head}")
plt.tight_layout()
plt.savefig("binding_to_binding_heatmap_label1.png")
plt.close()

plt.figure(figsize=(6, 5))
sns.heatmap(diff_heatmap, cmap="coolwarm", center=0)
plt.title(f"Difference: Bound - Unbound\nLayer {selected_layer}, Head {selected_head}")
plt.tight_layout()
plt.savefig("binding_to_binding_diff_heatmap.png")
plt.close()





















######################

binding_site_attention_by_label = {0: [], 1: []}
binding_heatmaps_by_label = {0: [], 1: []}

for idx, (seq, binding_sites, label) in enumerate(zip(sequences, binding_sites_list, labels)):
    if not binding_sites:  # Skip sequences with empty binding_sites
        print(f"Skipping sequence {idx} due to empty binding_sites.")
        continue
    print(f"\nSequence {idx}: {seq[:30]}... (length {len(seq)})")
    inputs = tokenizer(seq, return_tensors="pt")
    with torch.no_grad():
        outputs = model(**inputs, output_attentions=True)
        attentions = outputs.attentions  # list of (batch_size, num_heads, seq_len, seq_len)
    seq_len = inputs['input_ids'].shape[1]
    num_layers = len(attentions)
    # Compute average attention FROM binding sites to all positions (i.e., attention focus of binding site residues)
    avg_binding_site_attention = np.zeros((len(binding_sites), seq_len))
    for layer_idx, layer_attn in enumerate(attentions):
        layer_attn = layer_attn[0]  # (num_heads, seq_len, seq_len)
        avg_layer_attn = layer_attn.mean(dim=0)  # (seq_len, seq_len)
        avg_binding_site_attention += avg_layer_attn[binding_sites, :].cpu().numpy()
    avg_binding_site_attention /= num_layers  # normalize by number of layers
    binding_site_attention_by_label[label].append(avg_binding_site_attention)
    # Extract attention heatmap for the selected head/layer
    head_attn = attentions[selected_layer][0, selected_head]  # (seq_len, seq_len)
    binding_heatmap = head_attn[binding_sites, :].cpu().numpy()  # attention FROM binding sites
    # Pad heatmap to fixed width for comparison
    padded_heatmap = np.full((len(binding_sites), 400), np.nan)
    padded_heatmap[:, :seq_len] = binding_heatmap
    binding_heatmaps_by_label[label].append(padded_heatmap)

# Pad heatmaps to uniform shape
max_rows = max([h.shape[0] for h in binding_heatmaps_by_label[0] + binding_heatmaps_by_label[1]])
max_cols = max([h.shape[1] for h in binding_heatmaps_by_label[0] + binding_heatmaps_by_label[1]])

def pad_all_heatmaps(heatmaps, max_rows, max_cols):
    padded = []
    for h in heatmaps:
        new_h = np.full((max_rows, max_cols), np.nan)
        new_h[:h.shape[0], :h.shape[1]] = h
        padded.append(new_h)
    return np.array(padded)

heatmaps_0 = pad_all_heatmaps(binding_heatmaps_by_label[0], max_rows, max_cols)
heatmaps_1 = pad_all_heatmaps(binding_heatmaps_by_label[1], max_rows, max_cols)

# Ensure arrays have same shape before averaging
avg_heatmap_0 = np.nanmean(heatmaps_0, axis=0)
avg_heatmap_1 = np.nanmean(heatmaps_1, axis=0)
diff_heatmap = avg_heatmap_1 - avg_heatmap_0



# Plot heatmaps
plt.figure(figsize=(10, 6))
sns.heatmap(avg_heatmap_0, cmap="viridis")
plt.title(f"Avg Attention FROM Binding Sites (Unbound, Layer {selected_layer}, Head {selected_head})")
plt.xlabel("Sequence Positions")
plt.ylabel("Binding Site Residues")
plt.tight_layout()
plt.savefig("avg_attention_from_binding_heatmap_label0.png")
plt.close()

plt.figure(figsize=(10, 6))
sns.heatmap(avg_heatmap_1, cmap="viridis")
plt.title(f"Avg Attention FROM Binding Sites (Bound, Layer {selected_layer}, Head {selected_head})")
plt.xlabel("Sequence Positions")
plt.ylabel("Binding Site Residues")
plt.tight_layout()
plt.savefig("avg_attention_from_binding_heatmap_label1.png")
plt.close()

plt.figure(figsize=(10, 6))
sns.heatmap(diff_heatmap, cmap="coolwarm", center=0)
plt.title(f"Difference in Attention FROM Binding Sites (Bound - Unbound, Layer {selected_layer}, Head {selected_head})")
plt.xlabel("Sequence Positions")
plt.ylabel("Binding Site Residues")
plt.tight_layout()
plt.savefig("diff_attention_from_binding_heatmap.png")
plt.close()


#####################

residue_attention_by_label = {0: [], 1: []}
binding_heatmaps_by_label = {0: [], 1: []}

for idx, (seq, binding_sites, label) in enumerate(zip(sequences, binding_sites_list, labels)):
    print(f"\nSequence {idx}: {seq[:30]}... (length {len(seq)})")
    inputs = tokenizer(seq, return_tensors="pt")
    with torch.no_grad():
        outputs = model(**inputs, output_attentions=True)
        attentions = outputs.attentions  # list of (batch_size, num_heads, seq_len, seq_len)
    seq_len = inputs['input_ids'].shape[1]
    num_layers = len(attentions)
    num_heads = attentions[0].shape[1]
    # Aggregate attention to binding sites for each residue position
    residue_level_attention = np.zeros((seq_len,))
    for layer_idx, layer_attn in enumerate(attentions):
        layer_attn = layer_attn[0]  # (num_heads, seq_len, seq_len)
        avg_layer_attn = layer_attn.mean(dim=0)  # average across heads -> (seq_len, seq_len)
        residue_level_attention += avg_layer_attn[:, binding_sites].mean(dim=1).cpu().numpy()
    residue_level_attention /= num_layers  # normalize by number of layers
    residue_attention_by_label[label].append(residue_level_attention)
    # Extract attention heatmap to binding sites for selected (layer, head)
    head_attn = attentions[selected_layer][0, selected_head]  # (seq_len, seq_len)
    binding_heatmap = head_attn[:, binding_sites].cpu().numpy()  # attention to binding sites
    binding_heatmaps_by_label[label].append(binding_heatmap)


from scipy.linalg import norm

def pad_to_max_length(arrays, fill_value=np.nan):
    max_len = max(len(arr) for arr in arrays)
    padded = np.full((len(arrays), max_len), fill_value)
    for i, arr in enumerate(arrays):
        padded[i, :len(arr)] = arr
    return padded
    
padded_0 = pad_to_max_length(residue_attention_by_label[0])
padded_1 = pad_to_max_length(residue_attention_by_label[1])


max_len = max(padded_0.shape[1], padded_1.shape[1])

# Pad both to same final length
if padded_0.shape[1] < max_len:
	padded_0 = np.pad(padded_0, ((0, 0), (0, max_len - padded_0.shape[1])), constant_values=np.nan)
if padded_1.shape[1] < max_len:
	padded_1 = np.pad(padded_1, ((0, 0), (0, max_len - padded_1.shape[1])), constant_values=np.nan)
	
if len(residue_attention_by_label[0]) == 0 or len(residue_attention_by_label[1]) == 0:
    raise ValueError("One of the label groups has no sequences. Check your input CSV or labels.")

avg_residue_attn_0 = np.nanmean(padded_0, axis=0)
avg_residue_attn_1 = np.nanmean(padded_1, axis=0)
diff_residue_attn = avg_residue_attn_1 - avg_residue_attn_0

# Plot residue-level attention curves
plt.figure(figsize=(14, 3))
plt.plot(avg_residue_attn_0, label="Unbound (label=0)", color="orange")
plt.plot(avg_residue_attn_1, label="Bound (label=1)", color="green")
plt.title("Residue-Level Attention to Binding Sites")
plt.xlabel("Residue Index")
plt.ylabel("Avg Attention to Binding Sites")
plt.legend()
plt.tight_layout()
plt.savefig("residue_attention_binding_comparison.png")
plt.close()

plt.figure(figsize=(14, 3))
plt.plot(diff_residue_attn, color="blue")
plt.axhline(0, color='black', linestyle='--')
plt.title("Difference in Residue-Level Attention (Bound - Unbound)")
plt.xlabel("Residue Index")
plt.ylabel("Attention Difference")
plt.tight_layout()
plt.savefig("residue_attention_binding_difference.png")
plt.close()

# Compare average heatmaps to binding sites
avg_heatmap_0 = np.mean(binding_heatmaps_by_label[0], axis=0)
avg_heatmap_1 = np.mean(binding_heatmaps_by_label[1], axis=0)
diff_heatmap = avg_heatmap_1 - avg_heatmap_0

plt.figure(figsize=(10, 6))
sns.heatmap(avg_heatmap_0, cmap="viridis")
plt.title(f"Avg Attention to Binding Sites (Unbound, Layer {selected_layer}, Head {selected_head})")
plt.xlabel("Binding Site Residues")
plt.ylabel("Sequence Positions")
plt.tight_layout()
plt.savefig("avg_attention_binding_heatmap_label0.png")
plt.close()

plt.figure(figsize=(10, 6))
sns.heatmap(avg_heatmap_1, cmap="viridis")
plt.title(f"Avg Attention to Binding Sites (Bound, Layer {selected_layer}, Head {selected_head})")
plt.xlabel("Binding Site Residues")
plt.ylabel("Sequence Positions")
plt.tight_layout()
plt.savefig("avg_attention_binding_heatmap_label1.png")
plt.close()

plt.figure(figsize=(10, 6))
sns.heatmap(diff_heatmap, cmap="coolwarm", center=0)
plt.title(f"Difference in Attention (Bound - Unbound, Layer {selected_layer}, Head {selected_head})")
plt.xlabel("Binding Site Residues")
plt.ylabel("Sequence Positions")
plt.tight_layout()
plt.savefig("avg_attention_binding_heatmap_difference.png")
plt.close()
