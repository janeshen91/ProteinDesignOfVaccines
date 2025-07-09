import os
import torch
import numpy as np
import pandas as pd
import h5py
from esm import pretrained
from Bio.PDB import PDBParser

# === Load CSV ===
def load_sequence_pdb_csv(csv_path):
    df = pd.read_csv(csv_path,sep='\t')
    return [(str(row["ID"]), row["sequence"], row["pdb_file"]) for _, row in df.iterrows()]

# === Get attention maps ===
def get_attention_maps(seqs, model_name="esm2_t33_650M_UR50D", device="cuda"):
    model, alphabet = pretrained.load_model_and_alphabet(model_name)
    model = model.to(device)
    model.eval()
    batch_converter = alphabet.get_batch_converter()
    results = {}
    for pid, seq, _ in seqs:
        data = [(pid, seq)]
        _, _, tokens = batch_converter(data)
        tokens = tokens.to(device)
        with torch.no_grad():
            out = model(tokens, repr_layers=[], need_head_weights=True)
            attn = out["attentions"][0][:, :, :len(seq)+2, :len(seq)+2].cpu()  # BOS/EOS
        results[pid] = attn  # [layers, heads, L, L]
    return results

# === Compute contact map from AlphaFold PDB ===
def compute_contact_map(pdb_path, threshold=8.0):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_path)
    ca_coords = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if 'CA' in residue:
                    ca_coords.append(residue['CA'].get_coord())
    ca_coords = np.array(ca_coords)
    dist = np.linalg.norm(ca_coords[:, None, :] - ca_coords[None, :, :], axis=-1)
    contact_map = (dist < threshold).astype(np.uint8)
    np.fill_diagonal(contact_map, 0)
    return contact_map

# === Equation 1 ===
def compute_attention_alignment(attn_maps, contact_map, theta=0.3):
    L, H, N, _ = attn_maps.shape
    contact_tensor = torch.tensor(contact_map).bool().unsqueeze(0).unsqueeze(0)
    high_attn = attn_maps > theta
    numerator = (high_attn & contact_tensor).sum(dim=(-2, -1))
    denominator = high_attn.sum(dim=(-2, -1))
    return (numerator.float() / (denominator.float() + 1e-8)).numpy()  # [L, H]

# === Full pipeline: run and save to HDF5 ===
def run_pipeline_and_save(csv_path, out_path="protein_outputs.h5", theta=0.3, device="cuda"):
    protein_data = load_sequence_pdb_csv(csv_path)
    attn_dict = get_attention_maps(protein_data, device=device)
    with h5py.File(out_path, "w") as h5f:
        for pid, seq, pdb_path in protein_data:
            print(f"Processing {pid}...")
            if pid not in attn_dict or not os.path.exists(pdb_path):
                print(f"  Skipping {pid}: missing attention or PDB")
                continue
            try:
                attn = attn_dict[pid]
                contact = compute_contact_map(pdb_path)
                if contact.shape[0] + 2 != attn.shape[-1]:
                    print(f"  Skipping {pid}: length mismatch")
                    continue
                # Pad contact map for BOS/EOS
                contact_padded = np.pad(contact, pad_width=1, mode='constant')
                alignment = compute_attention_alignment(attn, contact_padded, theta)
                # Save everything
                grp = h5f.create_group(pid)
                grp.create_dataset("attention", data=attn.numpy(), compression="gzip")
                grp.create_dataset("contact_map", data=contact, compression="gzip")
                grp.create_dataset("alignment", data=alignment, compression="gzip")
                grp.attrs["sequence"] = seq
            except Exception as e:
                print(f"  Error processing {pid}: {e}")




import h5py
import matplotlib.pyplot as plt

with h5py.File("protein_outputs.h5", "r") as h5f:
    alignment = h5f["C109_166"]["alignment"][:]  # Replace 1abc with your protein ID

# Plot it
plt.imshow(alignment, aspect="auto", cmap="viridis")
plt.colorbar(label="p_alpha(f)")
plt.title("Equation 1 Alignment (Contact Attention)")
plt.xlabel("Head")
plt.ylabel("Layer")
plt.tight_layout()
plt.savefig("C109_166.png")


