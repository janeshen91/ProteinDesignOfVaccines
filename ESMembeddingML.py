import torch
import esm

# Load the pretrained model
model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
batch_converter = alphabet.get_batch_converter()
model.eval()  # disable dropout

# Prepare sequence
data = [("protein1", "MKTFFVLLL")]

batch_labels, batch_strs, batch_tokens = batch_converter(data)

# Extract per-residue embeddings
with torch.no_grad():
    results = model(batch_tokens, repr_layers=[33], return_contacts=False)
token_representations = results["representations"][33]

# Remove BOS/EOS tokens and get embedding per residue
sequence_embedding = token_representations[0, 1:len(data[0][1])+1]
