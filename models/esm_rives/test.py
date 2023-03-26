import torch
import esm
import time

#-------- general run test------------------
# start = time.time()
# model, alphabet = esm.pretrained.esm1_t6_43M_UR50S()
# batch_converter = alphabet.get_batch_converter()
# model.eval()


# data = [("protein1", "MKTVRQERLKSIVRILERSKEPVSGAQLA")]
# batch_labels, batch_strs, batch_tokens = batch_converter(data)

# with torch.no_grad():
#     results = model(batch_tokens, repr_layers=[6], return_contacts=False)

# token_representations = results["representations"][6]

# print(token_representations.shape)
# print(results["logits"][0].shape)
# end = time.time()
# print(f"time taken: {end-start} seconds")


# ---------------Test whether loaded model is pretrained and generate same outputs everytime-------------
# Protocol: runs twice a model and check visually if they are same
model, alphabet = esm.pretrained.esm2_t33_650M_UR50D() #esm1_t6_43M_UR50S, esm1v_t33_650M_UR90S, esm1b_t33_650M_UR50S
batch_converter = alphabet.get_batch_converter()
model.eval()

data = [("protein1", "GCTVEDRC")] #l=8
batch_labels, batch_strs, batch_tokens = batch_converter(data)
logits = model(batch_tokens)["logits"][0].detach().numpy() 
print(logits) 
print(logits.shape) 
# Result: esm1_t6_43M_UR50S, esm1v_t33_650M_UR90S, esm1b_t33_650M_UR50S, esm2_t33_650M_UR50D generates same result for every run

# for esm1v_t33_650M_UR90S: logits shape is l+2 x vocab_size=33
# for esm1b_t33_650M_UR50S: logits shape is l+2 x vocab_size=33
# for esm2_t33_650M_UR50D: logits shape is l+2 x vocab_size=33