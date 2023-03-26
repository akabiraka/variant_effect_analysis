
import time
import numpy as np
import torch

#-------- general run test------------------
# from tape import ProteinBertModel, TAPETokenizer
# print(f"Cuda availability: {torch.cuda.is_available()}")

# start = time.time()
# model = ProteinBertModel.from_pretrained('bert-base')
# tokenizer = TAPETokenizer(vocab='iupac')

# sequence = 'GCTVEDRCLIGMGAILLNGCVIGSGSLVAAGALITQ'
# token_ids = torch.tensor(np.array([tokenizer.encode(sequence)]))
# output = model(token_ids)
# print(output[0][0].detach().numpy().shape)

# # from tape import UniRepModel
# # model = UniRepModel.from_pretrained("babbler-1900")
# # tokenizer = TAPETokenizer(vocab='unirep')

# # Note: TAPE (pytorch) does not contain pretrained weights for biLSTM and ResNet model. So deprecated TAPE (TF) needs to be installed. 

# end = time.time()
# print(f"Time taken: {end-start} seconds") # 3.473545789718628 seconds



# ---------------Test whether loaded model is pretrained and generate same outputs everytime-------------
# Protocol: runs twice a model and check visually if they are same
from tape import ProteinBertForMaskedLM, UniRepForLM, TAPETokenizer
# tokenizer = TAPETokenizer(vocab='iupac')
# model = ProteinBertForMaskedLM.from_pretrained('bert-base')

tokenizer = TAPETokenizer(vocab='unirep')
model = UniRepForLM.from_pretrained('babbler-1900')

sequence = 'GCTVEDRC' #l=8
token_ids = torch.tensor(np.array([tokenizer.encode(sequence)]))
logits = model(token_ids)[0][0].detach().numpy() # l+1 x vocab_size=30
print(logits) 
print(logits.shape) 

# Result: ProteinBertForMaskedLM and UniRepForLM generates same result for every run