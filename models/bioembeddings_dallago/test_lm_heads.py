import sys
sys.path.append("../variant_effect_analysis")

import time
from models.bioembeddings_dallago.lm_heads.prottrans_lms_factory import load_prottrans_lm_model


# ---------------Test whether loaded model is pretrained and generate same outputs everytime-------------
# Protocol: runs twice a model and check visually if they are same
start = time.time()
# prottrans_bert_bfd, prottrans_albert_bfd, prottrans_xlnet_uniref100, prottrans_t5_bfd, prottrans_t5_uniref50, prottrans_t5_xl_u50
# plus_rnn
model = load_prottrans_lm_model("plus_rnn")  
print(model.name)
print(model._model_directory)
# print(model._model)
# print(model._tokenizer)

logits = model.embed("GCTVEDRC") # it return numpy array. #l=8
print(logits) 
print(logits.shape)
print(f"Time taken to load model: {time.time()-start} s")
# Result: plus_rnn, prottrans_bert_bfd, prottrans_albert_bfd, prottrans_xlnet_uniref100, prottrans_t5_bfd, prottrans_t5_uniref50, prottrans_t5_xl_u50 generates same result for every run

