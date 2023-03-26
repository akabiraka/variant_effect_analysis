
import sys
sys.path.append("../variant_effect_analysis")

import pickle
  

def save_as_pickle(data, path):
    with open(path, "wb") as f: pickle.dump(data, f)

def load_pickle(path):
    with open(path, "rb") as f: return pickle.load(f)