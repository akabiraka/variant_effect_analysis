

pmd_colors_dict = {"Effect":"darkred", "Knock-out":"red", "No-effect":"darkblue"}
pmd_colors_list = [v for k, v in pmd_colors_dict.items()]
pmd_class_order = [k for k, v in pmd_colors_dict.items()]

patho_colors_dict = {"Pathogenic":"coral", "Likely-pathogenic":"lightcoral", "Neutral":"gray"}
patho_colors_list = [v for k, v in patho_colors_dict.items()]
patho_class_order = [k for k, v in patho_colors_dict.items()]

popu_freq_colors_dict = {"Singleton":"maroon", "Ultra-rare":"chocolate", "Rare":"orange", "Common":"lightsteelblue"}
popu_freq_colors_list = [v for k, v in popu_freq_colors_dict.items()]
popu_freq_class_order = [k for k, v in popu_freq_colors_dict.items()]
popu_freq_class_order_wo_singleton = ["Ultra-rare", "Rare", "Common"]

# pmd_class_order = ["Effect", "Knock-out", "No-effect"]
# patho_class_order = ["Pathogenic", "Likely-pathogenic", "Neutral"]
# popu_freq_class_order = ["Singleton", "Ultra-rare", "Rare", "Common"]
