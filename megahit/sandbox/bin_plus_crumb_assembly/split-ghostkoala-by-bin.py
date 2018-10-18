import pandas as pd
        
hu = pd.read_table("all_bins_prokka.ko-ann-full.txt", header=None)
hu[0] = hu[0].replace("user:",  "", regex = True)
hu[6] = hu[0].replace("_[0-9]{5,6}", "", regex = True)

for bin, hu_bin in hu.groupby(6):
    hu_bin.to_csv(f"{bin}.ko-ann-full.txt", index = False, header = False, sep = "\t")
