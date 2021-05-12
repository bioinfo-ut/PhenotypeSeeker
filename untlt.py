

import pandas as pd

df = pd.DataFrame()

df["kmer1"] = [1,2,3,4,5]
df["kmer2"] = [5,4,3,2,1]
pvals = {"kmer1" : 6.66, "kmer2" : 9.99}

df2 = df.append(pvals, ignore_index=True)

print(df2)