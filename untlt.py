

import pandas as pd

df = pd.DataFrame()

df["kmer1"] = [1,2,3,4,5]
df["kmer2"] = [5,4,3,2,1]
df["kmer3"] = [4,3,3,6,1]
pvals = {"kmer1" : 6.66, "kmer2" : 9.99, "kmer3" : 2}

print(df)
df.append(pvals, ignore_index=True)
# df2.index = ['sam1', 'sam2', 'sam3', 'sam4', 'sam5', 'pval']
print(df)
# print(df2)

# df3 = df2.sort_values('pval', axis=1, ascending=True)

# print(df3)

# df4 = df3.drop('pval')

# print(df4)