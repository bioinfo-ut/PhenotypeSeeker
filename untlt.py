

import pandas as pd

df = pd.DataFrame()

import sys

# df["kmer1"] = [1,2,3,4,5]
# df["kmer2"] = [5,4,3,2,1]
# df["kmer3"] = [4,3,3,6,1]
# pvals = {"kmer1" : 6.66, "kmer2" : 9.99, "kmer3" : 2}

# df2 = df.append(pvals, ignore_index=True)
# df2.index = ['sam1', 'sam2', 'sam3', 'sam4', 'sam5', 'pval']
# print(df2)

# df3 = df2.sort_values('pval', axis=1, ascending=True)

# print(df3)

# df4 = df3.drop('pval')

# print(df4)

# def myfunc():
# 	print("helo")
# 	sys.exit()
# 	print("oleh")

# myfunc()

# print("ulluuka")

# pvals = {"kmer1" : 6.66, "kmer2" : 9.99, "kmer3" : 2}

# if "kmer1" in pvals:
# 	print("Olemas")


# listy = ['aba', 'kaba', 'naba']
# listx = list(range(3))

# listx[listy.index('aba')] = 'prabla'

# print(listx)

mynum = 3.0590232050182594e-07
mynum = '{:.2E}'.format(mynum)
print(mynum)