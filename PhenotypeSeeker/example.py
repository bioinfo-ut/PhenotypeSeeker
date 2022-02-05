import pandas as pd

mydf = pd.DataFrame().from_dict({'A': ['aa','aa','aa','bb','bb','bb','cc','cc','cc'], 'B':[1,2,3,4,5,6,7,8,9]})
print(mydf)

print(mydf.groupby('A').head(2))