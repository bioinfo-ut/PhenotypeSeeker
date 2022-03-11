import pandas as pd


myser=pd.Series(["a", "b", "c"], index=["aa","bb","cc"])
print(myser[myser == "b"].index.values)


lst1 = [1,2,3,4]
lst2 = [5,6,7,8]

print(lst1 + lst2)