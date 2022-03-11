import pandas as pd


myser=pd.Series(["a", "b", "c"], index=["aa","bb","cc"])
print(myser[myser == "b"].index.values)