import pandas as pd
import numpy as np


myser=pd.Series(["a", "b", "c"], index=["aa","bb","cc"])
print(myser[myser == "b"].index.values)


lst1 = np.array([1,2,3,4])
lst2 = [5,6,7,8]

print(np.append(lst1, lst2))