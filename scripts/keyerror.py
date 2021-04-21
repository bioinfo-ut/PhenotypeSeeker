import pandas as pd
import numpy as np

min_cv_inner = np.min(np.bincount([1, 1, 0, 0, 0]))
miin = np.min([
				[9, 1, 3], 
				[4, 5, 6],
				[7, 8, 2],
				], 1)
print(miin)