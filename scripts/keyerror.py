import pandas as pd

products = {'Product': ['Tablet','iPhone','Laptop','Monitor'],
            }

df = pd.DataFrame(products)
print(df)
products_list = df.Product.values.tolist()
print (products_list)