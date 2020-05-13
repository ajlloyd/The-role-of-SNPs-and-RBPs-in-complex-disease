import numpy as np
import os
import pandas as pd

# load in data in chunks:
mylist = []
for chunk in pd.read_csv("./RBP_SNP_table.tsv", sep="\t", error_bad_lines=False, low_memory=False, chunksize=100000):
    mylist.append(chunk)
df = pd.concat(mylist, axis= 0)

# RBP Binding region (strictly inside):
is_cvd = df["Inside"].str.contains("In")
df_cvd = df[is_cvd].set_index("Inside")
du = pd.crosstab(df_cvd["DISEASE/TRAIT"], df_cvd["Protein"], margins=True)
du.reset_index(level=0, inplace=True)
print(du)
du.to_csv("./RBP_matrix_inside.csv", index=False, encoding='utf-8')


# RBP Binding Region (distance restricted to 20nt, 25nt, and 50nt)
distances = [20, 25, 50]
for n in distances:
    print("Processing distance =", n)
    is_cvd = df["DELTA"] <= n
    df_cvd = df[is_cvd].set_index("DELTA")
    du = pd.crosstab(df_cvd["DISEASE/TRAIT"], df_cvd["Protein"], margins=True)
    du.reset_index(level=0, inplace=True)
    print(du)
    du.to_csv(f"./RBP_matrix_{n}.csv", index=False, encoding='utf-8')
