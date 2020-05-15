import numpy as np
import os
import pandas as pd
from scipy.stats import chi2_contingency, fisher_exact
import os

#-----------------------------------------------------------------------------------------------------------------------
# RUN THE MULTIPLE FISHER TEST USING THE CONTINGENCY TABLES:

# Directory to contingency tables:
root = "./contingency_tables"

for ctable in os.listdir(root):
    NAME = ctable.split(".")[0]                             # Name of contingency table
    path = os.path.join(root, ctable)                       # Path to contingency table in iteration
    contingency_df = pd.read_csv(path, index_col=0)         # Read contingency table in as .csv
    data = contingency_df.values                            # Convert pandas object to np.array

    # Values check:
    total_for_disease = data[0, :].sum()
    total_for_nondisease = data[1, :].sum()
    totals = np.array([[total_for_disease], [total_for_nondisease]])

    p_vals = []                                             # list of p-vals for each RBP iteration
    odds = []                                               # List of ORs for each RBP iteration
    for idx, RBP in enumerate(contingency_df.columns):
        print("Test_number:", idx + 1)
        print(RBP)
        RBP_x = contingency_df[f"{RBP}"].values.reshape(2, 1)   # Values for RBP-x CVD and non-CVD
        non_RBP_x = totals - RBP_x                              # Values for all-non-RBP-x for CVD and non-CVD
        full = np.c_[RBP_x, non_RBP_x]                      # Contingency table of RBP-x in iteration vs all-non-RBP-x
        print(full)
        o, p = fisher_exact(full)                           # Fisher test of RBP-x vs all-non-RBP-x
        print("Odds_Ratio:", o)
        odds.append(o)
        print("P_value:", p)
        p_vals.append(p)

    p_vals = np.array(p_vals)                               # Convert list objects to np.arrays
    odds = np.array(odds)
    all_p = np.c_[contingency_df.columns.values.reshape(-1, 1), p_vals.reshape(-1, 1), odds.reshape(-1, 1)]
    all_p = pd.DataFrame(all_p, columns=["RBP", "p_value", "odds_ratio"])       # Convert to pandas dataframe
    all_p.to_csv(f"./fisher_pvalues_{NAME}.csv", index=False, encoding='utf-8') # Save pandas dataframe to .csv



#-----------------------------------------------------------------------------------------------------------------------
