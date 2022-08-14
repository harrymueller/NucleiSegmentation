"""
    Gets means and STD of F1, IoU, Precision and Recall per image
"""

import numpy as np
import pandas as pd
import os

# DIR
DIR = sys.argv[1]#"/mnt/perkinsdata/tongue_STOmics/benchmarking/25_2k_4"
TOOLS = ["watershed", "cellpose", "deepcell"]
MEASURES = "{dir}/{tool}/results/stats/per_image"
OUTPUT = "{dir}/{tool}/results/stats/measures.csv"

COLS = ["ID", "F1_m", "F1_s", "IoU_m", "IoU_s", "Precision_m", "Precision_s", "Recall_m", "Recall_s"]
TYPES = {i: str if i == "ID" else "float32" for i in COLS}

# loop through tools
for tool in TOOLS:
    output_df = pd.DataFrame(columns = [])

    path = MEASURES.format(dir = DIR, tool = tool)
    measures_df = pd.DataFrame(columns = COLS)

    # loop through image results
    for f in os.listdir(path):
        df = pd.read_csv(os.path.join(path, f), index_col = 0)
        

        # calculate measures
        df["F1"] = 2.0 * df["AnB"] / (df["A"] + df["B"])
        df["IoU"] = df["AnB"] / df["AuB"]
        df["Precision"] = df["AnB"] / df["B"]
        df["Recall"] = df["AnB"] / df["A"]
        
        vals = [f.replace(".csv", "")]
        if len(vals[0]) == 10: vals[0] = vals[0] + " "
        for i in ["F1", "IoU", "Precision", "Recall"]:
            vals.append(np.nanmean(df[i]))
            vals.append(np.nanstd(df[i]))
        measures_df = pd.concat([measures_df, pd.DataFrame([vals], index = [int(vals[0].replace("segments_",""))], columns = COLS)])

    # change types and sort by index
    measures_df = measures_df.astype(TYPES).sort_index()

    measures_df.to_csv(OUTPUT.format(dir = DIR, tool = tool))