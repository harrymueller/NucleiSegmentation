"""
    Gets means and STD of F1, IoU, Precision and Recall per image
"""

import numpy as np
import pandas as pd
import os, sys

# DIR
DIR = sys.argv[1]#"/mnt/perkinsdata/tongue_STOmics/benchmarking/25_2k_4"
TOOLS = ["watershed", "cellpose", "deepcell"]
if len(sys.argv) > 2 and sys.argv[2] == "watershed":
    TOOLS = ["cellpose", "deepcell"]
MEASURES = "{dir}/{tool}/results/stats/per_image/{tf_ft}"
OUTPUT = "{dir}/{tool}/results/stats/measures_{tf_ft}.csv"

COLS = ["ID", "F1_m", "F1_s", "IoU_m", "IoU_s", "Precision_m", "Precision_s", "Recall_m", "Recall_s", "Entropy_w_m", "Entropy_w_s", "Entropy_wo_m", "Entropy_wo_s"]
TYPES = {i: str if i == "ID" else "float32" for i in COLS}

# loop through tools
for tool in TOOLS:
    for tf_ft in ["tf", "ft"]:
        output_df = pd.DataFrame(columns = [])

        path = MEASURES.format(dir = DIR, tool = tool, tf_ft = tf_ft)
        measures_df = pd.DataFrame(columns = COLS)

        # loop through image results
        for f in os.listdir(path):
            df = pd.read_csv(os.path.join(path, f), index_col = 0)

            # calculate measures
            df["F1"] = 2.0 * df["AnB"] / (df["A"] + df["B"])
            df["IoU"] = df["AnB"] / (df["A"] + df["B"] - df["AnB"])
            df["Precision"] = df["AnB"] / df["B"]
            df["Recall"] = df["AnB"] / df["A"]

            # calculated entropy weighted by area
            A_sum = np.sum(df["A"])
            
            vals = [f.replace(".csv", "")]
            if len(vals[0]) == 10: vals[0] = vals[0] + " "

            for i in ["F1", "IoU", "Precision", "Recall", "ent_w", "ent_wo"]:
                vals.append(np.nanmean(df[i]))
                vals.append(np.nanstd(df[i]))

            if False:
                for i in ["F1", "IoU", "Precision", "Recall"]:
                    vals.append(np.nanmean(df[i]))
                    vals.append(np.nanstd(df[i]))

                for i in ["Entropy_w", "Entropy_wo"]:
                    vals.append(np.nanstd(df[i])) # calc standard deviations prior to weighting 

                    df["Entropy_w"] = df["ent_w"] * (df["A"] / A_sum)
                    df["Entropy_wo"] = df["ent_wo"] * (df["A"] / A_sum)
                    vals.append(np.sum(df[i]))

            measures_df = pd.concat([measures_df, pd.DataFrame([vals], index = [int(vals[0].replace("measures_",""))], columns = COLS)])

        # change types and sort by index
        measures_df = measures_df.astype(TYPES).sort_index()

        measures_df.to_csv(OUTPUT.format(dir = DIR, tool = tool, tf_ft = tf_ft))