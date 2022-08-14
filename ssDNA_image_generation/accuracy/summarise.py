import numpy as np
import pandas as pd
import os.path, sys

# DIR
DIR = sys.argv[1]#"/mnt/perkinsdata/tongue_STOmics/benchmarking/25_2k_4"
TOOLS = ["watershed", "cellpose", "deepcell"]
MEASURES = "{dir}/{tool}/results/stats/measures.csv"
OUTPUT = "{dir}/accuracy.csv"

COLS = ["ID", "F1_m", "F1_s", "IoU_m", "IoU_s", "Precision_m", "Precision_s", "Recall_m", "Recall_s"]
TYPES = {i: str if i == "ID" else "float32" for i in COLS}

summary = pd.DataFrame(columns = COLS)

for tool in TOOLS:
    # load data, drop ID col, and convert to 2D np array
    df = np.array(pd.read_csv(MEASURES.format(dir = DIR, tool = tool), index_col = 0).drop(columns = ["ID"]))

    # get means
    means = np.mean(df[:,[0,2,4,6]], axis = 0)

    # combine stds by taking the sqrt of the mean of the variances
    stds = np.power(df[:,[1,3,5,7]], 2) # get variances
    stds = np.sqrt(np.sum(stds, axis = 0) / len(df[:,0]))

    vals = [tool, *np.array(list(zip(means,stds))).reshape(8)]
    summary = pd.concat([summary, pd.DataFrame([vals], columns = COLS)], ignore_index = True)

summary = summary.astype(TYPES)
summary.to_csv(OUTPUT.format(dir = DIR), index = False, float_format='%.4f')