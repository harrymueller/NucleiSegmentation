import numpy as np
import pandas as pd
import os.path, sys

# DIR
DIR = sys.argv[1]#"/mnt/perkinsdata/tongue_STOmics/benchmarking/25_2k_4"
TOOLS = ["watershed", "cellpose", "deepcell"]
MEASURES = "{dir}/{tool}/results/stats/measures_{tf_ft}.csv"
OUTPUT = "{dir}/accuracy.csv"

COLS = ["ID", "TF_FT", "F1_m", "F1_s", "IoU_m", "IoU_s", "Precision_m", "Precision_s", 
        "Recall_m", "Recall_s", "Entropy_w_m", "Entropy_w_s", "Entropy_wo_m", "Entropy_wo_s"]
TYPES = {i: str if i == "ID" or i == "TF_FT" else "float32" for i in COLS}
summary = pd.DataFrame(columns = COLS)

for tf_ft in ["tf", "ft"]:
    for tool in TOOLS:
        # load data, drop ID col, and convert to 2D np array
        df = np.array(pd.read_csv(MEASURES.format(dir = DIR, tool = tool, tf_ft = tf_ft), index_col = 0).drop(columns = ["ID"]))

        # get means
        means = np.mean(df[:,[0,2,4,6,8,10]], axis = 0)

        # combine stds by taking the sqrt of the mean of the variances
        stds = np.power(df[:,[1,3,5,7,9,11]], 2) # get variances
        stds = np.sqrt(np.sum(stds, axis = 0) / len(df[:,0]))

        #stds = np.hstack([stds, np.std(df[:,[8,9]], axis = 0)])
        vals = [tool, tf_ft, *np.array(list(zip(means,stds))).reshape(len(COLS)-2)]

        summary = pd.concat([summary, pd.DataFrame([vals], columns = COLS)], ignore_index = True)

summary = summary.astype(TYPES)
summary.to_csv(OUTPUT.format(dir = DIR), index = False, float_format='%.4f')