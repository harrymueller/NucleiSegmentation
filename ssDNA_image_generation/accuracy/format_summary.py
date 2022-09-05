import numpy as np
import pandas as pd
import sys

# DIR
DIR = sys.argv[1]
TOOLS = ["watershed", "cellpose", "deepcell"]
INPUT = "{dir}/accuracy.csv"
OUTPUT = "{dir}/formatted_accuracy_{type}.csv"

# new col names
NEW_COLS = ["ID", "TF||FT", "F1", "IoU", "Precision", "Recall", "Entropy w/", "Entropy w/o"]

# load
summary = pd.read_csv(INPUT.format(dir = DIR))

# capitalise tool names
def cap(s):
    return s.capitalize()
summary["ID"] = summary["ID"].apply(cap)

def rep_tf_ft(x):
    return "True->False" if x == "tf" else "False->True"
summary["TF_FT"] = summary["TF_FT"].apply(rep_tf_ft)

# split into mean and sd
p1 = summary.iloc[:,[0,1,2,4,6,8,10,12]]
p1.columns = NEW_COLS

p2 = summary.iloc[:,[0,1,3,5,7,9,11,13]]
p2.columns = NEW_COLS

# add column with headers
index_df = pd.DataFrame([NEW_COLS], columns = NEW_COLS)

# recombine
new = pd.concat([p1, p2]).sort_index()
new = new.loc[[0,3,1,4,2,5]]

# add mean||std headers
new.insert(2, "Stat", ["Mean", "Std"]*6)

# round to 3dp and add 0s
new = new.reset_index(drop = True).round(3).astype(str)

# add zeros
def add_zeros(x):
    return [a + "0"*(5-len(a)) if i > 2 else a for (i, a) in enumerate(x)]

new = new.apply(add_zeros, raw = True, axis = 1)

# subset and transpose
df_all = new.transpose()
df_tf = new[new["TF||FT"] == "True->False"].drop("TF||FT", axis = 1).transpose()
df_ft = new[new["TF||FT"] == "False->True"].drop("TF||FT", axis = 1).transpose()

# save
df_all.to_csv(OUTPUT.format(dir = DIR, type = "all"), index = True, header = False)
df_tf.to_csv(OUTPUT.format(dir = DIR, type = "tf"), index = True, header = False)
df_ft.to_csv(OUTPUT.format(dir = DIR, type = "ft"), index = True, header = False)