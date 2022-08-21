import numpy as np
import pandas as pd
import sys

# DIR
DIR = sys.argv[1]
TOOLS = ["watershed", "cellpose", "deepcell"]
INPUT = "{dir}/accuracy.csv"
OUTPUT = "{dir}/formatted_accuracy.csv"

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

# transpose
new = new.reset_index(drop = True).transpose()

# save
new.to_csv(OUTPUT.format(dir = DIR), index = True, header = False)