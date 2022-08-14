import numpy as np
import pandas as pd
import sys

# DIR
DIR = sys.argv[1]
TOOLS = ["watershed", "cellpose", "deepcell"]
INPUT = "{dir}/accuracy.csv"
OUTPUT = "{dir}/formatted_accuracy.csv"

# new col names
NEW_COLS = ["ID", "F1", "IoU", "Precision", "Recall"]

# load
summary = pd.read_csv(INPUT.format(dir = DIR))

# capitalise tool names
def cap(s):
    return s.capitalize()
summary["ID"] = summary["ID"].apply(cap)

# split into mean and sd
p1 = summary.iloc[:,[0,1,3,5,7]]
p1.columns = NEW_COLS

p2 = summary.iloc[:,[0,2,4,6,8]]
p2.columns = NEW_COLS

# add column with headers
index_df = pd.DataFrame([NEW_COLS], columns = NEW_COLS)

# recombine
new = pd.concat([index_df, p1,p2]).sort_index()
new.insert(1, "Stat", ["Stat"] + ["Mean", "Std"]*3)
new = new.transpose()

# save
new.to_csv(OUTPUT.format(dir = DIR), index = False, header = False)