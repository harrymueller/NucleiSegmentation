import cv2
import numpy as np
import pandas as pd
import os.path


DIR = "/home/harry/Documents"

# load true image
true_f = os.path.join(DIR, "segments.png")
true = cv2.imread(true_f)

# convert RGB colour to code
true = true.astype(np.uint32)
true = np.left_shift(true[:,:,0], 16) + np.left_shift(true[:,:,1], 8) + true[:,:,2]

# load segments
false_f = os.path.join(DIR, "segments.csv")
false = np.array(pd.read_csv(false_f, header = None))

# remove small colours


t_ids = np.unique(true)

A   = np.ones((len(t_ids)), dtype = np.float32)
B   = np.ones((len(t_ids)), dtype = np.float32)
AnB = np.ones((len(t_ids)), dtype = np.float32)

# loop through true segments
for (i, t_id) in enumerate(t_ids):
    # skip black
    if t_id == 0: continue

    # masks
    t_mask = true == t_id
    f_masked = false[t_mask]

    # false ids
    f_ids = np.unique(f_masked)

    # find the colour that occurs most in the mask
    if len(f_ids) > 1:
        tallies = [np.count_nonzero(f_masked == f_id) for f_id in f_ids]
        f_id = f_ids[np.argmax(tallies)]
    # otherwise use the only one
    else: f_id = f_ids[0]

    A[i] = np.count_nonzero(t_mask)
    B[i] = np.count_nonzero(false == f_id)
    AnB[i] = np.count_nonzero(np.bitwise_and(t_mask, false == f_id))
    
    #print(f_id)
    #dice[i-1] = 2.0 * AnB / (A + B)
    #ious[i-1] = AnB / AuB 

AuB = A + B - AnB
TP  = AnB
FP  = B - AnB
FN  = A - AnB

dice = 2.0 * AnB / (A + B)
ious = AnB / AuB 

precision = TP / (TP + FP)
recall    = TP / (TP + FN)
f1        = 2 * precision * recall / (precision + recall)

print(dice)
print(ious)
print(precision)
print(recall)
print(f1)