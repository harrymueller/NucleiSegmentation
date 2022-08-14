"""
    Get areas of intersections for for the major segment of all true ids
"""

import cv2
import numpy as np
import pandas as pd
import os.path, sys

# CONSTs
MIN_NUM = 25 # minimum area of a true ID

# DIR
DIR = sys.argv[1]
NAME = sys.argv[2]

INPUT_DIR = os.path.join(DIR, NAME, "results/segments")
TRUTH_DIR = os.path.join(DIR, "ground_truths_fixed")
OUTPUT_DIR = os.path.join(DIR, NAME, "results/stats")

# make sure directories made
for f in [OUTPUT_DIR, os.path.join(OUTPUT_DIR, "per_image")]:
    if not os.path.isdir(f): os.mkdir(f)

files = os.listdir(INPUT_DIR)

for (u, f) in enumerate(files):
    print(" > " + f)
    # load true image
    true_f = os.path.join(TRUTH_DIR, f.replace("segments", "truth").replace("csv", "png"))
    true = cv2.imread(true_f)

    # convert RGB colour to code
    true = true.astype(np.uint32)
    true = np.left_shift(true[:,:,0], 16) + np.left_shift(true[:,:,1], 8) + true[:,:,2]

    # load segments
    false_f = os.path.join(INPUT_DIR, f)
    false = np.array(pd.read_csv(false_f, header = None))

    if NAME == "watershed":
        false = false[:,:-1] # remove NaN col

    t_ids = np.unique(true)

    # required measures
    A   = np.zeros((len(t_ids)), dtype = np.float32)
    B   = np.zeros((len(t_ids)), dtype = np.float32)
    AnB = np.zeros((len(t_ids)), dtype = np.float32)

    # set to NaN
    A[:] = np.nan
    B[:] = np.nan
    AnB[:] = np.nan

    # loop through true segments
    for (i, t_id) in enumerate(t_ids):
        # skip black
        if t_id == 0: continue

        # masks
        t_mask = true == t_id

        # remove small colour segments
        if np.count_nonzero(t_mask) < MIN_NUM:
            continue

        f_masked = false[t_mask]

        # false ids
        f_ids = np.unique(f_masked)

        # find the colour that occurs most in the mask
        if len(f_ids) > 1:
            tallies = [np.count_nonzero(f_masked == f_id) if f_id != 1 else 0 for f_id in f_ids]
            f_id = f_ids[np.argmax(tallies)]

        # otherwise use the only one
        else: f_id = f_ids[0]
        
        if f_id == 1: continue # skip if background the major colour
        
        A[i] = np.count_nonzero(t_mask)
        B[i] = np.count_nonzero(false == f_id)
        AnB[i] = np.count_nonzero(np.bitwise_and(t_mask, false == f_id))

    # get the union
    AuB = A + B - AnB

    # save A,B,AnB,AuB to file
    df = pd.DataFrame(np.transpose(np.array([A,B,AnB,AuB])), columns = ['A','B','AnB','AuB'])
    df.to_csv(os.path.join(OUTPUT_DIR, "per_image", f.replace("nuclei", "measures")))