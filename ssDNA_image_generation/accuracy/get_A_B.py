"""
    Get areas of intersections for for the major segment of all true ids
"""

import cv2
import numpy as np
import pandas as pd
import os.path, sys
from math import log as ln

def log2(x):
    return ln(x) / ln(2)

# CONSTs
MIN_SIZE = 25 # minimum area of a true ID

# DIR
DIR = sys.argv[1]
NAME = sys.argv[2]

INPUT_DIR = os.path.join(DIR, NAME, "results/segments")
TRUTH_DIR = os.path.join(DIR, "ground_truths_fixed")
OUTPUT_DIR = os.path.join(DIR, NAME, "results/stats")

# make sure directories made
for f in [OUTPUT_DIR, os.path.join(OUTPUT_DIR, "per_image"), os.path.join(OUTPUT_DIR, "per_image", "tf"), os.path.join(OUTPUT_DIR, "per_image", "ft")]:
    if not os.path.isdir(f): os.mkdir(f)

files = os.listdir(INPUT_DIR)

def main():
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

        measures_tf = get_areas(true, false)
        measures_ft = get_areas(false, true)

        # save A,B,AnB,AuB to file
        df = pd.DataFrame(np.transpose(np.array(measures_tf)), 
                          columns = ['A', 'B', 'AnB', "ent_w", "ent_wo"])
        df.to_csv(os.path.join(OUTPUT_DIR, "per_image", "tf", f.replace("segments", "measures")))
        df = pd.DataFrame(np.transpose(np.array(measures_ft)), 
                          columns = ['A', 'B', 'AnB', "ent_w", "ent_wo"])
        df.to_csv(os.path.join(OUTPUT_DIR, "per_image", "ft", f.replace("segments", "measures")))
    


def get_areas(A_segments, B_segments):
    """
        Each segment in A is mapped to the largest segment of B in A

        @returns (AreaA, AreaB, AreaIntersection, EntropyWBackground, EntropyWOBackground)
    """
    A_ids = np.unique(A_segments)

    # required measures
    A   = np.zeros((len(A_ids)), dtype = np.float32)
    B   = np.zeros((len(A_ids)), dtype = np.float32)
    AnB = np.zeros((len(A_ids)), dtype = np.float32)

    # entropy per segment
    ent_w = np.zeros((len(A_ids)), dtype = np.float64) # with background
    ent_wo = np.zeros((len(A_ids)), dtype = np.float64) # without background

    # set to NaN
    A[:] = np.nan
    B[:] = np.nan
    AnB[:] = np.nan
    ent_w[:] = np.nan
    ent_wo[:] = np.nan
    
    # loop through segments
    for (i, A_id) in enumerate(A_ids):
        # skip black
        if A_id == 0: continue

        # masks
        A_mask = A_segments == A_id

        # remove small segments
        if np.count_nonzero(A_mask) < MIN_SIZE:
            continue

        B_masked = B_segments[A_mask]

        # false ids
        B_ids = np.unique(B_masked) # sorted

        # find the colour that occurs most in the mask
        if len(B_ids) > 1:
            # count number of pixels per B id
            tallies = [np.count_nonzero(B_masked == B_id) for B_id in B_ids]
            B_id = B_ids[np.argmax(tallies[1:])+1] # ignore the first tally for background

        # otherwise use the only one
        else: B_id = B_ids[0]
        
        if B_id == 1: continue # skip if background the major colour
        
        # area measures
        A[i] = np.count_nonzero(A_mask)
        B[i] = np.count_nonzero(B_segments == B_id)
        AnB[i] = np.count_nonzero(np.bitwise_and(A_mask, B_segments == B_id))

        # calculating entropy
        ent_w[i] = 0
        ent_wo[i] = 0

        for (j, b) in enumerate(B_ids):
            # calc p use MLE
            frac = tallies[j] / A[i]
            ent_w[i] += frac * log2(frac)
            if j != 1: ent_wo[i] += frac * ln(frac)

    # get the union
    return (A, B, AnB, ent_w*(-1), ent_wo*(-1))

if __name__ == "__main__":
    main()