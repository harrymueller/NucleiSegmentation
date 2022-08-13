import cv2
import numpy as np
import pandas as pd
import os.path, sys

# TODO
# watershed nan last col

# CONSTs
MIN_NUM = 25 # minimum area of a true ID
MAX_DIAM = 64 # maximum width || height of a segment

# DIR
DIR = "/mnt/perkinsdata/tongue_STOmics/benchmarking/25_2k_3"
NAME = sys.argv[1] #"cellpose"
print("#"*30)
print("# " + NAME)
print("#"*30)

INPUT_DIR = os.path.join(DIR, NAME, "results/segments")
TRUTH_DIR = os.path.join(DIR, "ground_truths_fixed")
OUTPUT_DIR = os.path.join(DIR, NAME, "results/stats")

# make sure directories made
for f in [OUTPUT_DIR, os.path.join(OUTPUT_DIR, "per_image")]:
    if not os.path.isdir(f): os.mkdir(f)

files = os.listdir(INPUT_DIR)

all_dice = np.zeros(len(files), dtype = np.float32)
all_ious = np.zeros(len(files), dtype = np.float32)
all_f1   = np.zeros(len(files), dtype = np.float32)

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
        
        if f_id == 1: continue
        A[i] = np.count_nonzero(t_mask)
        B[i] = np.count_nonzero(false == f_id)
        AnB[i] = np.count_nonzero(np.bitwise_and(t_mask, false == f_id))

    # other measures
    AuB = A + B - AnB
    TP  = AnB
    FP  = B - AnB
    FN  = A - AnB

    # final statistics
    dice = 2.0 * AnB / (A + B)
    ious = AnB / AuB 

    precision = TP / (TP + FP)
    recall    = TP / (TP + FN)
    f1        = 2 * precision * recall / (precision + recall)

    # save A,B,AnB,AuB to file
    df = pd.DataFrame(np.transpose(np.array([A,B,AnB,AuB])), columns = ['A','B','AnB','AuB'])
    df.to_csv(os.path.join(OUTPUT_DIR, "per_image", f.replace("nuclei", "measures")))

    all_dice[u] = np.nanmean(dice)
    all_ious[u] = np.nanmean(ious)
    all_f1[u]   = np.nanmean(f1)

df = pd.DataFrame(np.transpose(np.array([all_dice, all_ious, all_f1])), columns = ['Dice','IoUs','F1'])
df.to_csv(os.path.join(OUTPUT_DIR, "measures.csv"))

with open(os.path.join(OUTPUT_DIR, "statistics.csv"), "w") as f:
    f.write("Dice = {:.4f}\nIoU  = {:.4f}\nF1   = {:.4f}\n".format(np.nanmean(all_dice), np.nanmean(all_ious), np.nanmean(all_f1)))