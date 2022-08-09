import numpy as np
import cv2 as cv
import os

# consts
ORIG_DIR = "/mnt/perkinsdata/tongue_STOmics/benchmarking/25_2k_2/deepcell/orig"
MASK_DIR = "/mnt/perkinsdata/tongue_STOmics/benchmarking/25_2k_2/deepcell/masks"
OUTPUT_DIR = "/mnt/perkinsdata/tongue_STOmics/benchmarking/25_2k_2/deepcell/results"

RED = np.array([0,0,255])

for orig in os.listdir(ORIG_DIR):
    sample = orig.replace(".png", "")
    
    print("#"*20)
    print(">>> " + sample[7:])
    print("#"*20)

    # load mask array
    print("Reading mask file")
    x = cv.imread(os.path.join(MASK_DIR, sample + "_mask.tif"), cv.IMREAD_ANYDEPTH)

    # load orig image
    print("Loading original image")
    im = cv.imread(os.path.join(ORIG_DIR, sample + ".png"))

    # get outlines of segments
    #contours, hierarchy = cv.findContours(x, cv.RETR_FLOODFILL, cv.CHAIN_APPROX_SIMPLE)

    # apply outlines to im
    print("Creating outlines and fixing mask")
    mask = np.ones_like(x)
    unique_vals = np.unique(x)

    for i in range(x.shape[1]): # x
        for j in range(x.shape[0]): # y
            if x[j,i] != 0: # not background px
                # fix new mask
                #mask[j,i] = np.where(unique_vals == x[j,i])[0][0]
                #break
                if j > 0 and x[j-1][i] != x[j,i]: 
                    im[j,i] = RED
                elif i > 0 and x[j][i-1] != x[j,i]: 
                    im[j,i] = RED
                elif j < (x.shape[0] - 1) and x[j+1][i] != x[j,i]: 
                    im[j,i] = RED
                elif i < (x.shape[1] - 1) and x[j][i+1] != x[j,i]: 
                    im[j,i] = RED

    print("Saving outlines")
    cv.imwrite(os.path.join(OUTPUT_DIR, "outlines", sample + ".png"), im)

    # save masks as segments.csv
    print("Converting masks")
    # VERY INEFFICIENT BUT IT DOES WORK (probably)

    for i in range(1, len(unique_vals)):
        x[x == unique_vals[i]] = i

    x.dtype = np.uint16

    print("Saving masks")
    np.savetxt(os.path.join(OUTPUT_DIR, "segments", sample + ".csv"), mask, delimiter=",", fmt = "%d")