import cv2
import numpy as np
import sys
import os.path

DIR = "7x1"
INPUT = os.path.join(DIR, "ssDNA_stains_fuzzy")
OUTPUT = os.path.join(DIR, "ssDNA_stains_combined")

NUM_OUTPUTS = 10

# if output exists, clean, otherwise create output folder
if os.path.isdir(OUTPUT):
    for f in os.listdir(OUTPUT):
        os.remove("{0}/{1}".format(OUTPUT, f))
else: os.mkdir(OUTPUT)

# loop through groundtruths
ims = np.array([cv2.imread("{0}/{1}".format(INPUT, f)) for f in os.listdir(INPUT)[:7]])

# set seed for consistant image
np.random.seed(31415)

for r in range(NUM_OUTPUTS):
    # shuffle images
    np.random.shuffle(ims)

    # create im
    result = np.concatenate(ims, axis = 0)
    cv2.imwrite("{0}/ssDNA_{1}.png".format(OUTPUT, r), result)