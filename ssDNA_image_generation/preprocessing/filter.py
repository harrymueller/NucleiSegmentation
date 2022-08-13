import cv2
import numpy as np
import os.path

DIR = "5x5"
INPUT = os.path.join(DIR, "ground_truths")
OUTPUT = os.path.join(DIR, "ground_truths_fixed")

# if output exists, clean, otherwise create output folder
if os.path.isdir(OUTPUT):
    for f in os.listdir(OUTPUT):
        os.remove("{0}/{1}".format(OUTPUT, f))
else: os.mkdir(OUTPUT)

# loop through groundtruths
for f in os.listdir(INPUT):
    im = cv2.imread("{0}/{1}".format(INPUT, f))

    # any pixels less than 128 are multiplied by 2
    im[im<128] = im[im<128] * 2
    cv2.imwrite("{0}/{1}".format(OUTPUT, f), im)
