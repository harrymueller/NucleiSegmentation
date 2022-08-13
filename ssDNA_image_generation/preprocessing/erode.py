import cv2
import numpy as np
import sys
import os.path

DIR = "5x5"
INPUT = os.path.join(DIR, "ssDNA_stains_raw")
OUTPUT = os.path.join(DIR, "ssDNA_stains_fuzzy")

NOISE_KERNEL = (25,25)
NOISE_SIGMA = 5

BLUR_KERNEL = (5,5)
BLUR_SIGMA = 1

# if output exists, clean, otherwise create output folder
if os.path.isdir(OUTPUT):
	for f in os.listdir(OUTPUT):
		os.remove("{0}/{1}".format(OUTPUT, f))
else: os.mkdir(OUTPUT)

# loop through all images
for f in os.listdir(INPUT):
	im = cv2.imread("{0}/{1}".format(INPUT, f), 0)
	
	# apply blur
	noise = cv2.GaussianBlur(im, NOISE_KERNEL, NOISE_SIGMA, NOISE_SIGMA, 0)
	im = cv2.GaussianBlur(im, BLUR_KERNEL, BLUR_SIGMA, BLUR_SIGMA, 0)
	
	mask = noise > im # use noise image where noise > blur
	im[mask] = noise[mask]

	cv2.imwrite("{0}/{1}".format(OUTPUT, f), im)