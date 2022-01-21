import numpy as np
import cv2 as cv

from image_ops import mask_image, save_multiple_images

"""
    Apply a global threshold to the passed image

    @param val: Greyscale level to threshold at

    @return threshold array containing two values - 255 = incl && 0 = excl
"""
def global_thresholding(im, val = 26, filename = False):
    # apply threshold
    ret, thresh = cv.threshold(im, val, 255, 0)

    if filename:
        # save different images to debug the threshold
        save_multiple_images([thresh, mask_image(im, thresh == 255), mask_image(im, thresh == 0)],
                             ["Thresholded (Binary)", "Included Pixels", "Excluded Pixels"],
                             filename)

    return thresh

"""
    Apply a gaussian-weighted local-threshold to the passed image

    @param blockSize - size of neighbourhood
    @param C - offset 

    @return threshold array containing two values - 255 = incl && 0 = excl
"""
def gaussian_thresholding(im, extra_thresh = None, blockSize = 41, C = 0.03, filename = False):
    thresh = cv.adaptiveThreshold(im, 255, cv.ADAPTIVE_THRESH_GAUSSIAN_C, cv.THRESH_BINARY, blockSize, C)
    if extra_thresh is not None:
        thresh[extra_thresh == 0] = 0

    if filename:
        # save different images to debug the threshold
        save_multiple_images([thresh, mask_image(im, thresh == 255), mask_image(im, thresh == 0)],
                             ["Thresholded (Binary)", "Included Pixels", "Excluded Pixels"],
                             filename)
    
    return thresh