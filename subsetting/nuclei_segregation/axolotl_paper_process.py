import numpy as np
import cv2 as cv
import sys
from matplotlib import pyplot as plt

TESTING = True
DEBUGGING = True # whether to save lots of images throughout process

if len(sys.argv) != 3 and not TESTING:
    print("py test.py <input> <output>")
    print("  <input>\t\tPath to ssDNA stain image")
    print("  <output folder>\tPath to output folder")
    exit()
elif not TESTING:
    FILEPATH = sys.argv[1]
    OUTPUT = sys.argv[2]
else:
    #FILEPATH = "/mnt/perkinsdata/tongue_STOmics/original/image- ssDNA and H&E/Mouse-tongue-5-Library Construction/FP200000495BR_E5.tif"
    FILEPATH = "/mnt/perkinsdata/tongue_STOmics/discovery/nuclei_segregation/testing/FP200000495BR_E5_6000_3000.tif"
    OUTPUT = "/mnt/perkinsdata/tongue_STOmics/discovery/nuclei_segregation/testing"

# open a window to show the passed image
def show_image(im, title = ""):
    # show image
    cv.imshow(title, im)

    # Wait indefinitely until you push a key.  Once you do, close the windows
    cv.waitKey(0)
    cv.destroyAllWindows()

# write image to file
def write_image(im, filename):
    cv.imwrite("%s/%s" % (OUTPUT, filename), im)

# given a mask, create a duplicate array and add pixels values based on mask
def mask_image(im, mask):
    out = np.zeros_like(im)
    out[mask] = im[mask]
    return out

"""
    Saves a set of images together with titles

    Settings specifically for 6000x3000px image
"""
def save_multiple_images(images, titles, filename):
    grey = (192)

    # create vertical space
    v_space = np.zeros((280, images[0].shape[1] + 100, 1), np.uint8)
    v_space[:] = grey # set col to grey

    for i in range(len(images)):
        # add thin white border && add thicker grey border
        images[i] = cv.copyMakeBorder(images[i], 10, 10, 10, 10, cv.BORDER_CONSTANT, value = (255))
        images[i] = cv.copyMakeBorder(images[i], 20, 20, 40, 40, cv.BORDER_CONSTANT, value = grey)

        # add vertical grey space above image
        images[i] = cv.vconcat([v_space, images[i]]) 

        # add text
        cv.putText(images[i], titles[i], (80,220), cv.FONT_HERSHEY_SIMPLEX, 8, (0), 16)

    # concat all separate images together, then save to file
    images[i] = cv.vconcat(images) if images[0].shape[1] > images[0].shape[0] else cv.hconcat(images)
    write_image(images[i], filename) 

"""
    Apply a global threshold to the passed image

    @param val: Greyscale level to threshold at

    @return threshold array containing two values - 255 = incl && 0 = excl
"""
def global_thresholding(im, val = 26):
    # apply threshold
    ret, thresh = cv.threshold(im, val, 255, 0)

    if DEBUGGING:
        # save different images to debug the threshold
        save_multiple_images([thresh, mask_image(im, thresh == 255), mask_image(im, thresh == 0)],
                             ["Thresholded (Binary)", "Included Pixels", "Excluded Pixels"],
                             "global_thresholding.png")

    return thresh

"""
    Apply a gaussian-weighted local-threshold to the passed image

    @param blockSize - size of neighbourhood
    @param C - offset 

    @return threshold array containing two values - 255 = incl && 0 = excl
"""
def gaussian_thresholding(im, blockSize = 41, C = 0.03):
    thresh = cv.adaptiveThreshold(im, 255, cv.ADAPTIVE_THRESH_GAUSSIAN_C, cv.THRESH_BINARY, blockSize, C)

    if DEBUGGING:
        # save different images to debug the threshold
        save_multiple_images([thresh, mask_image(im, thresh == 255), mask_image(im, thresh == 0)],
                             ["Thresholded (Binary)", "Included Pixels", "Excluded Pixels"],
                             "gaussian_thresholding.png")
    
    return thresh


def main():
    # read image
    im = cv.imread(FILEPATH, cv.IMREAD_GRAYSCALE)

    # apply thresholds
    global_thresh = global_thresholding(im, val = 50)
    gaussian_tresh = gaussian_thresholding(mask_image(im, global_thresh == 255), blockSize = 101, C = 0)

    # apply Euclidean dist transformation

    # apply watershed algo
    

if __name__ == "__main__":
    main()