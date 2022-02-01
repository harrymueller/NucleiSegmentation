import numpy as np
import cv2 as cv
import sys
from matplotlib import pyplot as plt

# import local files
from image_ops import *
import thresholding
from watershed import *
from filtering import *

TESTING = True

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
    FILEPATH = "/mnt/perkinsdata/tongue_STOmics/discovery/nuclei_segregation/tongue-5_testing/FP200000495BR_E5.tif"
    OUTPUT = "/mnt/perkinsdata/tongue_STOmics/discovery/nuclei_segregation/tongue-5_testing"

def get_filepath(filename):
    return "%s/%s" % (OUTPUT, filename)

def main():
    print("Reading in image...")
    # read image
    im = cv.imread(FILEPATH, cv.IMREAD_GRAYSCALE)
    plot_histogram(im, get_filepath("histogram.png"))

    # apply thresholds
    print("Applying thresholds...")
    global_thresh = thresholding.global_thresholding(im, val = 100, filename = get_filepath("global_threshold.png"))
    gaussian_thresh = thresholding.gaussian_thresholding(mask_image(im, global_thresh == 255), extra_thresh = global_thresh, blockSize = 41, C = 0.03, filename = get_filepath("gaussian_threshold.png"))

    # determine regions highly likely to be either nuclei || background
    print("Determing sure fg and bg...")
    sure_fg, sure_bg = isolate_fg_bg(gaussian_thresh, filename = get_filepath("nuclei_background_isolation.png"))
    unknown = cv.subtract(sure_bg,sure_fg)

    # get markers -> "seed" positions
    print("Determining markers...")
    markers = get_markers(sure_fg, unknown, get_filepath("markers.png"))

    # apply watershed algo
    print("Applying watershed algorithm...")
    base_image = mask_image(im, gaussian_thresh == 255)
    base_image = cv.cvtColor(base_image, cv.COLOR_GRAY2BGR) # CV_8UC3
    markers = watershed(base_image, markers, im, get_filepath("watershed.png"))

    """
        Produces a border discrepancy in such the outer most layer of pixels is counted as a border
    """

    # save markers
    print("Saving labels as a csv...")
    np.savetxt(get_filepath("markers_pre_filtering.csv"), markers, delimiter = ",")

    print("Filtering labels by size...")
    markers = filter_by_size(markers, min_limit = 36, max_limit = 2**16 - 1, im = im, filename = get_filepath("size_filtering.png"))
    np.savetxt(get_filepath("markers_post_filtering.csv"), markers, delimiter = ",")

    plot_histogram(im[markers > 1], get_filepath("histogram_labels.png"))

    print("Done.")

if __name__ == "__main__":
    main()
