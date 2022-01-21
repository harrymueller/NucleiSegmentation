import numpy as np
import cv2 as cv

import thresholding
from image_ops import save_multiple_images, save_labels

"""
    Isolates fg and bg from the threshold

    Fg isolated via Euclidean dist transformation then gaussian local thresholding
    Bg isolated via dilating threshold 

"""
def isolate_fg_bg(thresh, filename = False):
    dist_transform = cv.distanceTransform(thresh, cv.DIST_L2, 5) # eucl dist transformation
    dist_transform = (dist_transform * (255 / dist_transform.max())).astype(np.uint8) # scales so max = 255
    
    # sure fg calculated by using a gaussian threshold
    sure_fg = thresholding.gaussian_thresholding(dist_transform, extra_thresh = dist_transform, blockSize = 21, C = -30)
    
    # sure bg calculated by dilate the threshold
    kernel = np.ones((3,3), np.uint8)
    sure_bg = cv.dilate(thresh, kernel, iterations=3)
        
    if filename is not False:
        save_multiple_images([dist_transform, sure_fg, sure_bg], ["Euclidean Distance Transformation", "Sure Foreground", "Sure Background"], filename)

    return (sure_fg, sure_bg)

"""
    Given the sure_fg and the unknown region, produces an image containing seed positions
        Label markers

    @param sure_fg, unknown, filename
    @return array of size sure_fg containing seed positions
"""
def get_markers(sure_fg, unknown, filename = False):
    ret, markers = cv.connectedComponents(sure_fg)

    # Makes known background label = 1, not 0
    markers = markers + 1

    # set label of unknown pixels to 0
    markers[unknown==255] = 0

    # if filename given, make a JET colormap, and save to file
    if filename:
        markers_out = markers.astype(np.uint8) # cv.applyColorMap reqs uint8 dtype
        markers_out = cv.applyColorMap(markers_out, cv.COLORMAP_JET)

        # change color of unknown region to light grey
        markers_out[unknown==255] = (192,192,192)
        save_multiple_images([markers_out], ["Markers"], filename)

    return markers

"""
    Apply the watershed algo to im using the seed positions labeled in markers

    @param im: image of type CV8UC3 used for watershed algo
    @param markers: of type CV_32SC1 -> seeds
    @param display_image: image to use to display
"""
def watershed(im, markers, display_image = False, filename = None):
    # watershed
    markers = cv.watershed(im, markers)

    if filename is not None:
        save_labels(display_image if display_image is not None else im, markers, filename)

    return markers