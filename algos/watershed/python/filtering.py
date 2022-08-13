import numpy as np
import cv2 as cv

from image_ops import save_multiple_images, plot_histogram, save_labels

"""
    If the given pixel is connected to a label, then returns True, else False
"""
def check_is_bordering_label(markers, x, y, max_x, max_y):
    # positions to check
    positions = []
    if x > 0: positions.append((x-1,y))
    if x < max_x: positions.append((x+1,y))
    if y > 0: positions.append((x,y-1))
    if y < max_y: positions.append((x,y+1))

    for pos in positions: 
        try:
            if markers[pos[0], pos[1]] > 1: return True
        except:
            print(max_x, max_y, x, y)

    return False

"""
    For any border pixels not connecting to a label, sets to background
"""
def check_all_borders_connected(markers):
    borders = np.column_stack(np.where(markers == -1))

    max_x = markers.shape[0] - 1
    max_y = markers.shape[1] - 1

    # loop through all border pixels, checking if connects to label
    for x,y in borders:
        if not check_is_bordering_label(markers, x, y, max_x, max_y):
            markers[x,y] = 1 

    return markers

"""
    Ensures markers are not too small or large
    Limits inclusive
"""
def filter_by_size(markers, min_limit = 25, max_limit = None, im = None, filename = None):
    unravelled = markers.ravel().astype(np.int32)

    # ignore borders (-1) and background (1)
    unravelled = unravelled[unravelled >= 1]

    # bin
    counts = np.bincount(unravelled)

    # set max_limit to upper conf bound
    if not max_limit:
        max_limit = np.mean(counts[2:]) + np.std(counts[2:]) * 1.96 # 95%
    
    # remove labels outside of bounds
    for i in np.where((counts <= min_limit) | (counts >= max_limit))[0]:
        markers[markers == i,] = 1

    markers = check_all_borders_connected(markers)

    if filename is not None and im is not None:
        save_labels(im, markers, filename)

    return markers

"""
    Ensures markers are not too small or large
    Limits inclusive
"""
def filter_by_brightness(im, markers, min_limit = None, max_limit = 255, filename = False):
    unravelled = im.ravel().astype(np.uint8)
    
    labels = np.int64(range(2, np.max(markers).astype(np.int64)))
    print("Num labels", labels)
    #means = np.zeros_like(labels)
    #medians = np.zeros_like(labels)
    #stds = np.zeros_like(labels)

    mean = np.mean(im[np.where(markers > 1)])
    std = np.std(im[np.where(markers > 1)])
    
    if min_limit is None:
        min_limit = mean - std * 1.96
    
    print("(%f %f)" % (min_limit, max_limit))

    # for each label
    for lab in labels:    
        # calc mean brightness
        m = np.mean(im[np.where(markers == lab)])
        if m < min_limit or m > max_limit:
            print("Removing", lab)
            markers[markers == lab] = 1
        elif lab % 200 == 0:
            print(lab)

    if filename is not None:
        save_labels(im, markers, filename)

    return markers