import numpy as np
import cv2 as cv
from matplotlib import pyplot as plt

# open a window to show the passed image
def show_image(im, title = ""):
    # show image
    cv.imshow(title, im)

    # Wait indefinitely until you push a key.  Once you do, close the windows
    cv.waitKey(0)
    cv.destroyAllWindows()

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
    # colour definitions
    white = (255, 255, 255)
    grey = (192, 192, 192)
    black = (0, 0, 0)

    # create vertical space
    v_space = np.zeros((280, images[0].shape[1] + 100, 3), np.uint8)
    v_space[:] = grey # set col to grey

    for i in range(len(images)):
        if len(images[i].shape) <= 2: # assumed 8UCV1
            images[i] = cv.cvtColor(images[i], cv.COLOR_GRAY2BGR)

        # add thin white border && add thicker grey border
        images[i] = cv.copyMakeBorder(images[i], 10, 10, 10, 10, cv.BORDER_CONSTANT, value = white)
        images[i] = cv.copyMakeBorder(images[i], 20, 20, 40, 40, cv.BORDER_CONSTANT, value = grey)

        # add vertical grey space above image
        images[i] = cv.vconcat([v_space, images[i]]) 

        # add text
        cv.putText(images[i], titles[i], (80,220), cv.FONT_HERSHEY_SIMPLEX, 8, black, 16)

    # concat all separate images together, then save to file
    if len(images) > 1:
        images[i] = cv.vconcat(images) if images[0].shape[1] > images[0].shape[0] else cv.hconcat(images)
    cv.imwrite(filename, images[i])


"""
    Plot a histogram of values for the given grayscale image
"""
def plot_histogram(im, filename, bins = 256, hist_range = (0, 256)):
    plt.figure()

    plt.hist(im.ravel(), bins, hist_range)
    plt.title('Histogram for gray scale image')
    plt.yscale('log')

    fig = plt.gcf()
    fig.set_size_inches(10, 6)

    plt.savefig(filename, dpi = 300, pad_inches = 0)
