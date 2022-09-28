import pandas as pd
import numpy as np
import os
import cv2

# file paths
DIR = "/mnt/perkinsdata/tongue_STOmics/discovery/gems"
image = os.path.join(DIR, "FP200000495BR_E5_regist.tif")
gems = {"introns": os.path.join(DIR, "tongue-5_introns.gem"),
        "full": os.path.join(DIR, "tongue-5_all.gem")}
temp = os.path.join(DIR, "temp.csv")
out = os.path.join(DIR, "plots")

# create output folder
if not os.path.isdir(out):
    os.mkdir(out) 

# get im
im = cv2.imread(image)
im[:,:,0:2] = 0

# get gem data
def get_gem(path):
    dat = pd.read_csv(path, skiprows = 6, sep = "\t")
    dat = dat.groupby(["x", "y"]).sum().MIDCount
    # likely better way to get back x and y
    dat.to_csv(temp)
    dat = pd.read_csv(temp)
    os.remove(temp)
    return dat

print("reading gems...")
intron_dat = get_gem(gems["introns"])
full = get_gem(gems["full"])

# creating images
def plot_data(dict_data, filename, base_im = None):
    if base_im is None:
        base_im = np.zeros_like(im)
    for (i, dat) in dict_data.items():
        if dat is None: # empty channel
            base_im[:, :, i] = 0
        else:
            base_im[dat.y, dat.x, i] = 255
    cv2.imwrite(os.path.join(out, filename), base_im)

print("plotting...")

if False:
    # plot unthresholded introns
    plot_data({0: intron_dat, 1: intron_dat, 2: intron_dat}, "unthresholded_introns.png")

    plot_data({0: full, 1: full, 2: full}, "full.png")

    # plot thresholded introns
    intron_dat_t = intron_dat[intron_dat.MIDCount > 2]
    plot_data({0: intron_dat_t, 1: intron_dat_t, 2: intron_dat_t}, "thresholded_introns.png")

    # plot non- and thresholded introns
    intron_dat_low = intron_dat[intron_dat.MIDCount <= 2] # different colours of more than 2 and less than 3
    plot_data({0: intron_dat[intron_dat.MIDCount <= 1], 1: intron_dat[intron_dat.MIDCount <= 2], 2: intron_dat[intron_dat.MIDCount > 2]}, "introns.png")

    # plot non- and thresholded introns against image
    plot_data({0: intron_dat_low, 1: intron_dat_t}, "introns_im.png", im)
    plot_data({0: None, 1: intron_dat_t}, "introns_t_im.png", im)

intron_dat_t = intron_dat[intron_dat.MIDCount > 2]

# FIG 1
# compare introns (thresholded and not) to full GEM
plot_data({0: full, 1: None, 2: intron_dat_t}, "gems_t.png")
plot_data({0: full, 1: None, 2: intron_dat}, "gems.png")

# FIG 2
# compare introns (thresholded and not) to image
plot_data({0: None, 1: intron_dat}, "introns_im.png", im)
plot_data({0: None, 1: intron_dat_t}, "introns_t_im.png", im) # FIG 3