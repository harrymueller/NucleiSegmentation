import numpy as np
import matplotlib.image as mpimg
import os

# consts
INPUT_FILE = "/mnt/perkinsdata/tongue_STOmics/discovery/cellpose/images/tongue-4_auto_seg.npy"
OUTPUT_DIR = "/mnt/perkinsdata/tongue_STOmics/discovery/cellpose/tongue-4"

PATH = os.path.join(OUTPUT_DIR, "02_cellpose")

# make folders
if not os.path.isdir(PATH): os.mkdir(PATH)

# load NPY array
print("Loading NPY object")
x = np.load(INPUT_FILE, allow_pickle = True).item()

# show segment outlines on ssDNA image
# derived from cellpose.plot.outline_view 
print("Segmentation outlines")
im = x["img"]
if len(im.shape)<3: im = np.stack([im]*3,axis=-1)

outY, outX = np.nonzero(x["outlines"])
im[outY, outX] = np.array([255,0,0])
mpimg.imsave(os.path.join(PATH, "cellpose_segmented.png"), im)

# save to masks
print("Saving masks")
masks = x["masks"] + 1
np.savetxt(os.path.join(PATH, "segments.csv"), masks, delimiter=",", fmt = "%d")