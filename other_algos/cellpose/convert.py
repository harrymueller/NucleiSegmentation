import numpy as np
import matplotlib.image as mpimg
import os

# consts
INPUT_DIR = "/mnt/perkinsdata/tongue_STOmics/benchmarking/25_2k_4/cellpose/pickles"
OUTPUT_DIR = "/mnt/perkinsdata/tongue_STOmics/benchmarking/25_2k_4/cellpose/results"

outputs = [os.path.join(OUTPUT_DIR, "outlines"), 
           os.path.join(OUTPUT_DIR, "segments")]

# make output folder if not made
if not os.path.isdir(OUTPUT_DIR): os.mkdir(OUTPUT_DIR)
if not os.path.isdir(outputs[0]): os.mkdir(outputs[0])
if not os.path.isdir(outputs[1]): os.mkdir(outputs[1])

for f in os.listdir(INPUT_DIR):
    i = f.split("_")[1] # get image number
    print(">>> " + i)
    
    # load NPY array
    print("Loading NPY object")
    x = np.load(os.path.join(INPUT_DIR, f), allow_pickle = True).item()

    # show segment outlines on ssDNA image
    # derived from cellpose.plot.outline_view 
    print("Segmentation outlines")
    im = x["img"]
    if len(im.shape)<3: im = np.stack([im]*3,axis=-1)

    outY, outX = np.nonzero(x["outlines"])
    im[outY, outX] = np.array([255,0,0])
    mpimg.imsave(os.path.join(outputs[0], "cellpose_segmented_%s.png" % (i)), im)

    # save to masks
    print("Saving masks")
    masks = x["masks"] + 1
    np.savetxt(os.path.join(outputs[1], "segments_%s.csv" % (i)), masks, delimiter=",", fmt = "%d")