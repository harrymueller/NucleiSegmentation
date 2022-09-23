import numpy as np
import pandas as pd
import cv2, os

INPUT = "/mnt/perkinsdata/tongue_STOmics/discovery/gems/tongue-5_all.gem"
OUTPUT = "/mnt/perkinsdata/tongue_STOmics/discovery/gems"
IM_OUTPUT = os.path.join(OUTPUT, "tongue-5.png")
GEM_OUTPUT = os.path.join(OUTPUT, "tongue-5_filtered.gem")

# positions
POS0 = np.array((6950, 15700)) #y,x
SIZE = np.array((15600, 8500))
POS1 = POS0+SIZE

# load data
data = pd.read_csv(INPUT, sep = "\t", skiprows = 6)

# subset data
new = data[np.array(data["x"] > POS0[1]) & np.array(data["x"] < POS1[1])]
new = new[np.array(new["y"] > POS0[0]) & np.array(new["y"] < POS1[0])]

# move back to the origin
new["x"] = new["x"] - POS0[1]
new["y"] = new["y"] - POS0[0]

# transpose x and y values
temp = new["y"]
new["y"] = new["x"]
new["x"] = temp

# create an image
print("Writing image")
maxes = new.max()
im = np.zeros(maxes[1:3]+100, np.uint8)
im[new["x"], new["y"]] += new["MIDCount"]
im = im / np.quantile(im, 0.999) * 255

cv2.imwrite(IM_OUTPUT, im)

# sort
print("sorting")
new = new.sort_values(["geneID","x","y"])

# save gem output
print("saving")
new.to_csv(GEM_OUTPUT, sep = "\t", index = False)