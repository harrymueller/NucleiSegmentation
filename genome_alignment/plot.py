import pandas as pd
import os
import cv2

# file paths
dir = "/mnt/perkinsdata/tongue_STOmics/discovery/gems"
image = os.path.join(dir, "FP200000495BR_E5_regist.tif")
gem = os.path.join(dir, "tongue-5_introns.gem")
temp = os.path.join(dir, "temp.csv")
out = os.path.join(dir, "introns.png")

# get im
im = cv2.imread(image)
im[:,:,0:2] = 0

# get ge data
dat = pd.read_csv(gem, skiprows = 6, sep = "\t")
dat = dat.groupby(["x", "y"]).sum().MIDCount

# better way to get back x and y
dat.to_csv(temp)
dat = pd.read_csv(temp)
os.remove(temp)

# threshold
dat = dat[dat.MIDCount > 2]
im[dat.y, dat.x, 1] = 255

# save
cv2.imwrite(out, im)