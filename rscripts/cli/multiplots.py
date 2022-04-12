import pandas as pd
import numpy as np
import cv2 as cv

# params
NAME = "log_csv"
#NAME = "temp"
ID = "tongue-5_bin10"

OUTPUT_DIR = "/mnt/data/R_analysis/count_feature_plots/gene_expression_plots/" + NAME + "/" + ID
INPUT = OUTPUT_DIR + "/expr.csv"

OFFSET = 100

# read csv
df = pd.read_csv(INPUT)

# divid cols [3,4,5] by max
df.iloc[:,2:5] = df.iloc[:,2:5] / np.max(df.iloc[:,2:5])

# get max x and y
y,x = np.max(df.iloc[:,0:2])

# create blank image
im = np.ones((y+OFFSET,x,3), np.uint8)*255

# add colours
def set_rgb(arr):
    arr = np.array(arr)
    bgr = np.array(arr[2:5] * 255, np.uint8)
    
    im[int(arr[0])-1+OFFSET][int(arr[1])-1] = bgr

df.apply(set_rgb, axis=1)

# add labels
label = "red: " + df.columns[4] + " - " + \
        "green: " + df.columns[3] + " - " + \
        "blue: " + df.columns[2]

cv.putText(im, label, 
    (20,50), 
    cv.FONT_HERSHEY_SIMPLEX, 
    1,
    (0,0,0),
    2,
    2)


# save image to file
cv.imwrite(OUTPUT_DIR + "/image.png", im)