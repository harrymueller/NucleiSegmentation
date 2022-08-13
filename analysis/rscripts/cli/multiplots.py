import pandas as pd
import numpy as np
import cv2 as cv

# params
NAME = "new_triple_1"
#NAME = "temp"
ID = "tongue-5_bin1"

APPLY_KERNEL = False
KERNEL_SIZE = 3

OUTPUT_DIR = "/mnt/data/R_analysis/count_feature_plots/gene_expression_plots/" + NAME + "/" + ID
INPUT = OUTPUT_DIR + "/expr.csv"

OFFSET = 100

NUM_PLOTS = 8

# read csv
df = pd.read_csv(INPUT)

max_cols = NUM_PLOTS*3 + 1

# divid cols [3,4,5] by max
df.iloc[:,2:max_cols] = df.iloc[:,2:max_cols] / np.max(df.iloc[:,2:max_cols])

# get max x and y
y,x = np.max(df.iloc[:,0:2])

I_OFFSET = 0

# set rgb func
def set_rgb(arr):
    arr = np.array(arr)

    if len(np.nonzero(arr[(I_OFFSET*3+2):(I_OFFSET*3+5)])[0]) > 0:
        bgr = np.array(np.ceil(arr[(I_OFFSET*3+2):(I_OFFSET*3+5)]) * 255, np.uint8)
        im[int(arr[0])-1+OFFSET][int(arr[1])-1] = bgr
    else:
        im[int(arr[0])-1+OFFSET][int(arr[1])-1] = np.array([240,240,240])

for i in range(NUM_PLOTS):
    print(i)
    # create blank image
    im = np.ones((y+OFFSET,x,3), np.uint8)*255
    I_OFFSET = i
    # add colours
    df.apply(set_rgb, axis=1)

    # apply kernel
    if (APPLY_KERNEL):
        kernel = np.ones((KERNEL_SIZE,KERNEL_SIZE))
        new_im = cv.filter2d(im, ddepth = -1, kernel = kernel)

    # add labels
    label = "red: " + df.columns[i*3+4] + " - " + \
            "green: " + df.columns[i*3+3] + " - " + \
            "blue: " + df.columns[i*3+2]

    cv.putText(im, label, 
        (20,50), 
        cv.FONT_HERSHEY_SIMPLEX, 
        1,
        (0,0,0),
        2,
        2)

    filename = df.columns[i*3+2] + "_" + df.columns[i*3+3] + "_" + df.columns[i*3+4] + ".png"

    # save image to file
    cv.imwrite(OUTPUT_DIR + "/" + filename, im)