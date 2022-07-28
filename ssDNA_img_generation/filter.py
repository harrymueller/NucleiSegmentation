import cv2
import numpy as np

im = cv2.imread("truth_0.png")
MINCOUNT = 25

# create cluster ids
segments = np.ones(im.shape[0:2], dtype = np.uint32)
counter = 2
for pos in np.dstack(np.where(np.sum(im, axis = 2) != 0))[0]: # loop through all positions
    # check not done
    if np.sum(im[pos[0], pos[1]]) == 0: continue

    # -15:35
    col = im[pos[0], pos[1]]
    colour_mask = np.sum(im == col, axis = 2) == 3

    # limits size
    if pos[0] > 16: colour_mask[:pos[0]-26,:] = False
    if pos[0] < im.shape[0]-56: colour_mask[pos[0]+55:,:] = False
    if pos[1] > 26: colour_mask[:,:pos[1]-26] = False
    if pos[1] < im.shape[1]-56: colour_mask[:,pos[1]+55:] = False

    if np.sum(colour_mask) >= MINCOUNT:
        segments[colour_mask] = counter
        counter+=1

    im[colour_mask] = 0

np.savetxt("segments.csv", segments, delimiter=",", fmt = "%d")

#cv2.imwrite("temp0.png", segments)

"""
def get_values(i, j):
    vals = np.array([im[i-1, j-1], im[i-1, j], im[i-1, j+1],
                     im[i, j-1],               im[i, j+1],
                     im[i+1, j-1], im[i+1, j], im[i+1, j+1]], dtype = np.int16)
    return vals
cv2.imwrite("truth_0_filtered.png", cv2.cvtColor(im, cv2.COLOR_BGR2GRAY))
exit()

# ensure all colours the same
for pos in np.dstack(np.where(np.sum(im, axis = 2) != 0))[0].astype(np.int16):
    vals = get_values(pos[0], pos[1])
    if pos[1] == 27 and pos[0] == 16:
        print(vals)
        print(np.abs(im[pos[0], pos[1]]*2 - vals))
        print(np.sum(np.abs(im[pos[0], pos[1]]*2 - vals), axis = 1))
    #print(vals)
    diffs = np.sum(np.abs(im[pos[0], pos[1]]*2 - vals), axis = 1)
    i = np.where(diffs<4)[0]
    if i.size > 0:
        im[pos[0], pos[1]] = vals[i[0]]

cv2.imwrite("truth_0_filtered.png", im.astype(np.uint8))
"""