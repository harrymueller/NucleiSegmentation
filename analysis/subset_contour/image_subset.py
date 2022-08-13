import numpy as np
import cv2
import sys

if len(sys.argv) < 3:
    print("py image_subset.py <filename> <num bins> [output dir]")
    print("  <filename>\ttif fake image")
    print("  <num bins>\tNumber of bins to include from boundary")
    print("  [output dir]\tDirectory to save image and bin coords to - optional")
    exit()

filepath = sys.argv[1]
filename = filepath.split("/")[-1]
num_bins = int(sys.argv[2])
if len(sys.argv) == 4:
    output_dir = sys.argv[3]

# read image, threshold, and find countours
im = cv2.imread(filepath, cv2.COLOR_BGR2GRAY)
ret, thresh = cv2.threshold(im, 0, 255, 0)
thresh = thresh.astype(np.uint8)
contours, hierarchy = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)

# create mask and draw contours onto mask
mask = np.zeros_like(im)
cv2.drawContours(mask, contours, -1, 255, num_bins*2+1)

# create clean image
out = np.zeros_like(im)
# where contour is drawn, add real values
out[mask == 255] = im[mask==255] 

if len(sys.argv) == 4:
    # write image
    cv2.imwrite(("%s/%s_subset%d.%s" % (output_dir, filename.split(".")[0], num_bins, filename.split(".")[-1])).replace("UMIGreyScaleFakeIMG_", ""), out)

    x = np.where(out != 0)
    x = np.column_stack(x)

    f = open("%s/%s_subset%d.tsv" % (output_dir, filename.split(".")[0], num_bins), "w")
    f.write("y_coord\tx_coord\n")
    for i in x:
        f.write("%s\t%s\n" % (i[0], i[1]))
    f.close()
else:
    # show image
    cv2.imshow("Subset of %s with %d bins" % (filename, num_bins), out)

    # Wait indefinitely until you push a key.  Once you do, close the windows
    cv2.waitKey(0)
    cv2.destroyAllWindows()
