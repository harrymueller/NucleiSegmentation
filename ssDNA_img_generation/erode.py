import cv2
import numpy as np
import sys

if __name__ == "__main__":
	filename = sys.argv[1]
	im = cv2.imread(filename, 0)

	kernel = np.ones((5,5), np.uint8)

	#img_erosion = cv2.erode(im, kernel, iterations=1)
	#img_dilation = cv2.dilate(im, kernel, iterations=1)
 
	#cv2.imwrite(filename.replace(".png", "_eroded.png"), img_erosion)
	#cv2.imwrite(filename.replace(".png", "_dilated.png"), img_dilation)
	#gim = im
	# for i in range(5):
	# 	#im = cv2.filter2D(im, -1, kernel)
	im = cv2.GaussianBlur(im, (5,5), 0)
	# im = im.astype(np.uint16) + gim
	# im[im > 255] = 255
	# im = im.astype(np.uint8)

	# blurred_img = cv2.GaussianBlur(im, (7, 7), 0)
	# mask = np.zeros(im.shape, np.uint8)

	# #gray = cv2.cvtColor(im, cv2.COLOR_BGR2GRAY)
	# #thresh = cv2.threshold(im, 60, 255, cv2.THRESH_BINARY)[2]
	# contours, hierarchy = cv2.findContours(im, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

	# cv2.drawContours(mask, contours, -1, 255,5)
	# output = np.where(mask==255, blurred_img, im)

	cv2.imwrite(filename.replace(".png", "_dilated.png"), im)