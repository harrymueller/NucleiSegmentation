import numpy as np

# minPos = (sizes - sphereSizes) / -2.0f + borderSizes;
# maxPos = (sizes - sphereSizes) /  2.0f - borderSizes;
xsize = 192
maxscale = 1
bordersize = 2
# xmin, xmax
x = [-93.5, 93.5]
# rmin, rmax
r = x
# proportions
proportions = [0.25, 0.50]

rmax = 200
pos1 = (r[0], x[0]) # rmin, xmin
pos2 = (proportions[0] * (r[1] - r[0]) + r[0], proportions[1] * (x[1] - x[0]) + x[0]) # a,b
pos3 = (r[1], r[1]) # rmax, xmax

print(pos1)
print(pos2)
print(pos3)

X = np.matrix([[pos1[0]**2, pos1[0], 1],
               [pos2[0]**2, pos2[0], 1],
               [pos3[0]**2, pos3[0], 1]])
y = np.matrix([[pos1[1]], [pos2[1]], [pos3[1]]])

sol = np.linalg.inv(X) * y
print(sol)
print(X * sol)