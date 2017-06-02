
import cv2

import numpy as np
import matplotlib.pyplot as plt
from skimage import filters

im = cv2.imread('../data/slither1.png')

img_gray = cv2.cvtColor( im , cv2.COLOR_RGB2GRAY)

ret, th = cv2.threshold( img_gray , 0 , 255 , cv2.THRESH_BINARY+cv2.THRESH_OTSU)


plt.imshow( th , cmap='gray' )
plt.show()

