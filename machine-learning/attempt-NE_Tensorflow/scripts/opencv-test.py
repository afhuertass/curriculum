
import numpy as np

import cv2

from matplotlib import pyplot as plt

im = cv2.imread('../data/slither2.jpg', 0 )

#plt.imshow(im),plt.show()
orb = cv2.ORB_create()

kp = orb.detect(im , None)

kp , des = orb.compute(im , kp )

dummy = np.zeros( (1,1) )
img2 = cv2.drawKeypoints(im,kp,color=(0,255,0), outImage=dummy , flags=0  )
plt.imshow(img2),plt.show()

###################33

img1 = cv2.imread('../data/slither.png' , 0)
img2 = cv2.imread('../data/spark.png' , 0 )

orb = cv2.ORB_create()

kp1, des1 = orb.detectAndCompute(img1,None)
kp2, des2 = orb.detectAndCompute(img2,None)

bf = cv2.BFMatcher(cv2.NORM_HAMMING, crossCheck=True)

matches = bf.match( des1 , des2 )

matches = sorted(matches, key = lambda x:x.distance)

img3 = cv2.drawMatches(img1,kp1,img2,kp2,matches[:10],outImg=dummy  ,flags=2)

plt.imshow(img1),plt.show()
