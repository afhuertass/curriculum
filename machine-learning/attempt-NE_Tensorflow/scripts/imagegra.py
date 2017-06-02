from __future__ import division

import pyscreenshot as ImageGrab
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image 

import time


prop = 0.5
n = 100
t0 = time.clock()
for i in range(n):
    print(i)
    
    #im = ImageGrab.grab()
    im = ImageGrab.grab(bbox=(10,10,810,810))
    nx , ny = im.size

    im2 = im.resize( (int(0.3*nx) , int(0.3*ny) ) ,  Image.ANTIALIAS )

    #plt.imshow( np.array(im2) )

tf = time.clock()

tt = tf - t0
print("Tiempo total:" + str(tt))

avg = n/tt

print("Capturas por segundo:" + str(avg))


im = ImageGrab.grab(bbox=(10,10,810,810))
nx , ny = im.size

im2 = im.resize( (int(0.3*nx) , int(0.3*ny) ) ,  Image.ANTIALIAS )

plt.imshow( np.array(im2) )
plt.show()
