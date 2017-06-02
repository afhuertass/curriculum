# -*- coding: utf-8 -*-
#import scipy.ndimage as nd
#filtos:http://docs.scipy.org/doc/scipy/reference/ndimage.html
#plot: http://matplotlib.org/users/image_tutorial.html
#import numpy as np
#import pylab as pl
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
#Lee la imagien------

#image=nd.imread("Sal-3-2.png").astype("double")
img=mpimg.imread('Sal-3-2.png').astype("double")
imgplot = plt.imshow(img)
plt.show(img)
lum_img = img[:,:,0]
plt.show(lum_img)

#Cambia el color-----

#R = image[:,:,0]  
#G = image[:,:,1]  
#B = image[:,:,2]  

#imGris=0.21*R + 0.71*G + 0.07*B  

#Filtro Gaussiano----

#n=10
#l=256
#im_filtrada=nd.gaussian_filter(imGris,sigma=l/(4.*n))

#Binarizar la imagen------

#[M,N]=imGris.shape

#for i in range (M):
#    for j in range(N):
#        if(im_filtrada[i,j]<120):
#            im_filtrada[i,j]=255
#        else:
#            im_filtrada[i,j]=1

#mask= im_filtrada > im_filtrada.mean()

#label_im,nb_labels=nd.label(mask)

#El area-----------------

#suma=0.0
#for i in range (M):
#    for j in range(N):
#        if(label_im[i,j]>1):
#            suma+=1.0

#pixeles=M*N
#porcentaje=(suma/pixeles)*100

#pl.figure()
#pl.imshow(label_im,cmap="gray",origin="upper")


#print "El numero total de objetos es {}".format(nb_labels-1)
#print "El area total de los objetos es {}".format(suma)
#print "El porcentaje de ocupacion de la imagen es {:4.2f} % ".format(porcentaje)

#pl.show()
