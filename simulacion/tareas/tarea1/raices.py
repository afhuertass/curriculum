## programa para encontrar las raices de la funcion f(lambda)

import numpy as np

data = np.recfromcsv("2bessel.dat");

#data[0][1-max] son los valores de lamba
#data[0][1-max] son los valores de la funcion 
size =  len(data)
print "N Lambda"
raices = 0
for i in xrange(0,size-1):
    
    lam=data[i][0]
    lam_sig=data[i+1][0]
    f = data[i][1]
    f_sig = data[i+1][1]
    if (f*f_sig < 0):
        #la raiz esta entre los dos 
        lam_raiz = (lam+lam_sig)/2;
        print str(raices) + " " +str(lam_raiz)
        raices = raices+1
