# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>
###########################
## Lsd project 4 part 1
## Zhengjian Song
###########################
## please check the analysis in the markdown below 
import DAL
import numpy as np

import time
import string
import re
import fractions
import math
import scipy
from numpy import arange,array,ones,linalg
from pylab import plot,show
from __future__ import division
from IPython.parallel import Client
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt

lightcurves = DAL.create('lightcurves')
s = lightcurves.subsets()
'''
dat = []
for i in lightcurves.iter(s[101]):
    name = i['id']
    #print type(name),name == 3096237
    if name == 3096237:
        break
    lc = i['data']
    time = lc[:int(len(lc)/2)]  #modified - first half time
    flux = lc[int(len(lc)/2):]  #modified - first half flux
    dat.append( (name, time, flux) ) #modified
'''
print len(s)
print "files:"
print s

# <codecell>

# kernel functions
def boxcar_kernel(x):
    ret = 0
    if( x <= 1 and x>=-1):
        ret = 1
    return ret
def gaussian_kernel(x):
    import math
    return math.exp(-0.5*x*x)/math.sqrt(2*math.pi)
def epan_kernel(x):
    ret = 0
    if( x <= 1 and x >= -1):
        ret = 0.75*(1-x*x)
    return ret
def tric_kernel(x):
    ret = 0
    if( x <= 1 and x >= -1):
        ret =  70*(1- x*x*abs(x))/80
    return ret

# <codecell>

sample_list = [(126,2042713),(245,1255171),(231,485585)]
sample_data = []
IDS = {}
isbreak = False
s_index = 0
for sample in sample_list:
    for si in s:
        if si[:-9] == str(sample[0]):
            break
    print si
    for i in lightcurves.iter(si):
        name = i['id']
        IDS[name] = s_index
        if name == sample[1]:
            lc = i['data']
            time = lc[:int(len(lc)/2)]  #modified - first half time
            flux = lc[int(len(lc)/2):]  #modified - first half flux
            sample_data.append((name, time, flux))
            break

# <codecell>

#nadaraya-watson kernel estimator
def kernel_estimator(X, Y, h, kernel_type):#, kernel_type):
    import numpy as np
    N = len(X)
    L = np.zeros((N, N))
    #default kernel_type = "b"
    L_diag = np.zeros(N)
    
    if(kernel_type[0] == "g"):
        for i in xrange(N):
            #print i
            for j in xrange(i,N):
                L[i,j] = gaussian_kernel((X[i] - X[j])/h)
                L[j,i] = L[i,j]
                if L[i,j] <= 0.01:
                    break
            #scaling
            L[i] = L[i] / L[i].sum()
            L_diag[i] = L[i,i]
        Y_hat = np.dot(L, Y)
        
    elif(kernel_type[0] == "e"):
        for i in xrange(N):
            #print i
            for j in xrange(i,N):
                L[i,j] = epan_kernel((X[i] - X[j])/h)
                L[j,i] = L[i,j]
                if L[i,j] <= 0.01:
                    break
            #scaling
            L[i] = L[i] / L[i].sum()
            L_diag[i] = L[i,i]
        Y_hat = np.dot(L, Y)
        
    elif(kernel_type[0] == "t"):
        for i in xrange(N):
            #print i
            for j in xrange(i,N):
                L[i,j] = tric_kernel((X[i] - X[j])/h)
                L[j,i] = L[i,j]
                if L[i,j] <= 0.01:
                    break
            #scaling
            L[i] = L[i] / L[i].sum()
            L_diag[i] = L[i,i]
        Y_hat = np.dot(L, Y)
        
    else: # boxcar kernel_as default
        #print "b"
        for i in xrange(N):
            #print i
            for j in xrange(i,N):
                L[i,j] = boxcar_kernel((X[i] - X[j])/h)
                L[j,i] = L[i,j]
                if L[i,j] <= 0.01:
                    break
            #scaling
            L[i] = L[i] / L[i].sum()
            L_diag[i] = L[i,i]
        Y_hat = np.dot(L, Y)
    return (Y_hat, L_diag) 

# <codecell>

import time
kernel_list = ["boxcar", "epanechnikov", "tricube", "guassian"]

sample = sample_data[0]
print "-"*10
print "id:", sample[0]
X = np.array(sample[1], dtype = 'float32')
Y = np.array(sample[2], dtype = 'float32')
X = X[np.logical_not(np.isnan(Y))]
Y = Y[np.logical_not(np.isnan(Y))]
t1 = time.clock()
res = kernel_estimator(X, Y, 0.1, "b")
Y_hat = res[0]
#print res[1]
print "time spent ", time.clock() - t1
plt.plot(X, Y)
plt.plot(X, Y_hat)
title_string = "id:" + str(sample[0]) +", kernel regression,"+ kernel_list[0] +" , bandwidth="+str(0.1)
plt.title(title_string)
plt.show()
print "-"*10

t1 = time.clock()
res = kernel_estimator(X, Y, 0.5, "b")
Y_hat = res[0]
#print res[1]
print "time spent ", time.clock() - t1
plt.plot(X, Y)
plt.plot(X, Y_hat)
title_string = "id:" + str(sample[0]) +", kernel regression,"+ kernel_list[0] +" , bandwidth="+str(0.5)
plt.title(title_string)
plt.show()
print "-"*10

t1 = time.clock()
res = kernel_estimator(X, Y, 1.0, "b")
Y_hat = res[0]
#print res[1]
print "time spent ", time.clock() - t1
plt.plot(X, Y)
plt.plot(X, Y_hat)
title_string = "id:" + str(sample[0]) +", kernel regression,"+ kernel_list[0] +" , bandwidth="+str(1.0)
plt.title(title_string)
plt.show()
print "-"*10

t1 = time.clock()
res = kernel_estimator(X, Y, 1.0, "e")
Y_hat = res[0]
#print res[1]
print "time spent ", time.clock() - t1
plt.plot(X, Y)
plt.plot(X, Y_hat)
title_string = "id:" + str(sample[0]) +", kernel regression,"+ kernel_list[1] +" , bandwidth="+str(1.0)
plt.title(title_string)
plt.show()
print "-"*10

t1 = time.clock()
res = kernel_estimator(X, Y, 1.0, "t")
Y_hat = res[0]
#print res[1]
print "time spent ", time.clock() - t1
plt.plot(X, Y)
plt.plot(X, Y_hat)
title_string = "id:" + str(sample[0]) +", kernel regression,"+ kernel_list[2] +" , bandwidth="+str(1.0)
plt.title(title_string)
plt.show()
print "-"*10


t1 = time.clock()
res = kernel_estimator(X, Y, 1.0, "g")
Y_hat = res[0]
#print res[1]
print "time spent ", time.clock() - t1
plt.plot(X, Y)
plt.plot(X, Y_hat)
title_string = "id:" + str(sample[0]) +", kernel regression,"+ kernel_list[3] +" , bandwidth="+str(1.0)
plt.title(title_string)
plt.show()
print "-"*10

# <codecell>

import time
kernel_list = ["boxcar", "epanechnikov", "tricube", "guassian"]

sample = sample_data[1]
print "-"*10
print "id:", sample[0]
X = np.array(sample[1], dtype = 'float32')
Y = np.array(sample[2], dtype = 'float32')
X = X[np.logical_not(np.isnan(Y))]
Y = Y[np.logical_not(np.isnan(Y))]
t1 = time.clock()
res = kernel_estimator(X, Y, 0.1, "b")
Y_hat = res[0]
#print res[1]
print "time spent ", time.clock() - t1
plt.plot(X, Y)
plt.plot(X, Y_hat)
title_string = "id:" + str(sample[0]) +", kernel regression,"+ kernel_list[0] +" , bandwidth="+str(0.1)
plt.title(title_string)
plt.show()
print "-"*10

t1 = time.clock()
res = kernel_estimator(X, Y, 0.5, "b")
Y_hat = res[0]
#print res[1]
print "time spent ", time.clock() - t1
plt.plot(X, Y)
plt.plot(X, Y_hat)
title_string = "id:" + str(sample[0]) +", kernel regression,"+ kernel_list[0] +" , bandwidth="+str(0.5)
plt.title(title_string)
plt.show()
print "-"*10

t1 = time.clock()
res = kernel_estimator(X, Y, 1.0, "b")
Y_hat = res[0]
#print res[1]
print "time spent ", time.clock() - t1
plt.plot(X, Y)
plt.plot(X, Y_hat)
title_string = "id:" + str(sample[0]) +", kernel regression,"+ kernel_list[0] +" , bandwidth="+str(1.0)
plt.title(title_string)
plt.show()
print "-"*10

t1 = time.clock()
res = kernel_estimator(X, Y, 1.0, "e")
Y_hat = res[0]
#print res[1]
print "time spent ", time.clock() - t1
plt.plot(X, Y)
plt.plot(X, Y_hat)
title_string = "id:" + str(sample[0]) +", kernel regression,"+ kernel_list[1] +" , bandwidth="+str(1.0)
plt.title(title_string)
plt.show()
print "-"*10

t1 = time.clock()
res = kernel_estimator(X, Y, 1.0, "t")
Y_hat = res[0]
#print res[1]
print "time spent ", time.clock() - t1
plt.plot(X, Y)
plt.plot(X, Y_hat)
title_string = "id:" + str(sample[0]) +", kernel regression,"+ kernel_list[2] +" , bandwidth="+str(1.0)
plt.title(title_string)
plt.show()
print "-"*10


t1 = time.clock()
res = kernel_estimator(X, Y, 1.0, "g")
Y_hat = res[0]
#print res[1]
print "time spent ", time.clock() - t1
plt.plot(X, Y)
plt.plot(X, Y_hat)
title_string = "id:" + str(sample[0]) +", kernel regression,"+ kernel_list[3] +" , bandwidth="+str(1.0)
plt.title(title_string)
plt.show()
print "-"*10

# <codecell>

import time
kernel_list = ["boxcar", "epanechnikov", "tricube", "guassian"]

sample = sample_data[2]
print "-"*10
print "id:", sample[0]
X = np.array(sample[1], dtype = 'float32')
Y = np.array(sample[2], dtype = 'float32')
X = X[np.logical_not(np.isnan(Y))]
Y = Y[np.logical_not(np.isnan(Y))]
t1 = time.clock()
res = kernel_estimator(X, Y, 0.1, "b")
Y_hat = res[0]
#print res[1]
print "time spent ", time.clock() - t1
plt.plot(X, Y)
plt.plot(X, Y_hat)
title_string = "id:" + str(sample[0]) +", kernel regression,"+ kernel_list[0] +" , bandwidth="+str(0.1)
plt.title(title_string)
plt.show()
print "-"*10

t1 = time.clock()
res = kernel_estimator(X, Y, 0.5, "b")
Y_hat = res[0]
#print res[1]
print "time spent ", time.clock() - t1
plt.plot(X, Y)
plt.plot(X, Y_hat)
title_string = "id:" + str(sample[0]) +", kernel regression,"+ kernel_list[0] +" , bandwidth="+str(0.5)
plt.title(title_string)
plt.show()
print "-"*10

t1 = time.clock()
res = kernel_estimator(X, Y, 1.0, "b")
Y_hat = res[0]
#print res[1]
print "time spent ", time.clock() - t1
plt.plot(X, Y)
plt.plot(X, Y_hat)
title_string = "id:" + str(sample[0]) +", kernel regression,"+ kernel_list[0] +" , bandwidth="+str(1.0)
plt.title(title_string)
plt.show()
print "-"*10

t1 = time.clock()
res = kernel_estimator(X, Y, 1.0, "e")
Y_hat = res[0]
#print res[1]
print "time spent ", time.clock() - t1
plt.plot(X, Y)
plt.plot(X, Y_hat)
title_string = "id:" + str(sample[0]) +", kernel regression,"+ kernel_list[1] +" , bandwidth="+str(1.0)
plt.title(title_string)
plt.show()
print "-"*10

t1 = time.clock()
res = kernel_estimator(X, Y, 1.0, "t")
Y_hat = res[0]
#print res[1]
print "time spent ", time.clock() - t1
plt.plot(X, Y)
plt.plot(X, Y_hat)
title_string = "id:" + str(sample[0]) +", kernel regression,"+ kernel_list[2] +" , bandwidth="+str(1.0)
plt.title(title_string)
plt.show()
print "-"*10


t1 = time.clock()
res = kernel_estimator(X, Y, 1.0, "g")
Y_hat = res[0]
#print res[1]
print "time spent ", time.clock() - t1
plt.plot(X, Y)
plt.plot(X, Y_hat)
title_string = "id:" + str(sample[0]) +", kernel regression,"+ kernel_list[3] +" , bandwidth="+str(1.0)
plt.title(title_string)
plt.show()

print "-"*10

# <codecell>


# <markdowncell>

# We chose 3 light curves [(126,2042713),(245,1255171),(231,485585)] and fitted three curves using different kernels (["boxcar", "epanechnikov", "tricube", "guassian"]) and  bandwidths [0.1 to 1.0]. As we did not have much time, we only plotted a couple of the more interesting plots.
# From the couple of plots above, we see that the smaller the bandwidth, the spikier and wigglier the fitted curve. As for the effect of the kernels, we see that in general, for a fixed bandwidth, the boxcar kernel produces curves which are more wiggly compared to the epanechnikov or gaussian curves. Nevertheless, the choice of kernel is not too important since estimates obtained by using different kernels are numerically very similar. Theoretically, the risk should be rather insensitive to the choice of kernels (see more details in Scott (1992)).
# NOTE: In this part, we fitted the curves together with the "outliers" so the residuals are not that great. We will improve this in part 2

# <codecell>


# <codecell>


# <codecell>


# <codecell>


