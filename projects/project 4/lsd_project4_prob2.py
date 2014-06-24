# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

## Lsd project 4 part 2
## Zhengjian Song
## please check the analysis in the markdown below 
import DAL
import numpy as np
from DAL.datasets.checkpoint import Checkpoint
from IPython.parallel import Client

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
rc = Client()
dview = rc[:]

checkpoint = Checkpoint()

print "files:"
print s

# <codecell>


# <codecell>

# we omit the constants since they cancel
def bkernel(x,y,h):
    m = len(x)
    K = np.repeat(x, m).reshape([m,m])
    K = np.absolute(K - K.T) / h
    L = (1.0)*(K<=1)
    L = (L.T / np.sum(L, axis=1)).T
    y_hat=np.dot(L, y.T)
    print 
    return (x,y,y_hat,L)

def gkernel(x,y,h):
    m = len(x)
    K = np.repeat(x, m).reshape([m,m])
    K = np.absolute(K - K.T) / h
    L = np.exp(-0.5 * K**2)
    L = (L.T / np.sum(L, axis=1)).T
    y_hat=np.dot(L, y.T)
    return (x,y,y_hat)

def ekernel(x,y,h):
    m = len(x)
    K = np.repeat(x, m).reshape([m,m])
    K = np.absolute(K - K.T) / h
    L = (1.0 - K**2) * (K<=1)
    L = (L.T / np.sum(L, axis=1)).T
    y_hat=np.dot(L, y.T)
    return (x,y,y_hat)

def tkernel(x,y,h):
    m = len(x)
    K = np.repeat(x, m).reshape([m,m])
    K = np.absolute(K - K.T) / h
    L = (1.0 - K**3) * (K<=1)
    L = (L.T / np.sum(L, axis=1)).T
    y_hat=np.dot(L, y.T)
    return (x,y,y_hat)

## NW estimator
## We chose the boxcar kernel
def nad_wat(X, Y, h):
    N = len(X)
    #m = len(z)
    L_diag = L = np.zeros(N,dtype='float32')
    Yhat = np.zeros(N, dtype='float32') 
    for i in range(N):
        x = (X[i] - X) / h
        # Boxcar
        K = np.array(abs(x) <= 1, dtype='float32') * 1
        L = K / np.sum(K)
        Yhat[i] = np.dot(L, Y)
        L_diag[i] = L[i]
    return (Yhat, L_diag)

def nad_wat_with_L(X, Y, h):
    N = len(X)
    #m = len(z)
    L = np.zeros((N,N),dtype='float32')
    Yhat = Lt =np.zeros(N, dtype='float32') 
    for i in xrange(N):
        x = (X[i] - X) / h
        # Boxcar
        K = np.array(abs(x) <= 1, dtype='float32') * 1
        Lt = K / np.sum(K)
        Yhat[i] = np.dot(Lt, Y)
        L[i] = Lt
    return (Yhat, L)

# <codecell>

# nadaraya-watson estimator, with cross validation
def kernel_regress_cross_validation(X, Y):#, kernel_type):
    import numpy as np
    #default kernel_type = "b"
    N = len(X)
    L = np.zeros((N, N))
    H = np.arange(0.3, 1.1, 0.3) #bandwidth values, here only 0.1 and 0.5, together with 1.0 outside
    Risk_loo = np.zeros(len(H))
    best_hat = np.zeros(N)
    best_risk = 1e15
    best_h = 0.1
    is_too_small = False
    
    for k in xrange(len(H)):
        h = H[k]
        result = nad_wat(X, Y, h)
        Y_hat = result[0]
        L_diag = result[1]
        risk = (Y-Y_hat)/(1-L_diag)
        n_empi_risk = np.sum(risk**2)
        
        Risk_loo[k] = n_empi_risk
        if(best_risk > n_empi_risk):
            best_risk = n_empi_risk
            best_hat = Y_hat
            best_h = H[k]
    '''     
    #plt.plot(Risk_loo)
    #plt.show()
    fig = plt.figure(figsize=(10,6))
    plt.plot(X, Y, '.k')
    plt.plot(X, best_hat, '-b')
    title_string = " bandwidth="+str(best_h)
    plt.title(title_string)
    plt.show()
    '''
    
    return (best_hat, best_risk, best_h)

# <codecell>


# <codecell>

# small demo to show the fit of several curves
'''
from scipy.stats.kde import gaussian_kde
from scipy.stats import norm
from numpy import linspace,hstack
from pylab import plot,show,hist

sample_list = [(189,914193), (126,2042713), (276,6852096), \
(37,3972533), (245,1255171), (211,5004731), (199,4848370), (64,6258272),\
(102,3096237), (313,7412246)]
norm_ones = []
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
            print name
            lc = i['data']
            time = lc[:int(len(lc)/2)]  #modified - first half time
            flux = lc[int(len(lc)/2):]  #modified - first half flux
            sample_data.append((name, time, flux))
            break


import time
for sample in sample_data:
    print "-"*10
    print "id:", sample[0]
    X = np.array(sample[1], dtype = 'float32')
    Y = np.array(sample[2], dtype = 'float32')
    X = X[np.logical_not(np.isnan(Y))]
    Y = Y[np.logical_not(np.isnan(Y))]
    t1 = time.clock()
    h = 0.1
    #res = nad_wat(X, Y, 0.1)
    res = kernel_regress_cross_validation(X,Y)
    Y_hat = res[0]

    print "time spent ", time.clock() - t1
    fig = plt.figure(figsize=(10,6))
    plt.plot(X, Y, '.k')
    plt.plot(X, Y_hat, '-b')
    title_string = "id:" + str(sample[0]) +", kernel regression, boxcar , bandwidth="+str(h)
    plt.title(title_string)
    plt.show()
    
    ###################################
    resids = Y- Y_hat
    sigma = 1.4826* np.median(abs(resids- np.median(resids))) # 1.4826 * MAD
    resids_standard = (resids - resids.mean() ) / sigma
    beta = math.sqrt(2*math.log(len(X)))
    resids_truc = resids_standard[np.logical_not(resids_standard >=  -1*beta)]
    norm_one_sum = np.linalg.norm(resids_truc, ord=1)
    print "norm1 of residuals ", norm_one_sum
    norm_ones.append(norm_one_sum)
    fig = plt.figure(figsize=(10,6))
    plt.plot(X, resids_standard, '.k')
    plt.plot(X, np.zeros(len(X)) - beta, '-r' )
    plt.show()
    
    # plotting the result
    my_pdf = gaussian_kde(resids_standard)
    x = np.linspace(resids_standard.min(), resids_standard.max(), 100)
    plot(x,my_pdf(x),'r') # distribution function
    hist(resids_standard,normed=1,alpha=.3, bins = 100) # histogram
    show()
    
    ###############################
print zip(sample_list,norm_ones)
'''
print " "

# <codecell>


# <codecell>


# <codecell>

def proc_light_curve_pl(subset_index):
    import DAL
    #from DAL.datasets.checkpoint import Checkpoint
    lightcurves = DAL.create('lightcurves')
    #checkpoint = Checkpoint()
    if(subset_index >= len(s)):
        return []
    subset_filename = s[subset_index]
    subset_number = subset_filename[:-9]
    print subset_filename , subset_number
    ret = []
    ccc = 0
    for i in lightcurves.iter(subset_filename):
        ccc += 1
        #if(ccc >= 5):
        #    break
        #print "--"*10
        #t1 = time.clock()
        name = i['id']
        lc = i['data']
        time_x = lc[:int(len(lc)/2)]  #modified - first half time
        flux_y = lc[int(len(lc)/2):]  #modified - first half flux
        X = np.array(time_x, dtype = 'float32')
        Y = np.array(flux_y, dtype = 'float32')
        X = X[np.logical_not(np.isnan(Y))]
        Y = Y[np.logical_not(np.isnan(Y))]
        
        res = kernel_regress_cross_validation(X, Y)
        Y_hat = res[0]
        Y_band = res[2]
        
        resids = Y- Y_hat
        sigma = 1.4826* np.median(abs(resids- np.median(resids))) # 1.4826 * MAD
        resids_standard = (resids - resids.mean() ) / sigma
        beta = math.sqrt(2*math.log(len(X)))
        resids_truc = resids_standard[np.logical_not(resids_standard >=  -1*beta)]
        norm_one_sum = np.linalg.norm(resids_truc, ord=1)
        ret.append((name, norm_one_sum))
        '''
        print subset_index, name, Y_band, norm_one_sum
        print "temp_spent ",time.clock() - t1
        
        fig = plt.figure(figsize=(10,6))
        plt.plot(X, Y, '.k')
        plt.plot(X, Y_hat, '-b')
        title_string = "id:" + str(sample[0]) +", kernel regression, boxcar , bandwidth="+str(h)
        plt.title(title_string)
        plt.show()
        '''
        
    #checkpoint.store("proj4_1_"+str(subset_index), obj=ret) #also supports fp=<file pointer> and s=<string>
    return ret

# <codecell>

for i in checkpoint.list():
    print i

# <codecell>

dview['s'] = s
dview['nad_wat'] =nad_wat 
dview['nad_wat_with_L'] = nad_wat_with_L

dview['kernel_regress_cross_validation'] = kernel_regress_cross_validation
dview.execute("import math, time, scipy")
dview.execute("import numpy as np")

project_name = "proj4_2_test1_"


# <codecell>

#checkpoint is used in case of kernel crash
start = 0
try:
    start = checkpoint.load( project_name+"last_index_processed",t ="obj") + 1
    print "continued ", start
except:
    print start
# if start_over, start = 0, for duty memory checkpoint
#start = 0
print start

# <codecell>

#Use the map_sync multi time, each time parallelize the len(rc) files to slave nodes. 
#Then store the result to the checkpoints.

processed_file_index = []
merged_result = []

m = int(np.ceil(len(s)/float( len(rc))))
print "total number of parallelization " + str(m) + " with "+str(len(rc))+" slave nodes."

for i in range(start,m):
    print  "---", i, " of ", m-1
    temp=  i*len(rc) + np.array(range(len(rc)))
    if (i == m-1):
        temp = range(i*len(rc), len(s))
    t1 = time.clock()
    result_p = dview.map_sync(proc_light_curve_pl, temp)
    
    print "  to process subset() ", temp
    
    #use checkpoint to store the result
    for j in range(len(result_p)):
        re = result_p[j]
        try:
            file_name = project_name + s[temp[j]][:-9]
            #processed_file_name.append(file_name)
            checkpoint.store(file_name, obj = re) ## check point the result
            for k in re:
                merged_result.append(k)
            processed_file_index.append(temp[j])
        except:
            break
    # save the list of precessed file name
    #checkpoint.store( project_name+"processed", obj = processed_file_name)
    checkpoint.store( project_name+"last_index_processed", obj = i)
    

# <codecell>

end = checkpoint.load( project_name+"last_index_processed",t ="obj")

# <codecell>

print "all files finished? ", end == m-1
#for i in merged_result:
#    print i

# <codecell>



# <codecell>

# after deal all light curves
# the the following program
processed_file_name = []
for i in range(len(s)):
    file_name = project_name + s[i][:-9]
    processed_file_name.append(file_name)
#print processed_file_name
merged_result_checkpoint = []
for file_name in processed_file_name:
    #print file_name
    re = checkpoint.load(file_name, t = "obj")
    for j in re:
        #print j
        merged_result_checkpoint.append(j)
#print merged_result_checkpoint

print len(merged_result_checkpoint)

# <codecell>

from operator import itemgetter, attrgetter
RESULT = sorted(merged_result_checkpoint, key = itemgetter(1), reverse = True)
with open(project_name+'.txt', 'w+') as f:
    for Re in RESULT:
        #print Re
        f.write(str(Re[0])+'\n')

# <codecell>

print RESULT

# <codecell>

norm_ones = []
for i in RESULT:
    norm_ones.append(i[1])
#plot only 1%
plt.plot(range(1700),norm_ones[:1700])
plt.show()

# <markdowncell>

# We also tried norm 2, but the result is not better than norm1

