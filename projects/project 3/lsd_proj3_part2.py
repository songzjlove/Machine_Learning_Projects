# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

# STAT 37600, proj 3, part 2
# Zhengjian Song

# this program improved the result in part1 with best bandwidth for each region (not one set bandwidth for all regions) and seasonal factor added.

# import the packages
import DAL
import time

import string
import re
import fractions
import numpy
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
import numpy as np


crime = DAL.create('crime')
crime_list = crime.get_crime_list()
# the following can be slow; do this once at the beginning
# of the program and use this data structure throughout
crime_counts = crime.get_crime_counts()
region_list = crime.get_region_list()
K = 10 #number of crime types
N = int(len(crime_counts)/K) #number of regions
T = len(crime_counts.get((100,0)))

print K,N,T
print crime_list


proj3_1_filename = []
proj3_2_filename = []
for i in range(10):
    proj3_1_filename.append("lsd_proj3_prob1_["+ str(i)+ string.lower(crime_list.get(i)) + "]rim.txt")
    proj3_2_filename.append("lsd_proj3_prob2_["+ str(i)+ string.lower(crime_list.get(i)) + "]rim.txt")
print proj3_1_filename
print proj3_2_filename

rc = Client()
dview = rc[:]

# <codecell>


C = 0 # type of crime {0, 1,2,3,4,..9}.K = 10}
#creating a crime count matrix for all districts over all time for a fixed crime type C
V= np.zeros((N,T)) #creating a NxT 0 matrix
for j in xrange(N):
    V[j] = np.array(crime_counts.get((j,C)))


# <codecell>

def loglikelihood_poisson(weeks, lamb, N, C): # threshold = 1e-06
    import math
    # lamb N*K
    loglikelihood = 0
    for i in xrange(N): #
        lam = lamb[i]
        if lam < 1e-12:
            lam = 1e-12
        counts = V[i]
        for j in weeks:
            loglikelihood += -1* lam + counts[j] * math.log(lam) # + math.log(math.factorial(counts[j])) # delete to let it faster
    return loglikelihood

def loglikelihood_poisson_each(weeks, lamb, N, C):
    import math
    import numpy as np
    # lamb N*K
    loglikelihood = np.zeros(N)
    for i in xrange(N): #
        lam = lamb[i]
        if lam < 1e-12:
            lam = 1e-12
        counts = V[i]
        temp_loglikelihood = 0
        for j in weeks:
            temp_loglikelihood += -1* lam + counts[j] * math.log(lam) # + math.log(math.factorial(counts[j])) # delete to let it faster
        loglikelihood[i] = temp_loglikelihood
    return loglikelihood
#----------------------------------------------------

# define the kernel functions:
def K_gaussian(x, y, h):
    import math
    norm2 = 0
    #check the length, then
    for i in xrange(len(x)):
        norm2 += math.pow((x[i]-y[i]),2)
    return math.exp(- norm2 / (2*h*h)) / (h *math.sqrt(2*math.pi))

#print K_gaussian(region_list.get(2),region_list.get(0),h)

def K_epan(x, y, h):
    import math
    norm2 = 0
    #check the length, then
    for i in xrange(len(x)):
        norm2 += math.pow((x[i]-y[i]),2)
    ret = 0
    #print norm2
    if math.sqrt(norm2) <=h :
        return 0.75*(1 - norm2/(h*h))
    else:
        return 0
    return 0

#print K_epan(region_list.get(2),region_list.get(0),h)

def K_boxcar(x,y, h):
    import math
    norm2 = 0
    #check the length, then
    for i in xrange(len(x)):
        norm2 += math.pow((x[i]-y[i]),2)
    ret = 0
    if math.sqrt(norm2) <= h:
        return 1/float(2*h)
    else:
        return 1
    return 1
#print K_boxcar((1,0),(0,0), 2)

#----------------------------------------------------

# <codecell>

# this function chooses the best pair of bandwidth for all chicago regions, based on the results of predicting the neighboring time windows
# 2 changes: the first one is that, now we assume the best bandwidth set vary from region to another, so for each regions, we get best bandwidth.
# the seconde one is that, we added the h_year parameter, to add the seasonal factor. 
#          to do this, we re- sharp the one dimensional time series into two-dimensional, where one-dimension is along the weeks, the other is along years
def search_bandwidth_each(index):
    
    loglike_result = []
    N = 2985
    T = 643
    K_dist = np.zeros((N,N))  
    position = [] #all the regions
    for i in xrange(N):
        position.append(region_list.get(i))
    
    #search for the bandwidth for distance
    for hd in [index]:
        h_dist = math.exp(H_dist[hd]) #check this... should it be log instead of exp?
        #get the Kernel matrix of distance (outputs a symmetric matrix)
        for i in xrange(N):
            for j in xrange(i,N):
                K_dist[i][j]= K_gaussian(position[i],position[j], h_dist)
                K_dist[i][j]= K_dist[j][i]
        
        ## search for the bandwidth for time
        for hy in xrange(len(H_year)):
            
            h_year = math.exp(H_year[hy])
            
            for ht in xrange(len(H_time)):
                h_time = math.exp(H_time[ht])
                # moving windows predicition
                # for last a few years, on this time
                pred_days = [ T-52*1-0, T-52*1-1, T-52*1+1,\
                              T-52*2-0, T-52*2-1, T-52*2+1,\
                              T-52*3-0, T-52*3-1, T-52*3+1,\
                              T-52*4-0, T-52*4-1, T-52*4+1,\
                              T-52*5-0, T-52*5-1, T-52*5+1,\
                              T-52*0-2, T-52*0-1]
                
                total_loglike = 0
                loglikelihood_array = np.zeros((len(pred_days), N))

                for ii in xrange(len(pred_days)):
                    day = pred_days[ii]
                    #for each week to predict, estimate the lamb_hat
                    
                    K_time =  np.zeros(day)
                    mapping = []
                    year = 0
                    for s in xrange(day + 1):
                        week = s % 52
                        if week == 0:
                            year += 1
                        mapping.append((year*h_year, week))
                    
                    for s in xrange(day):
                        K_time[s] = K_gaussian(mapping[s],mapping[day],h_time)
                    K_time_grid = np.zeros((N,day))
                    for j in xrange(N):
                        K_time_grid[j] = K_time
                        
                    #smoothing 
                    lamb_hat = np.dot(np.dot(K_dist, V[:,0:day]), K_time)
                    scales =  np.dot(np.dot(K_dist, K_time_grid), np.ones(day)) #scaling, denom
                    for i in xrange(N):
                        lamb_hat[i] = lamb_hat[i] / scales[i]
                        if lamb_hat[i] <= 1e-3 or math.isnan(lamb_hat[i]):
                            lamb_hat[i] = 1e-3
                    
                    loglikelihood_array[ii] =  loglikelihood_poisson_each(range(day, day+1), lamb_hat, N, C)
                    
                mean_loglikehood = loglikelihood_array.mean(axis = 0) # for each square, mean of the 11 weeks
                loglike_result.append( ((hd, ht, hy), mean_loglikehood, lamb_hat))
            
    return loglike_result


# <codecell>

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
# train the mle as base line

pred_days = [ T-52*1-0, T-52*1-1, T-52*1+1,\
                          T-52*2-0, T-52*2-1, T-52*2+1,\
                          T-52*3-0, T-52*3-1, T-52*3+1,\
                          T-52*4-0, T-52*4-1, T-52*4+1,\
                          T-52*5-0, T-52*5-1, T-52*5+1,\
                          T-52*0-2, T-52*0-1]
mean_li = []
lamb_mle = np.zeros((K,N))
for day in pred_days:
    # get a max_likelihood (approximately )on training set
    for i in xrange(N):
        counts = V[i,0:day]
        lamb_mle[C][i] = mean(counts)
    min_lamb = 10000
    for i in xrange(N):
        if lamb_mle[C][i] <= 1e-3:
            continue
        else:
            if min_lamb > lamb_mle[C,i]:
                min_lamb = lamb_mle[C,i]
    for i in xrange(N):
        if lamb_mle[C,i] <= 1e-3:
            lamb_mle[C,i] = min_lamb
            
    #print lamb_mle[C,:]
    t1 = time.clock()
    mean_li.append(loglikelihood_poisson(range(day,day+1), lamb_mle[C], N, C))
#print "the log_likehood baseline", mean(mean_li)


# <codecell>

# Simulation to choose lamb for each square

# got -5 to + 5 is the region for guassian of distance 
# got -2 to + 11 is the region for guassian of time 
# got -4 to +3 is the region for epan of distance 
# got 0 to 10 is the region for epan of time

H_dist = np.arange(-5, 4, 1) #range of bandwidth values for distance, that is enough
H_time = np.arange(-2, 8, 1) #range of bandwidth values for time, that is enough
H_year = np.arange(-5, 1, 1)
print len(H_dist)
print len(H_time)
dview['H_dist'] = H_dist
dview['H_time'] = H_time
dview['H_year'] = H_year

dview['V'] = V
dview['region_list'] = region_list

dview['K_boxcar'] = K_boxcar
dview['K_gaussian'] = K_gaussian
dview['K_epan'] = K_epan
#dview['loglikelihood_poisson'] = loglikelihood_poisson
dview['loglikelihood_poisson_each'] = loglikelihood_poisson_each

dview['C'] = C

dview.execute("import math, time, scipy")
dview.execute("import numpy as np")

# <codecell>

re = dview.map_sync(search_bandwidth_each, range(len(H_dist)))

# <codecell>

# merging the results

MIN = -1e10
lamb_hat = np.zeros(N)
lamb_max_loglik = np.zeros(N) + MIN
lamb_h_dist = np.zeros(N)
lamb_h_time = np.zeros(N)
lamb_h_year = np.zeros(N)


jj = 2450
surface_jj = np.zeros((len(H_dist), len(H_time)))

for i in re:
    for j in i:
        #for each pair of bandwidth:
        temp_mle = j[1]
        temp_lamb_hat = j[2]
        for k in range(N):
            if lamb_max_loglik[k] <  temp_mle[k] :
                lamb_max_loglik[k] = temp_mle[k]
                lamb_hat[k] = temp_lamb_hat[k]
                lamb_h_dist[k] = j[0][0]
                lamb_h_time[k] = j[0][1]
                lamb_h_year[k] = j[0][2]
            if k == jj:
                surface_jj[j[0][0], j[0][1]] = temp_mle[jj]
                
# the lamb_hat done
print C
print list(sort(lamb_hat))
plt.plot(lamb_hat)
plt.show()
print "test result: ",np.sum(loglikelihood_poisson_each(range(T-1,T), lamb_hat, N, C))
#print list(sort(lamb_hat))
plt.plot(lamb_mle[C])
print "mle baseline: ",np.sum(loglikelihood_poisson_each(range(T-1,T), lamb_mle[C], N, C))
#print list(sort(lamb_mle[C]))

# <codecell>

# check how the h_dist and h_time vary for different squares
print lamb_h_dist[jj], lamb_h_time[jj]
print lamb_h_dist.max(), lamb_h_dist.min()
print lamb_h_time.max(), lamb_h_time.min()

r = np.histogram(lamb_h_dist, bins = range(len(H_dist)))
print r
plt.plot(r[0])
plt.show()

r = np.histogram(lamb_h_time, bins = range(len(H_time)))
print r
plt.plot(r[0])
plt.show()

r = np.histogram(lamb_h_year, bins = range(len(H_time)))
print r
plt.plot(r[0])
plt.show()


# <codecell>

print list(sort(lamb_hat))

# <codecell>

## writing out to file 2
print C
v=np.arange(0, 21, 1)
print v

def pois_cdf(k,lamb_hat):
    f = math.exp(-lamb_hat)
    l=0
    for i in range(k+1):
        l+=math.pow(lamb_hat,i)/math.factorial(i)
    g=l*f
    return g

s=""
for i in range(N):
    cdf_list=[]
    for k in range(len(v)):
        #print k
        temp = pois_cdf(k, lamb_hat[i])
        cdf_list.append(temp)
    #print len(cdf_list)
    f = str(i)
    for k in range(len(cdf_list)):
        f += ","+str(cdf_list[k])
    f += '\n'
    s+=f
#print s

with open(proj3_2_filename[C], 'w') as g:
    g.write(s)
g.close

# <codecell>

print C
print lamb_hat.mean()
print lamb_hat.max()
print "test result: ",np.sum(loglikelihood_poisson_each(range(T-1,T), lamb_hat, N, C))
print "mle baseline: ",np.sum(loglikelihood_poisson_each(range(T-1,T), lamb_mle[C], N, C))

plt.plot( lamb_mle[C])
plt.show()
plt.plot(lamb_hat)
plt.show()


