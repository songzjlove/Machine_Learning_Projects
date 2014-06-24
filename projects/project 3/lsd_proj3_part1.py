# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

# STAT 376
# project 3 
# Zhengjian Song

# import package
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

# import the dataset
crime = DAL.create('crime')
crime_list = crime.get_crime_list()

crime_counts = crime.get_crime_counts()
region_list = crime.get_region_list()
K = 10 #number of crime types
N = int(len(crime_counts)/K) #number of regions
T = len(crime_counts.get((100,0)))
print K,N,T

proj3_1_filename = []
proj3_2_filename = []
for i in range(10):
    proj3_1_filename.append("lsd_proj3_prob1_["+ str(i)+ string.lower(crime_list.get(i)) + "]epan.txt")
    proj3_2_filename.append("lsd_proj3_prob2_["+ str(i)+ string.lower(crime_list.get(i)) + "]epan.txt")
print proj3_1_filename
print proj3_2_filename

rc = Client()
dview = rc[:]



# <codecell>

def loglikelihood_poisson(weeks, lamb, N, C): # threshold = 1e-06
    import math
    # lamb N*K
    loglikelihood = 0
    for i in xrange(N): #
        lam = lamb[i]
        if lam < 1e-10:
            lam = 1e-10
        counts = V[i]
        for j in weeks:
            loglikelihood += -1* lam + counts[j] * math.log(lam) # + math.log(math.factorial(counts[j])) # delete to let it faster
    return loglikelihood

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

# <codecell>

C = 0 # type of crime {0, 1,2,3,4,...9 }K = 10
#creating a crime count matrix for all districts over all time for a fixed crime type C
V= np.zeros((N,T)) #creating a NxT 0 matrix
for j in xrange(N):
    V[j] = np.array(crime_counts.get((j,C)))

# <codecell>


# this function chooses the best pair of bandwidth based on the results of predicting the neighboring time windows
def search_bandwidth(index):
    
    loglike_result = []
    N = 2985
    T = 643
    
    # this is the Kernel matrix for 
    K_dist = np.zeros((N,N))  
    # the regions list, easy to use
    position = [] 
    for i in xrange(N):
        position.append(region_list.get(i))
    # the neighboring days
    pred_days = [ T-52*1-0, T-52*1-1, T-52*1+1,\
                              T-52*2-0, T-52*2-1, T-52*2+1,\
                              T-52*3-0, T-52*3-1, T-52*3+1,\
                              T-52*4-0, T-52*4-1, T-52*4+1,\
                              T-52*5-0, T-52*5-1, T-52*5+1,\
                              T-52*0-2, T-52*0-1]
    
    #search for the bandwidth for distance
    for hd in [index]:
        h_dist = math.exp(H_dist[hd]) #check this... should it be log instead of exp?
        #get the Kernel matrix of distance (outputs a symmetric matrix)
        for i in xrange(N):
            for j in xrange(i,N):
                K_dist[i][j]= K_gaussian(position[i],position[j], h_dist)
                K_dist[i][j]= K_dist[j][i]
        
        ## search for the bandwidth for time
        for ht in xrange(len(H_time)):
            h_time = math.exp(H_time[ht])
            # moving windows predicition
            # for last a few years, on this time
            
            total_loglike = 0
            loglikelihood_array = []
            for day in pred_days:
                #for each week to predict, estimate the lamb_hat
                K_time =  np.zeros(day)
                for s in xrange(day):
                    K_time[s] = K_gaussian([s],[day],h_time)
                K_time_grid = np.zeros((N,day))
                for j in xrange(N):
                    K_time_grid[j] = K_time
                    
                # smoothing & scaling
                lamb_hat = np.dot(np.dot(K_dist, V[:,0:day]), K_time)
                scales =  np.dot(np.dot(K_dist, K_time_grid), np.ones(day)) #scaling, denom
                for i in xrange(N):
                    lamb_hat[i] = lamb_hat[i] / scales[i]
                    # the thresholding may be different for different crime, but we set it roughly here
                    if lamb_hat[i] <= 1e-3 or math.isnan(lamb_hat[i]):
                        lamb_hat[i] = 1e-3
                loglikelihood_array.append( loglikelihood_poisson(range(day, day+1), lamb_hat, N, C))
            # get the mean of log likelihood for all neighboring time windows.
            mean_loglikehood = np.mean(np.array(loglikelihood_array))
            loglike_result.append( ((hd, ht), mean_loglikehood))
    return loglike_result

# <codecell>

# Simulation to choose lamb for each square
# since we the data radius has been fixed, so 
# (-5, + 5) is the region for guassian of distance 
# (-2, + 11) is the region for guassian of time 
# (-4, + 3) is the region for epan of distance 
# (0,  + 10) is the region for epan of time

H_dist = np.arange(-5, 4, 1) #range of bandwidth values for distance, that is enough
H_time = np.arange(-2, 11, 1) #range of bandwidth values for time, that is enough

# passing the parameters and functions, get prepared for parallelization
dview['H_dist'] = H_dist
dview['H_time'] = H_time

dview['V'] = V
dview['region_list'] = region_list
dview['K_boxcar'] = K_boxcar
dview['K_gaussian'] = K_gaussian
dview['K_epan'] = K_epan
dview['loglikelihood_poisson'] = loglikelihood_poisson
dview['C'] = C

dview.execute("import math, time, scipy")
dview.execute("import numpy as np")


# <codecell>

# para
re = dview.map_sync(search_bandwidth, range(len(H_dist)))

# <codecell>

#print re

# <codecell>

max = -10000000
max_i = 0
max_j = 0
log_likelihood_surf = np.zeros((len(H_dist),len(H_time)))
for i in re:
    for j in i:
        index = j[0]
        log_likelihood_surf[index[0], index[1]] = j[1]
        if(j[1] > max):
            max = j[1]
            max_i = index[0]
            max_j = index[1]
#print log_likelihood_surf


from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure()
ax = fig.gca(projection='3d')
X = H_time
Y =  H_dist
X, Y = np.meshgrid(X, Y)
Z = log_likelihood_surf
surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
        linewidth=0, antialiased=False)
#ax.set_zlim(-1.01, 1.01)

ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_wireframe(X, Y, Z, rstride=1, cstride=1)

plt.show()

# <codecell>

print max_i, max_j, max

# <codecell>

h_dist_op = math.exp(H_dist[max_i])
h_time_op = math.exp(H_time[max_j])


#train the lambda
K_dist = np.zeros((N,N))
position = [] #all the regions
for i in xrange(N):
    position.append(region_list.get(i))


h_dist = math.exp(H_dist[max_i]) #check this... should it be log instead of exp?
#get the Kernel matrix of distance (outputs a symmetric matrix)
for i in xrange(N):
    for j in xrange(i,N):
        K_dist[i][j]= K_gaussian(position[i],position[j], h_dist)
        K_dist[i][j]= K_dist[j][i]

h_time = math.exp(H_time[max_j])

K_time =  np.zeros(T)
for s in xrange(T):
    K_time[s] = K_gaussian([s],[T],h_time)
K_time_grid = np.zeros((N,T))
for j in xrange(N):
    K_time_grid[j] = K_time
    
#smoothing 
thred = 1e-3
lamb_hat = np.dot(np.dot(K_dist, V), K_time)
scales =  np.dot(np.dot(K_dist, K_time_grid), np.ones(T)) #scaling, denom
for i in xrange(N):
    if lamb_hat[i] == 0 or math.isnan(lamb_hat[i]):
        lamb_hat[i] = 1e-3
    elif scales[i] != 0:
        lamb_hat[i] = lamb_hat[i] / scales[i]
    else :
        lamb_hat[i] = 1e-3
print list(lamb_hat)
    
## writing out to file



# <codecell>


# <codecell>

print C
s=""
for i in range(N):
    f = str(i)+","+str(lamb_hat[i])+'\n'
    s+=f
#print s
    
with open(proj3_1_filename[C], 'w') as g:
    g.write(s)
g.close

# <codecell>




# <codecell>


