# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# This file contains the program for lsd_project4 part 3 and part 4. Considering that the part 4 is to improve the result of part 2, and our algorithm for part 3 is based on the result from part2, we just combine part 3 and part 4 into this program. 
# 
#     

# <markdowncell>

# Brief introduction:
#     We only re-rank the top curves in the ranking from part 2. What we did is to tell three types of curves from these curves: the curves with white noise or non-period signals, the one with only one level strength of signals, and the ones with more than 2 levels of strength of signals. One basic hypothesis is that most stars with one planet have one level strength of signal, while eclipsing binary stars have more than more levels of signals. This may be not always true, but considering only ten curves are given as examples, we can not capture all differences of these two types of star systems. When we try to detect the levels of signals, another job in the part 4 is done, that is to pick out the most possible curves of not these two types. We place curves with one level on the top of new ranking, which are followed by the curves with more than two levels of signals, and finally put all those with noise or non-period curves in the end of new ranking. 
#     To implement this idea. 
#     First, we choose another kind of fit instead of local average using Boxcar kernel. We use the local median to fit the curves so that we get fitted curves that are more robust to extreme responses, Y s. 
#     Second, we take outliers with large residuals, as raw signals to process. Note that here the threshold we used is not as tight as part 2 to let more potential signals. 
#     Third, after scaling and other simple processing, we re-sharp the signals into square waves using a filter. In this filter, we just choose 0111 as open-signal and start to perceive the signals after one open-signal, and then close it till a close-signal as 1110. We also use the max-value as the magnitude of the wave. In this way, we filter most noises out and get neat sharps of square signals. Usually, the after filtering, there would be not discrete noises. However, we find one special situation where there is a period of explosive noise or something else on only one side of curve. If this leads to most false positive cases, we need to classify such situation into the white noise types. 
#     Fourth, we detect a safe zone (band) between two levels of signals to get the upper signals, and re-scale the left signals, and detect another safe zone to get another level out. 
#     Finally, get a new ranking of ids. 
#     
#     Attention: this algorithm works only when the different levels of strength singals are good features to tell the conf and ep. Bearing in mind that there is a limited number of labelled curves, we think our hypothesis requires further consideration.

# <codecell>

## Lsd project 4 part 3+4
## Zhengjian Song
## please check the analysis in the markdown above to see our main ideas
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
# we just go through the data set get the reversed dictionary id:subset_index
reversed_dict = checkpoint.load("reversed_dict", t = "obj")

# <codecell>

#This function fit the curve using local median, instead of local average
#So it will be more robust than local average to fit these light curves.
def nad_wat_robust(X, Y, h):
    #def nad_wat_robust(X, Y, h):
    N = len(X)
    Y_med = np.zeros(N, dtype='float32')
    for i in range(N):
        x = (X[i] - X) / h
        # Boxcar
        K = np.array(abs(x) <= 1, dtype='float32') * 1
        #print np.sum(K)
        Y_local_median = np.median(Y[np.logical_not(K < 1)])
        #print len(Y[np.logical_not(K > 1)])
        Y_med[i] = Y_local_median
        
    return (Y_med, 0)

# <codecell>

#print reversed_dict

# <codecell>


# <codecell>


# <markdowncell>

# The following window shows a small demo of testing the wave_filter function to filter the levels of signals strength. 
# As we can see, it can properly tell curves with different levels of signalss, without checking the frequencies of signals of curves. At least, it works for this purpose.

# <codecell>

# small demo result
# we can see that we successfully tell the curves only with one wave from white noise and more than one wave

# <codecell>


# <codecell>


# <codecell>

## small demo to show the fit of several curves

import time
def wave_filter(t_id):

    import DAL
    from DAL.datasets.checkpoint import Checkpoint
    #from DAL.datasets.checkpoint import Checkpoint
    lightcurves = DAL.create('lightcurves')
   
    try:
        reversed_dict
    except:
        checkpoint = Checkpoint()
        reversed_dict = checkpoint.load("reversed_dict", t = "obj")
    
    #get the date by id
    sample_data = []
    subset_index = reversed_dict.get(t_id)
    if(subset_index < 0 or subset_index>= len(s)):
        return subset_index
    for i in lightcurves.iter(s[subset_index]):
        name = i['id']
        if name == t_id:
            print name
            lc = i['data']
            time = lc[:int(len(lc)/2)]  #modified - first half time
            flux = lc[int(len(lc)/2):]  #modified - first half flux
            sample_data.append((name, time, flux))
            break
    
    
    for sample in [sample_data[0]]:
        print "-"*10
        print "id:", sample[0]
        X = np.array(sample[1], dtype = 'float32')
        Y = np.array(sample[2], dtype = 'float32')
        X = X[np.logical_not(np.isnan(Y))]
        Y = Y[np.logical_not(np.isnan(Y))]
        h = 0.8 # be carefully chosen
        #res = nad_wat(X, Y, 0.1)
        #fit the curve
        res = nad_wat_robust(X,Y, h)
        Y_hat = res[0]

        resids = Y- Y_hat
        sigma = 1.4826* np.median(abs(resids- np.median(resids))) # 1.4826 * MAD
        resids_standard = (resids - resids.mean() ) / sigma
        beta = math.sqrt(2*math.log(len(X)))
        
        
        ## detect the different levels, if there are more levels, just think they are binary stars.
        wave_level = 0
        outliers = np.logical_not(resids_standard >=  -1*0.9*beta) # let more potential signals in
        
        #scale the outlier residuals to range(0,1)
        mask = (np.zeros(len(X))+1)* outliers # len (mask) == len(X) , if signal, mask[i] == 1,, not signal, mask[i] == 0
        
        if(mask.sum() <= 3): # no signal, it is very possible to be outliers from white noise.
            # mark this to be white noise
            print "LEVEL: ", wave_level
            print "==="*10
            return wave_level
            
        # get the signal out
        signal = resids_standard * mask

        # rescale the signal to process
        signal = -1 * signal # get them reflected above the x-axis
        s_min = signal[ signal > 0].min() # the min non-zero 
        s_max = signal.max() # the max
        
        if(abs(s_min - s_max) <= 0.0001):
            #it is almost impossible to happen, just in case
            print "LEVEL: ", wave_level
            print "==="*10
            return wave_level
        
        #scaled signal
        signal = signal - s_min
        signal = signal / abs(s_max - s_min)
        

        #stat the area in each cell.
        #resharp the signals into rect
        total_num_signal = np.sum(mask)
        re_signal = np.zeros(len(signal))
        
        is_open = False
        signal_strength = 0
        for i in range(4, len(X)-4):
            # start to receive
            if( mask[i-1] == 0 and mask[i] == 1 and mask[i+1] == 1 and mask[i+2] == 1):
                is_open = True
                
            if(is_open == True):
                if(signal[i] > signal_strength):
                    signal_strength = signal[i]
            # start to refuse
            if( mask[i-2] == 1 and mask[i-1] == 1 and mask[i] == 1 and mask[i+1] == 0):
                j = i
                while(j >= 0 and mask[j] > 0):
                    re_signal[j] = signal_strength
                    j = j -1
                signal_strength = 0
                is_open = False
        
        if(re_signal.max() <= 0.0001): # no signal, it is very possible to be outliers from white noise.
            # mark this to be white noise
            print "no signal after filtering"
            print "LEVEL: ", wave_level
            print "==="*10
            
            return wave_level
            
        #re-scale after filter
        re_signal = re_signal / abs(re_signal.max() - re_signal.min())
        
        left_sum = 0
        right_sum = 0
        left_sum = np.sum(re_signal[:int(len(X)/2)])+1
        right_sum = np.sum(re_signal[int(len(X)/2)+1:])+1
        if(left_sum > 2* right_sum or right_sum > 2* left_sum ):
            print "unbalanced signal"
            print "LEVEL: ", wave_level
            print "==="*10
            return wave_level
        
        
        # confirm that there is at least one level of signals
        wave_level = 1 
        
        ########################################################
        #get the upper level of signal
        
        #find the safe zone between waves
        zone_width = 0.1
        up_bound = np.arange(1.0, 0.12, -0.025)
        down_bound = up_bound - zone_width
        
        up_zone = 1
        down_zone = 1
        for i in range(len(up_bound)):
            ones = np.zeros(len(X))+1
            is_in_zone = np.logical_and(re_signal <= up_bound[i], re_signal >= down_bound[i])
            num_in_zone = np.sum(is_in_zone * 1)
            if(num_in_zone == 0):
                up_zone = up_bound[i]-0.01
                down_zone = down_bound[i]+0.01
                break

        
        is_upper = np.logical_not(re_signal <= up_zone) # let more potential signals in
        #scale the outlier residuals to range(0,1)
        mask = (np.zeros(len(X))+1)* is_upper # len (mask) == len(X) , if signal, mask[i] == 1,, not signal, mask[i] == 0
        
        up_signal = re_signal * mask
        down_signal = re_signal - up_signal
        
        
        up_is_true_signal = True
        down_is_true_signal = True 
        # testing the upside is singel
        left_sum = 0
        right_sum = 0
        left_sum = np.sum(up_signal[:int(len(X)/2)])+1
        right_sum = np.sum(up_signal[int(len(X)/2)+1:])+1
        if(left_sum > 2* right_sum or right_sum > 2* left_sum ):
            print "unbalanced up signal"
            up_is_true_signal = False
            return wave_level
            
        # testing the downside is singel
        if(down_signal.max() <= 0.0001): # no signal, it is very possible to be outliers from white noise.
            # mark this to be white noise
            print "no signal down_side "
            down_is_true_signal = False
            return wave_level
    
        down_signal = down_signal / abs(down_signal.max())
        left_sum = 0
        right_sum = 0
        left_sum = np.sum(down_signal[:int(len(X)/2)])+1
        right_sum = np.sum(down_signal[int(len(X)/2)+1:])+1
        if(left_sum > 2* right_sum or right_sum > 2* left_sum ):
            print "unbalanced up signal"
            down_is_true_signal = False
            return wave_level
     
        if(up_is_true_signal and down_is_true_signal):
            wave_level += 1
        
        
        #############################################
        #continue if you want to get more levels
        #############################################
        re_signal = down_signal
        
        #re-scale after filter
        re_signal = re_signal / abs(re_signal.max() - re_signal.min())
        
    
        ########################################################
        #get the upper level of signal
        
        #find the safe zone between waves
        zone_width = 0.1
        up_bound = np.arange(1.0, 0.12, -0.025)
        down_bound = up_bound - zone_width
        
        up_zone = 1
        down_zone = 1
        for i in range(len(up_bound)):
            ones = np.zeros(len(X))+1
            is_in_zone = np.logical_and(re_signal <= up_bound[i], re_signal >= down_bound[i])
            num_in_zone = np.sum(is_in_zone * 1)
            if(num_in_zone == 0):
                up_zone = up_bound[i]-0.01
                down_zone = down_bound[i]+0.01
                break

        is_upper = np.logical_not(re_signal <= up_zone) # let more potential signals in
        #scale the outlier residuals to range(0,1)
        mask = (np.zeros(len(X))+1)* is_upper # len (mask) == len(X) , if signal, mask[i] == 1,, not signal, mask[i] == 0
        
        up_signal = re_signal * mask
        down_signal = re_signal - up_signal
        
        
        up_is_true_signal = True
        down_is_true_signal = True 
        # testing the upside is singel
        left_sum = 0
        right_sum = 0
        left_sum = np.sum(up_signal[:int(len(X)/2)])+1
        right_sum = np.sum(up_signal[int(len(X)/2)+1:])+1
        if(left_sum > 2* right_sum or right_sum > 2* left_sum ):
            print "unbalanced up signal"
            up_is_true_signal = False
            return wave_level
            
        # testing the downside is singel
        if(down_signal.max() <= 0.0001): # no signal, it is very possible to be outliers from white noise.
            # mark this to be white noise
            print "no signal down_side "
            down_is_true_signal = False
            return wave_level
    
        down_signal = down_signal / abs(down_signal.max())
        left_sum = 0
        right_sum = 0
        left_sum = np.sum(down_signal[:int(len(X)/2)])+1
        right_sum = np.sum(down_signal[int(len(X)/2)+1:])+1
        if(left_sum > 2* right_sum or right_sum > 2* left_sum ):
            print "unbalanced up signal"
            down_is_true_signal = False
            return wave_level
     
        if(up_is_true_signal and down_is_true_signal):
            wave_level += 1
        
        
        return wave_level
    #end of for sample

# <codecell>

#parallelizaiton
dview['s'] = s

dview['reversed_dict'] = reversed_dict
dview['nad_wat_robust'] =nad_wat_robust 
dview['wave_filter'] = wave_filter
dview.execute("import math, time, scipy")
dview.execute("import numpy as np")

f = open('proj4_2_test1_.txt', 'r') # read ranking list from result in part 2
ranking = []
for i in f:
    ranking.append(int(i))
#print ranking

id_list = ranking
project_name = "proj4_3_4_test2_"

# <codecell>

top_wave0 = []
top_wave1 = []
top_wave2 = []
finished_top_ids = 0

# <codecell>

# check point is used to restart.
start = 0
try:
    start = checkpoint.load( project_name+"last_index_processed",t ="obj") + 1
    print "continued ", start
except:
    print start
# if start_over, start = 0, for duty memory checkpoint
start = 0
print start

# <codecell>

TopN = int(900)
step = len(rc)*18
m = int(np.ceil(TopN/float(step)))
print "total number of parallelization " + str(m) + " with "+str(len(rc))+" slave nodes."
for i in range(start,m):
    print  "---", i, " of ", m-1

    temp=  i*step + np.array(range(step))
    if (i == m-1):
        temp = range(i*step, TopN)
        
    print temp
    
    ids_list = ranking[temp[0]: temp[len(temp)-1]+1]
    
    print len(ids_list)
    
    result_p = dview.map_sync(wave_filter, ids_list)
    
    merge_result = []
    for re in result_p:
        merge_result.append(re)
            
    for j in range(len(ids_list)):
        if(merge_result[j] == 0):
            top_wave0.append(ids_list[j])
        elif(merge_result[j] == 1 or merge_result[j] == 3 ):
            top_wave1.append(ids_list[j])
        elif(merge_result[j] == 2):
            top_wave2.append(ids_list[j])
        else:
            top_wave1.append(ids_list[j])
    
    checkpoint.store( project_name+"_wave0", obj = top_wave0)
    checkpoint.store( project_name+"_wave1", obj = top_wave1)
    checkpoint.store( project_name+"_wave2", obj = top_wave2)
    finished_top_ids += len(id_list)
    print finished_top_ids
    checkpoint.store( project_name+"_finished_top_ids", obj = finished_top_ids)
    
    print "  to process curve() ", ids_list
    #print result_p
    #use checkpoint to store the result
    #checkpoint.store( project_name+"processed", obj = processed_file_name)
    checkpoint.store( project_name+"last_index_processed", obj = i)
    

# <codecell>

#top_wave0 = checkpoint.load( project_name+"_wave0", t= "obj")
print top_wave0

# <codecell>

#top_wave1 = checkpoint.load( project_name+"_wave1", t= "obj")
print top_wave1

# <codecell>

#top_wave2 = checkpoint.load( project_name+"_wave2", t= "obj")
print top_wave2

# <codecell>


# <codecell>


# <codecell>


# <codecell>

wave0 = top_wave0
wave1 = top_wave1
wave2 = top_wave2

new_wave2 = []
for w in wave2:
    if w in new_wave2:
        pass
    else:
        new_wave2.append(w)
print len(new_wave2)

new_wave0 = []
for w in wave0:
    if w in new_wave0:
        pass
    else:
        new_wave0.append(w)
print len(new_wave0)

new_wave1 = []
for w in wave1:
    if w in new_wave1:
        pass
    else:
        new_wave1.append(w)
print len(new_wave1)

# <codecell>

print len(new_wave1 + new_wave2) + len(new_wave0)

# <codecell>


# <codecell>

f = open('proj4_2_test1_.txt', 'r')
ranking = []
for i in f:
    ranking.append(int(i))
#print ranking

end_of_processed = len(new_wave1 + new_wave2) + len(new_wave0)

new_ranking = []
left_ranking = ranking[end_of_processed:]
new_ranking.extend(new_wave1)
new_ranking.extend(new_wave2)
new_ranking.extend(new_wave0)
new_ranking.extend(left_ranking)


# <codecell>


# <codecell>

lightcurves.score(new_ranking)

# <codecell>

from operator import itemgetter, attrgetter
RESULT = new_ranking
with open(project_name+'_f.txt', 'w+') as f:
    for Re in RESULT:
        f.write(str(Re)+'\n')
f.close()

# <codecell>

## codes for the small demo above

import time
def wave_filter(t_id):

    import DAL
    #from DAL.datasets.checkpoint import Checkpoint
    lightcurves = DAL.create('lightcurves')
    
    
    #get the date by id
    sample_data = []
    subset_index = reversed_dict.get(t_id)
    if(subset_index < 0 or subset_index>= len(s)):
        return -1
    for i in lightcurves.iter(s[subset_index]):
        name = i['id']
        if name == t_id:
            print name
            lc = i['data']
            time = lc[:int(len(lc)/2)]  #modified - first half time
            flux = lc[int(len(lc)/2):]  #modified - first half flux
            sample_data.append((name, time, flux))
            break
    
    
    for sample in [sample_data[0]]:
        print "-"*10
        print "id:", sample[0]
        X = np.array(sample[1], dtype = 'float32')
        Y = np.array(sample[2], dtype = 'float32')
        X = X[np.logical_not(np.isnan(Y))]
        Y = Y[np.logical_not(np.isnan(Y))]
        h = 0.8 # be carefully chosen
        #res = nad_wat(X, Y, 0.1)
        #fit the curve
        res = nad_wat_robust(X,Y, h)
        Y_hat = res[0]
        
        ###################################
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
        
        #################################################
        resids_truc = resids_standard[np.logical_not(resids_standard >=  -1*beta)]
        norm_one_sum = np.linalg.norm(resids_truc, ord =2 ) # norm 2
        print "norm1 of residuals ", norm_one_sum
        #norm_ones.append(norm_one_sum)
        fig = plt.figure(figsize=(10,6))
        plt.plot(X, resids_standard, '.k')
        plt.plot(X, np.zeros(len(X)) - beta, '-r' )
        plt.show()
        ###############################
        
        
        ## detect the different levels, if there are more levels, just think they are binary stars.
        wave_level = 0
        outliers = np.logical_not(resids_standard >=  -1*0.9*beta) # let more potential signals in
        
        #scale the outlier residuals to range(0,1)
        mask = (np.zeros(len(X))+1)* outliers # len (mask) == len(X) , if signal, mask[i] == 1,, not signal, mask[i] == 0
        
        if(mask.sum() <= 3): # no signal, it is very possible to be outliers from white noise.
            # mark this to be white noise
            print "LEVEL: ", wave_level
            print "==="*10
            return wave_level
            
        # get the signal out
        signal = resids_standard * mask

        # rescale the signal to process
        signal = -1 * signal # get them reflected above the x-axis
        s_min = signal[ signal > 0].min() # the min non-zero 
        s_max = signal.max() # the max
        
        if(abs(s_min - s_max) <= 0.0001):
            #it is almost impossible to happen, just in case
            print "LEVEL: ", wave_level
            print "==="*10
            return wave_level
        
        #scaled signal
        signal = signal - s_min
        signal = signal / abs(s_max - s_min)
        
        ##########################
        plt.plot(X, signal,'.k')
        plt.show()
        #########################
        
        #stat the area in each cell.
        #resharp the signals into rect
        total_num_signal = np.sum(mask)
        re_signal = np.zeros(len(signal))
        
        is_open = False
        signal_strength = 0
        for i in range(4, len(X)-4):
            # start to receive
            if( mask[i-1] == 0 and mask[i] == 1 and mask[i+1] == 1 and mask[i+2] == 1):
                is_open = True
                
            if(is_open == True):
                if(signal[i] > signal_strength):
                    signal_strength = signal[i]
            # start to refuse
            if( mask[i-2] == 1 and mask[i-1] == 1 and mask[i] == 1 and mask[i+1] == 0):
                j = i
                while(j >= 0 and mask[j] > 0):
                    re_signal[j] = signal_strength
                    j = j -1
                signal_strength = 0
                is_open = False
        
        if(re_signal.max() <= 0.0001): # no signal, it is very possible to be outliers from white noise.
            # mark this to be white noise
            print "no signal after filtering"
            print "LEVEL: ", wave_level
            print "==="*10
            
            return wave_level
            
        #re-scale after filter
        re_signal = re_signal / abs(re_signal.max() - re_signal.min())
        
        left_sum = 0
        right_sum = 0
        left_sum = np.sum(re_signal[:int(len(X)/2)])+1
        right_sum = np.sum(re_signal[int(len(X)/2)+1:])+1
        if(left_sum > 2* right_sum or right_sum > 2* left_sum ):
            print "unbalanced signal"
            print "LEVEL: ", wave_level
            print "==="*10
            return wave_level
        
        ########################################################
        plt.ylim(-0.2, 1,2)
        plt.plot(X, re_signal,'.g')
        ########################################################
        
        # confirm that there is at least one level of signals
        wave_level = 1 
        
        ########################################################
        #get the upper level of signal
        
        #find the safe zone between waves
        zone_width = 0.1
        up_bound = np.arange(1.0, 0.12, -0.025)
        down_bound = up_bound - zone_width
        
        up_zone = 1
        down_zone = 1
        for i in range(len(up_bound)):
            ones = np.zeros(len(X))+1
            is_in_zone = np.logical_and(re_signal <= up_bound[i], re_signal >= down_bound[i])
            num_in_zone = np.sum(is_in_zone * 1)
            if(num_in_zone == 0):
                up_zone = up_bound[i]-0.01
                down_zone = down_bound[i]+0.01
                break
                
        ########################################################
        plt.plot(X, np.zeros(len(X)) + up_zone, '-y')
        plt.plot(X, np.zeros(len(X)) + down_zone, '-y')
        ########################################################
        
        is_upper = np.logical_not(re_signal <= up_zone) # let more potential signals in
        #scale the outlier residuals to range(0,1)
        mask = (np.zeros(len(X))+1)* is_upper # len (mask) == len(X) , if signal, mask[i] == 1,, not signal, mask[i] == 0
        
        up_signal = re_signal * mask
        down_signal = re_signal - up_signal
        
        ########################################################
        plt.plot(X,up_signal, ".b")
        plt.plot(X,down_signal, ".r")
        ########################################################
        
        up_is_true_signal = True
        down_is_true_signal = True 
        # testing the upside is singel
        left_sum = 0
        right_sum = 0
        left_sum = np.sum(up_signal[:int(len(X)/2)])+1
        right_sum = np.sum(up_signal[int(len(X)/2)+1:])+1
        if(left_sum > 2* right_sum or right_sum > 2* left_sum ):
            print "unbalanced up signal"
            up_is_true_signal = False
            return wave_level
            
        # testing the downside is singel
        if(down_signal.max() <= 0.0001): # no signal, it is very possible to be outliers from white noise.
            # mark this to be white noise
            print "no signal down_side "
            down_is_true_signal = False
            return wave_level
    
        down_signal = down_signal / abs(down_signal.max())
        left_sum = 0
        right_sum = 0
        left_sum = np.sum(down_signal[:int(len(X)/2)])+1
        right_sum = np.sum(down_signal[int(len(X)/2)+1:])+1
        if(left_sum > 2* right_sum or right_sum > 2* left_sum ):
            print "unbalanced up signal"
            down_is_true_signal = False
            return wave_level
     
        if(up_is_true_signal and down_is_true_signal):
            wave_level += 1
        
        #continue if you want to get more levels
        re_signal = down_signal
        
        plt.show()
        
        print "LEVEL: ", wave_level
        print "==="*10
        
        return wave_level
    #end of for sample

# <codecell>

f = open('proj4_2_test1_.txt', 'r')
ranking = []
for i in f:
    ranking.append(int(i))
#print ranking

#sample_list = [(189,914193), (126,2042713), (276,6852096), (102,3096237), (313,7412246), (231, 485585)]
#id_list = ranking[0:50]
id_list = [914193, 2042713, 6852096, 3972533, 1255171, 5004731, 4848370, 6258272, 3096237, 7412246]
#id_list = top_wave2
for t_id in id_list:
    print "WAVE-LEVEL-------", wave_filter(t_id)

# <codecell>


