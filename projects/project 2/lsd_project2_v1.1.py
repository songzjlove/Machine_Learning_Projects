# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>
##########################
# LSD Projecct 2
# Zhengjian Song
# Part 1
##########################
import DAL
import time
import string
import re
import fractions
import numpy
import scipy
import time
from numpy import arange,array,ones,linalg
from pylab import plot,show
from __future__ import division
from IPython.parallel import Client
import matplotlib.pyplot as plt 

rc = Client()
dview = rc[:]
wishes=DAL.create('wishes')
data = wishes.subsets()[13:18] #use only one recent week of twitter data to get the volcabulary,  

# Clean the raw text data with filters
# Create a very long string comprising of the first 5 days of twitter data
dictionary= {}
for i in range(len(data)):
    print 'day', i
    text = ""
    for tweet in wishes.iter(data[i]):
        if tweet.has_key('text'):
            lower = tweet['text'].lower()
            text += lower
        else:
            #print i
            break
    text = re.sub(r'@[a-zA-Z0-9_]+', "", text) #removes the @bullshit
    text = re.sub(r'[^\x00-\x7F]', "", text) #removes weird shit
    text = re.sub(r'&[a-zA-Z0-9_]+', "", text)#removes &amp etc.
    text = re.sub(r'\W', " ", text) #removes anything which is not alphanumeric
    text = re.sub(r'[0-9]', "", text) #removes #s
    text = re.sub(r"\s+", " ", text) #removing white space
    #print text
    #splitting text into words
    tt = text.split()
    for word in tt:
        if len(word) >= 2:
            if dictionary.has_key(word):
                dictionary[word] = dictionary[word] + 1
            else:
                dictionary.setdefault(word, 1)
    
#print len(dictionary.keys()) #
#Trim the dictionary, each word has at least 25 occurances.
V = []
for word in dictionary.keys():
    if dictionary[word] > 25:
        V.append(word)
M = len(V)
print 'get the V(representation) and its size equals ', len(V)
print V

# <markdowncell>

# The feature representation we chose is a vocabulary of words. We used the tweets from the week before the testing data (i.e. Dec 20th onwards). We used the wishes data to create the dictionary since in the original printout of the instructions, we were told to look at the wishes data. We also chose to use just the tweets from only the week before the testing data since we did not feel that using the entire set of training data would be useful. Thus we will be looking at each of the tweets from Dec 20 onwards and then cleaning it up and checking of the invidual words which comprise the tweets can be found in the vocabulary. If so, we update the feature vector by 1 in that spot taken by the word and 0 if not. 

# <codecell>

# LSD Projecct 2
# Michelle Yeo and Zhengjian Song
# Part 2

# L2-regularization online update beta function. 
def online_update_beta_l2(y_label, mu, xvec, beta, b0, eta, rambda):
    new_b0 = 0
    new_beta = np.zeros(len(beta))
    pp = 1/(1+exp(-1* mu))
    new_b0 = b0 + eta*((y_label - pp) - rambda)
    new_beta = beta + eta * ((y_label - pp )* xvec - rambda)
    return (new_b0,new_beta)

# Different learning rates
def learning( ty, t):
    # Switch
    if ty == 0:
        return 0.1
    elif ty == 1:
        return 0.5
    elif ty == 2:
        return 1.0
    #
    elif ty == 3:
        return 0.3/t
    elif ty == 4:
        return 1/t
    elif ty == 5:
        return 2/t
    #
    elif ty == 6:
        return 0.3/sqrt(t)
    elif ty == 7:
        return 1/sqrt(t)
    elif ty == 8:
        return 5/sqrt(t)
    
    else:
        return 1

# <markdowncell>

# NOTE: The training is fast enough, but loading the data will take a while.
# One possible way is to load all data at first, and then train models online with different learning rate, but this will lead to the kernel crash (out of memory) problem.
# Or one could train the model with different learning rate with loading the data only once, but the code will be ugly with 9 sets of error counts variable. Besides, you will see the results together in the end not one by one, which means, it will really take you a while to see the result to check the codes.
# 
# Wanna be faster? please try the code:
# test = lwishes.subsets()[20:24]

# <codecell>

# this function trains one set of online data using stochastic gradient descent, it updates the parameters that we have gotten already
def online_ridge_logistic(learning_type, g_beta, g_b0, g_t, f):
    #parameters, L, l, M --- tweets of one new day.
    #parameters, V --- vocabulary
    #parameters, f --- file to write predicted label for unlabeled tweets
    #parameters, errors stuff --- stat the errors
    #parameters, g_beta, g_b0 --- global parameter
    
    # the rambda here usually should be from CV, for online training?
    rambda = 0.001 #lambda penalty for ridge logistic regression
    
    errors = 0
    #using SGD to update betas
    for k in xrange(L):
    #for k in range(20):
        #normalising tweet i.e cleaning it up
        temp = l[k]['text']
        temp = re.sub(r'@[a-zA-Z0-9_]+', "", temp)
        temp = re.sub(r'[^\x00-\x7F]', "", temp)
        temp = re.sub(r'&[a-zA-Z0-9_]+', "", temp)
        temp = re.sub(r'\W', " ", temp)
        temp = re.sub(r'[0-9]', "", temp) #removes #s
        temp = temp.lower()
        temp = re.sub(r"\s+", " ", temp) #removing white space
        temp = temp.split()
        #print temp
        #creating betas
        x = np.zeros(M)
        
        for word in temp:
            if word in V:
                j = V.index(word)
                x[j] = 1
        mu = sum(g_beta * x) + g_b0
        
        '''
        
        ind = [] #list of indices eg. [34,78,90]
        for word in temp:
            if word in V:
                j = V.index(word)
                ind.append(j)
                x[j] = 1
        betas=[]
        if not ind: #checking for empty list
            continue
        else:
            for i in range(len(ind)):
                betas.append(g_beta[ind[i]])
        #print 'betas are', betas
        #prediction and writing out files
        mu = sum(betas) + g_b0
        '''
        pi  = exp(mu)/ (1 + exp(mu))
        y_predict = 0
        if pi > 0.5:
            y_predict = 1
        #y_predict = get_y(pi)
        
        
        if l[k]['label'] != '?':
            ## the tweet is labeled, then
            y_label = l[k]['label']
            if y_predict != y_label:
                errors = errors + 1
            else:
                pass # correct
            # no matter what we predict, we update the parameters
            # update step
            g_t = g_t + 1
            eta = learning(learning_type, g_t)
            # the function to online update the parameters
            updated_beta = online_update_beta_l2(y_label, mu, x, g_beta, g_b0, eta, rambda)
            
            #print 'true label is ', l[k]['label'] 
            #print temp
            #print updated_beta[1][ind] - beta[ind]
    
            g_b0 = updated_beta[0]
            g_beta = updated_beta[1]
            
        else: ## for unlabeled data
            id_and_label = "<"+str(l[k]['id'])+"> <"+str(y_predict)+"> "
            f.write("%s" % id_and_label)
        
        #print 'true label is', l[k]['label']
        #print id_and_label
        #print 'error count is', error_count
    return [g_beta, g_b0, g_t, errors]





# <codecell>

# input the data:
lwishes=DAL.create('wishes-labelled')
##testing the other data
test = lwishes.subsets()[20:22] 
print 'The number of days to predict on is', len(test)
i = 0
g_beta = np.zeros(M)
g_b0 = 0
g_error_vec = []
g_t = 0 #counting the number of iterations
g_ones = 0# count of labelled one data
f = open('predict_unlabelled.txt','w')

learn_type = 1 # 1,2,3,4,5,6,7,8,9, different learning rate respectively
print 'Here day 0 -- Dec 20, 2012, day 1 --Dec 21, 2012, ....'

error_rate_temp = []
total_errors = 0.0
total_tweet = 0
l = []
t1 = time.clock()
for day in test: # 44 days in total
    i = i + 1
    for tweet in lwishes.iter(day):
        l.append(tweet)
        total_tweet += 1
        if total_tweet%1000 == 0:
            #print '--'
            #print total_tweet
            L = len(l)
            result = online_ridge_logistic(learn_type, g_beta, g_b0, g_t, f)
            g_beta = result[0]
            g_b0 = result[1]
            g_t = result[2]
            total_errors = total_errors + result[3]
            error_rate_temp.append(total_errors / g_t)
            #print total_errors / g_t
            l = []
    print 'after day ' + str(i) + " " + str(total_tweet) + " got trained" 
# for the last part less then 1000
L = len(l)
result = online_ridge_logistic(learn_type, g_beta, g_b0, g_t, f)
g_beta = result[0]
g_b0 = result[1]
g_t = result[2]
total_errors = total_errors + result[3]
t2 = time.clock()
print t2 - t1
f.close()
print 'The average error rate for labeled data is', total_errors / g_t
plot(range(int(total_tweet / 1000)), error_rate_temp)
show()


# <codecell>


# <codecell>



# <codecell>


# <codecell>



# <codecell>



# <codecell>


# <codecell>



# <codecell>



# <codecell>


# <markdowncell>

# Short commentary:
#     
#     Based on our results, it seems like the best step size is 0.3/sqrt(t), which gives us an error rate of 0.030. This was only slightly better than most of the other rates which were around that level as well. The worst rates we got were from step sizes of 1, 0.3/t and 1/t which were near 0.08. 

# <markdowncell>

# 
#     

