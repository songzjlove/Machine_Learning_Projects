##############################
## Project 1 problem 2 + 3
## Stat 376
## Zhengjian Song
##############################
## problem 2
## Bible: If you can handle the program in 5 nodes cluster, never use the bullshit 10 nodes one!!!
#verify cluster connection
from IPython.parallel import Client
import DAL

rc = Client()
#print len(rc)
dview = rc[:]

tinyimages = DAL.create('tinyimages')

car_ids = tinyimages.search('car', 2000)
bicycle_ids = tinyimages.search('bicycle', 2000)
#ground truth, constructed manually
#you will have 100 in each category
car_true = [12025562,12025563,12025564,12025565,12025567,12025568,12025569,12025571,12025572,12025573,12025574,12025576,12025580,12025583,12025584,12025585,12025586,12025587,12025588,12025589,12025590,12025591,12025592,12025593,12025594,12025597,12025599,12025601,12025602,12025603,12025604,12025605,12025606,12025611,12025615,12025618,12025620,12025624,12025627,12025628,12025630,12025631,12025632,12025633,12025634,12025635,12025636,12025661,12025663,12025664,12025668,12025670,12025671,12025674,12025675,12025676,12025677,12025678,12025679,12025682,12025683,12025684,12025686,12025687,12025688,12025689,12025690,12025691,12025692,12025693,12025701,12025702,12025706,12025707,12025709,12025711,12025712,12025713,12025717,12025718,12025725,12025726,12025727,12025747,12025749,12025752,12025753,12025757,12025768,12025769,12025771,12025772,12025773,12025774,12025776,12025778,12025779,12025780,12025781,12025782]
bicycle_true = [7112211,7112212,7112213,7112214,7112215,7112216,7112218,7112219,7112220,7112223,7112224,7112225,7112226,7112227,7112228,7112229,7112231,7112232,7112234,7112235,7112237,7112239,7112240,7112241,7112243,7112244,7112245,7112246,7112247,7112248,7112249,7112250,7112251,7112253,7112259,7112260,7112261,7112262,7112263,7112265,7112266,7112268,7112270,7112272,7112273,7112275,7112276,7112281,7112285,7112287,7112297,7112300,7112301,7112302,7112303,7112304,7112309,7112310,7112312,7112320,7112325,7112328,7112329,7112330,7112335,7112339,7112341,7112342,7112343,7112344,7112347,7112348,7112349,7112350,7112351,7112352,7112355,7112358,7112359,7112360,7112361,7112362,7112365,7112366,7112371,7112374,7112380,7112381,7112386,7112394,7112403,7112409,7112413,7112414,7112416,7112421,7112424,7112426,7112427,7112431]

#@dview.remote(block = TRUE)
def gist(index):
    import DAL
    import scipy
    import leargist
    import math
    import numpy as np
    tinyimages=DAL.create('tinyimages')
    img=scipy.misc.toimage( \
        tinyimages.byid(index).reshape(32,32,3, order="F").copy())
    #return leargist.color_gist(img)
    a = leargist.color_gist(img)
    vec = np.array(a)
    sd = math.sqrt(np.var(vec))
    mu = np.mean(vec)
    #normalization
    standarized_gist =  (vec- mu) / sd 
    #In this case the gist of 960
    #ret = [[index]]
    return list(standarized_gist)
################################ 
################################
################################
first_ids_list = car_ids
second_ids_list =  bicycle_ids
################################
# be aware that map_sync function returns in the same order of first_order
gists_car_ids=dview.map_sync(gist,first_ids_list)
gists_bicycle_ids = dview.map_sync(gist, second_ids_list)

################################
## The results are in gists_car_ids and gists_bicycle_ids
#print len(gists_bicycle_ids)
print 'problem 2'
print 'the results are in gists_car_ids, with each mean = 0, var = 1'
print len(gists_car_ids)
print np.mean(np.array(gists_car_ids[0]))
print np.var(np.array(gists_car_ids[0]))

gists = gists_car_ids + gists_bicycle_ids
gists_ids = first_ids_list + second_ids_list
dview['gists']= gists
#dview['ids'] = gists_ids

###########################################################################################################
# problem 3
## Codes without parallelization are not desirable, because the computation complexity is O(N^2 * M), 
## I scatter the works to engines

Num = len(rc) # number of engines
N = len(gists)
M = len(gists[1]) # gist vector size
#print N
size = int(N / float(Num) + 0.999999) # ceiling of the float
down_size = int(size / float(2) + 0.999999)
up_size = size - down_size
#print size
# the assignment list of ranges for engines,

#print li
#Since the scatter funciton of dview just allocate the even parameters in order. 
#So just recoder to assign "evenly" heavy burdens for engines
reordered_li = [] 
down = N
up = 0
for i in range(Num):
    reordered_li.append( (down - down_size, down) )
    down = down - down_size
    reordered_li.append( (up, up + up_size))
    up = up + up_size
#print reordered_li
print reordered_li
dview.scatter('from_ids_range' , reordered_li)
#print N

@dview.remote(block = True)
def wc():
    import numpy as np
    #distance_map = {}
    ret = []
    from_ids_index = []
    N = len(gists)
    for i in range(len(from_ids_range)):
        start = from_ids_range[i][0]
        end = from_ids_range[i][1]
        if start < 0:
            start = 0
        if end > N:
            end = N
        from_ids_index = from_ids_index + range(start, end)
    
    for i in from_ids_index:
        for j in range(N):
            if i < N and i <= j:
                dist = 0
                #Euclidian Distance:
                diff = np.array(gists[i]) - np.array(gists[j])
                dist = sum( diff * diff)
                ret.append([(i,j), dist])
                #distance_map.setdefault((i,j), dist) !!Crap to build Hash table for large data
                #distance_map.setdefault((j,i), dist)
    #end of for i in from...
    
    return ret
    #return from_ids_index
    
wcs = wc()

N = len(gists_ids) # the indexing for images id ---- index in Distance_Matrix
Distance_Matrix = [] # or in P_Dist_Array to save memory(i <= j)
for i in range(N):
    Distance_Matrix.append([0.0 for i in range(N)])
    
for i in range(len(wcs)):
    List = wcs[i]
    for j in range(len(List)):
        J = List[j]
        x = J[0][0]
        y = J[0][1]
        if x >= N or y >= N:
            print '~~~~~~~'
            break       
        Distance_Matrix[x][y] = J[1]
        Distance_Matrix[y][x] = J[1]
    #end of j
#end of j
###########


k = 10
import heapq  
import scipy
def smallest(n, data):
    return heapq.nsmallest(n, data)
kNN_list = []
for i in range(N):
    s_zip =zip(Distance_Matrix[i], range(N))
    kNN_list.append(  smallest(k+1, s_zip)[1:] )
#end of for i 


##
seq = range(N)
find_id = dict(zip(gists_ids, seq)) # posting table of #gists_ids


#test ten cars and ten bicycles
test_cars_ids = car_true[0:10]
test_bicy_ids = bicycle_true[0:10] # id
test_cars_index = [] # index
test_bicy_index = []
for i in range(10):
    test_cars_index.append(find_id.get(car_true[i]))
    test_bicy_index.append(find_id.get(bicycle_true[i]))


Test_Dict1 = []
print 'cars'
for i in range(10):
    index = test_cars_index[i]
    from_id = gists_ids[index]
    #print from_id
    to_id = []
    for j in kNN_list[index]:
        to_id.append( gists_ids[j[1]])
    Test_Dict1.append((from_id, to_id))
    print '-'*10
    print from_id
    images = tinyimages.byid([from_id])
    tinyimages.display(images)
    print to_id
    images = tinyimages.byid(to_id)
    tinyimages.display(images)

Test_Dict2 = []
print 'bicycles'
for i in range(10):
    index = test_bicy_index[i]
    from_id = gists_ids[index]
    #print from_id
    to_id = []
    for j in kNN_list[index]:
        to_id.append( gists_ids[j[1]])
    Test_Dict2.append((from_id, to_id))
    print '-'*10
    print from_id
    images = tinyimages.byid([from_id])
    tinyimages.display(images)
    print to_id
    images = tinyimages.byid(to_id)
    tinyimages.display(images)