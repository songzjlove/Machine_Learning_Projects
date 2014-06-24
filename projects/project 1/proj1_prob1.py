###############################
## Project 1 problem 1
## Stat 376
## Zhengjian Song
###############################

# Verify cluster connection and connect to the data set
import DAL
from IPython.parallel import Client
rc = Client()
#print len(rc)
dview = rc[:]
tinyimages = DAL.create('tinyimages')

# Search the image set by keywords: car and bicycle
car_ids = tinyimages.search('car', 2000)
bicycle_ids = tinyimages.search('bicycle', 2000)

# Ground truth, constructed manually
# 100 in each category
car_true = [12025562,12025563,12025564,12025565,12025567,12025568,12025569,12025571,12025572,12025573,12025574,12025576,12025580,12025583,12025584,12025585,12025586,12025587,12025588,12025589,12025590,12025591,12025592,12025593,12025594,12025597,12025599,12025601,12025602,12025603,12025604,12025605,12025606,12025611,12025615,12025618,12025620,12025624,12025627,12025628,12025630,12025631,12025632,12025633,12025634,12025635,12025636,12025661,12025663,12025664,12025668,12025670,12025671,12025674,12025675,12025676,12025677,12025678,12025679,12025682,12025683,12025684,12025686,12025687,12025688,12025689,12025690,12025691,12025692,12025693,12025701,12025702,12025706,12025707,12025709,12025711,12025712,12025713,12025717,12025718,12025725,12025726,12025727,12025747,12025749,12025752,12025753,12025757,12025768,12025769,12025771,12025772,12025773,12025774,12025776,12025778,12025779,12025780,12025781,12025782]
bicycle_true = [7112211,7112212,7112213,7112214,7112215,7112216,7112218,7112219,7112220,7112223,7112224,7112225,7112226,7112227,7112228,7112229,7112231,7112232,7112234,7112235,7112237,7112239,7112240,7112241,7112243,7112244,7112245,7112246,7112247,7112248,7112249,7112250,7112251,7112253,7112259,7112260,7112261,7112262,7112263,7112265,7112266,7112268,7112270,7112272,7112273,7112275,7112276,7112281,7112285,7112287,7112297,7112300,7112301,7112302,7112303,7112304,7112309,7112310,7112312,7112320,7112325,7112328,7112329,7112330,7112335,7112339,7112341,7112342,7112343,7112344,7112347,7112348,7112349,7112350,7112351,7112352,7112355,7112358,7112359,7112360,7112361,7112362,7112365,7112366,7112371,7112374,7112380,7112381,7112386,7112394,7112403,7112409,7112413,7112414,7112416,7112421,7112424,7112426,7112427,7112431]
 
# Print the results
images = tinyimages.byid(car_ids[0:300])
print "\"cars\""
tinyimages.display(images)
images = tinyimages.byid(car_true)
print "verified cars"
tinyimages.display(images)
print "\n"

images = tinyimages.byid(bicycle_ids[0:300])
print "\"bicycles\""
tinyimages.display(images)
images = tinyimages.byid(bicycle_true)
print "verified bicycles"
tinyimages.display(images)
print "\n"

images = tinyimages.byid([7112211])
tinyimages.display(images)

# Store the result in files
with open('car.txt', 'w') as f:
    f.write("%s" % car_ids)
with open('true_car.txt', 'w') as f:
    f.write("%s" % car_true)
with open('bicycle.txt', 'w') as f:
    f.write("%s" % bicycle_ids)
with open('true_bicycle.txt', 'w') as f:
    f.write("%s" % bicycle_true)
