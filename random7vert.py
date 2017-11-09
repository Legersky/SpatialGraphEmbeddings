from algRealEmbeddings import *
from graphEmbedding import *
from random import uniform

def dist( u, v):
    return float(np.sqrt( (u[0]-v[0])**2 + (u[1]-v[1])**2 + (u[2]-v[2])**2))

a = -0.05
b = 0.05

#for i in range(0, 10):
#v1 = [uniform(a, b), 0.5 + uniform(a, b), -1.0+ uniform(a, b)]
#v2 = [-1.0+uniform(a, b), uniform(a, b), uniform(a, b)]
#v3 = [1.0+uniform(a, b), uniform(a, b), uniform(a, b)]
#v4 = [-1.0+uniform(a, b), uniform(a, b), uniform(a, b)]
#v5 = [1.0+uniform(a, b), uniform(a, b), uniform(a, b)]
#v6 = [uniform(a, b), 1.0+uniform(a, b), uniform(a, b)]
#v7 = [uniform(a, b),0.5 +  uniform(a, b), 1.0+uniform(a, b)]

#v1 = [0.6 + uniform(a, b), 0.6 + uniform(a, b), -3.0+ uniform(a, b)]
#v2 = [-1.0+uniform(a, b), uniform(a, b), uniform(a, b)]
#v3 = [1.0+uniform(a, b), uniform(a, b), uniform(a, b)]
#v4 = [-1.0+uniform(a, b), uniform(a, b), uniform(a, b)]
#v5 = [1.0+uniform(a, b), uniform(a, b), uniform(a, b)]
#v6 = [1.0+uniform(a, b), 1.0+uniform(a, b), uniform(a, b)]
#v7 = [0.6 + uniform(a, b), 0.6 + uniform(a, b), 1.0+uniform(a, b)]

v1 = [0.6, 0.6, -3.0]
v2 = [-1.0+uniform(a, b), uniform(a, b), uniform(a, b)]
v4 = [-1.0+uniform(a, b), uniform(a, b), uniform(a, b)]
v3 = [1.0+uniform(a, b), uniform(a, b), uniform(a, b)]
v5 = [1.0+uniform(a, b), uniform(a, b), uniform(a, b)]
v6 = [1.0, 1.0, 0]
v7 = [0.6, 0.6, 1.0]


lengths = {(1, 2) : dist(v1,v2),
           (1, 3) : dist(v1,v3),
           (1, 4) : dist(v1,v4),
           (1, 5) : dist(v1,v5),
           (1, 6) : dist(v1,v6),
           (7, 2) : dist(v7,v2),
           (7, 3) : dist(v7,v3),
           (7, 4) : dist(v7,v4),
           (7, 5) : dist(v7,v5),
           (7, 6) : dist(v7,v6),
            (2, 3) : dist(v2,v3), 
            (3, 4) : dist(v3,v4), 
            (4, 5) : dist(v4,v5), 
            (5, 6) : dist(v5,v6), 
            (2, 6) : dist(v2,v6), 
            }

#
#lengths = {(1, 2): 2.515622197507309, (2, 7): 15.537215472935126, (4, 7): 16.731525140348634, (2, 6): 3.5052783567172123, (6, 7): 13.357673687397279, 
#(5, 6): 1.767285448511689, (5, 7): 12.19172535570313, (1, 4): 17.34066043914577, (1, 5): 2.9172680066093295, (1, 3): 2.6518433506436776, 
#(1, 6): 1.8382523757977096, (4, 5): 16.295943525182871, (3, 7): 15.867510999992877, (3, 4): 18.684441964056564, (2, 3): 1.1492992481675544}
alg = AlgRealEmbeddings('Max7vertices',  name='random7vert_tbr')
#    alg.findMoreEmbeddings_tree(lengths, onlyOne=True)
alg.findMoreEmbeddings(lengths)










