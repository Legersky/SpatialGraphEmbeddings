from algRealEmbeddings import *
#import time
from random import uniform

v2, v3, v1 = [0, -2.3418234786589167, 0], [0, 6.327508996139455, 0], [3.697464444401205, 0.017900844941222793, 0]
v4 = [3.2930168732119, -0.157944798200452, -1.95429223125632]
v5 = [2.82039958291193, 0.679610905369317, -1.67465679235168]
v6 = [2.47296030555702, 0.342686069384568, -1.46758104891375]
v7 = [-2.63035589169772, 0, -10.6657888907962]
v8 = [1.0, 1.0, 1.0]
a = -0.1
b = 0.1

#v1 = [uniform(a, b), uniform(a, b), uniform(a, b)]
#v7 = [uniform(a, b), uniform(a, b), uniform(a, b)]
#
#v2 = [1.0+uniform(a, b), uniform(a, b), uniform(a, b)]
#v5 = [1.0+uniform(a, b), uniform(a, b), uniform(a, b)]
#
#v3 = [uniform(a, b), 1.0+uniform(a, b), uniform(a, b)]
#
#v6 = [uniform(a, b), uniform(a, b), 1.0+uniform(a, b)]
#
#v4 = [1.0+uniform(a, b), 1.0+uniform(a, b), uniform(a, b)]
#
#v8 = [ 1.0 + uniform(a, b),uniform(a, b), 1.0+uniform(a, b)]


#v1 = [uniform(a, b), 1.0 + uniform(a, b), uniform(a, b)]
#v8 = [uniform(a, b), 1.0 + uniform(a, b), uniform(a, b)]

#v4 = [-1.0+uniform(a, b), uniform(a, b), uniform(a, b)]
#v2 = [-1.0+uniform(a, b), uniform(a, b), uniform(a, b)]

#v3 = [1.0+uniform(a, b), uniform(a, b), uniform(a, b)]
#v5 = [1.0+uniform(a, b), uniform(a, b), uniform(a, b)]

#v6 = [uniform(a, b), 0.5 + uniform(a, b), 1.0+uniform(a, b)]

#v7 = [uniform(a, b), 0.5 + uniform(a, b), -1.0+uniform(a, b)]

def dist( u, v):
    return float(np.sqrt( (u[0]-v[0])**2 + (u[1]-v[1])**2 + (u[2]-v[2])**2))

lengths = {
           (2, 1) : dist(v2,v1),
            (2, 7) : dist(v2,v7),
            (2, 6) : dist(v2,v6),
            (3, 1) : dist(v3,v1),
            (3, 7) : dist(v3,v7),
            (3, 2) : dist(v3,v2),
            (4, 1) : dist(v4,v1),
            (4, 7) : dist(v4,v7),
            (4, 3) : dist(v4,v3),
            (5, 1) : dist(v5,v1),
            (5, 7) : dist(v5,v7),
            (5, 4) : dist(v5,v4),
            (6, 1) : dist(v6,v1),
            (6, 5) : dist(v6,v5),
            (5, 8) : dist(v5,v8),
            (6, 8) : dist(v6,v8),
            (7, 8) : dist(v7,v8),
            (2, 8) : dist(v2,v8),
           }

#lengths = {(1, 2): 4.386290254987835, (4, 7): 10.5579448289947, (2, 6): 2.1695833997369136, (4, 5): 1.6940062148563753, (2, 8): 16.83405853179157,
#           (1, 4): 1.8084396572052768, (1, 3): 11.026883517351973, (1, 6): 5.8450130495737325, (3, 7): 13.787489684821617, (5, 8): 13.075209742298672, 
#           (2, 7): 11.232184203671995, (6, 8): 18.368267879186273, (1, 5): 2.00289249524296, (7, 8): 18.762183274653285, (5, 7): 10.536273659997796, 
#           (5, 6): 5.918182159410676, (2, 3): 11.29840753979362, (3, 4): 10.81796534240893}
#start = time.time()
#G = GraphEmbedding(lengths, 'Max8vertices')
#sols = G.findEmbeddings()
#print len(sols['real'])
#print len(sols['complex'])
#end = time.time()
#print 'Time: ' , end - start
alg = AlgRealEmbeddings('Max8vertices', name='8vert')
alg.findMoreEmbeddings(lengths)












