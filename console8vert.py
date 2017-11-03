from algRealEmbeddings import *
import time

v2, v3, v1 = [0, -2.3418234786589167, 0], [0, 6.327508996139455, 0], [3.697464444401205, 0.017900844941222793, 0]
v4 = [3.2930168732119, -0.157944798200452, -1.95429223125632]
v5 = [2.82039958291193, 0.679610905369317, -1.67465679235168]
v6 = [2.47296030555702, 0.342686069384568, -1.46758104891375]
v7 = [-2.63035589169772, 0, -10.6657888907962]
v8 = [1.0, 1.0, 1.0]

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
           
#start = time.time()
#G = GraphEmbedding(lengths, 'Max8vertices')
#sols = G.findEmbeddings()
#print len(sols['real'])
#print len(sols['complex'])
#end = time.time()
#print 'Time: ' , end - start
alg = AlgRealEmbeddings('Max8vertices', num_phi=16, num_theta=16,  name='8vert')
alg.findMoreEmbeddings(lengths)












