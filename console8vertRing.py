from algRealEmbeddings import *
import time
#from random import uniform

lengths_7 = {(1, 2): 6.274873914397456, (2, 7): 12.093170928755793, (4, 7): 10.52881316000657, (2, 6): 10.479941653322763, (6, 7): 8.11321631042846, 
             (5, 6): 9.459324566316262, (5, 7): 10.73159467196586, (1, 4): 1.9761239449935557, (1, 5): 2.8483166740669206, (1, 3): 1.9564840951854308, 
             (1, 6): 9.4046159048512159, (4, 5): 2.313950594269001, (3, 7): 10.525064001461296, (3, 4): 0.47038007668728976, (2, 3): 5.952813550679546}
#/home/jan/Research/MaximumReal3dEmbeddings/CouplerCurve/important/7vert/48_embd_from_28_7vert_53129440.txt
G = GraphEmbedding(lengths_7, 'Max7vertices')
v1, v2, v3, v4, v5, v6, v8 = G.getEmbedding()
v7 = [0, 0, 0]

def dist( u, v):
    return float(np.sqrt( (u[0]-v[0])**2 + (u[1]-v[1])**2 + (u[2]-v[2])**2))

lengths = {
            (1, 2) : dist(v1,v2),
            (1, 3) : dist(v1,v3),
            (1, 4) : dist(v1,v4),
            (1, 5) : dist(v1,v5),
            (1, 6) : dist(v1,v6),
            (1, 7) : dist(v1,v7),
            (2, 3) : dist(v2,v3),
            (2, 7) : dist(v2,v7),
            (2, 8) : dist(v2,v8),
            (3, 4) : dist(v3,v4),
            (3, 8) : dist(v3,v8),
            (4, 5) : dist(v4,v5),
            (4, 8) : dist(v4,v8),
            (5, 6) : dist(v5,v6),
            (5, 8) : dist(v5,v8),
            (6, 8) : dist(v6,v8),
            (7, 8) : dist(v7,v8),
            (6, 7) : dist(v6,v7),
           }

#start = time.time()
#G = GraphEmbedding(lengths, 'Ring8vertices')
#sols = G.findEmbeddings()
#print len(sols['real'])
#print len(sols['complex'])
#end = time.time()
#print 'Time: ' , end - start


alg = AlgRealEmbeddings('Ring8vertices', name='8Ring')
alg.findMoreEmbeddings(lengths)












