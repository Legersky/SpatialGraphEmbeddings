from algRealEmbeddings import *
from graphEmbedding import *
from random import uniform

def dist( u, v):
    return float(np.sqrt( (u[0]-v[0])**2 + (u[1]-v[1])**2 + (u[2]-v[2])**2))

a = -10
b = 10

for i in range(0, 10):
    v1 = [uniform(a, b), uniform(a, b), uniform(a, b)]
    v2 = [uniform(a, b), uniform(a, b), uniform(a, b)]
    v3 = [uniform(a, b), uniform(a, b), uniform(a, b)]
    v4 = [uniform(a, b), uniform(a, b), uniform(a, b)]
    v5 = [uniform(a, b), uniform(a, b), uniform(a, b)]
    v6 = [uniform(a, b), uniform(a, b), uniform(a, b)]
    lengths = {(1, 2) : dist(v1,v2), 
                (2, 3) : dist(v2,v3), 
                (3, 4) : dist(v3,v4), 
                (4, 5) : dist(v4,v5), 
                (5, 6) : dist(v5,v6), 
                (1, 6) : dist(v1,v6), 
                (1, 3) : dist(v1,v3), 
                (2, 4) : dist(v2,v4), 
                (3, 5) : dist(v3,v5), 
                (4, 6) : dist(v4,v6), 
                (1, 5) : dist(v1,v5), 
                (2, 6) : dist(v2,v6)}

    alg = AlgRealEmbeddings('Max6vertices', num_phi=5, num_theta=5,  name='random6vert')#, choice_from_clusters='closestToAverageLength')
#    alg.findMoreEmbeddings_tree(lengths, onlyOne=True)
    alg.findMoreEmbeddings(lengths)










