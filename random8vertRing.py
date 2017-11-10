from algRealEmbeddings import *
from graphEmbedding import *
from random import uniform

def dist( u, v):
    return float(np.sqrt( (u[0]-v[0])**2 + (u[1]-v[1])**2 + (u[2]-v[2])**2))

a = -0.05
b = 0.05

v1 = [0.6 + uniform(a, b), 0.6 + uniform(a, b), -3.0+ uniform(a, b)]

v2 = [-1.0+uniform(a, b), uniform(a, b), uniform(a, b)]
v4 = [-1.0+uniform(a, b), uniform(a, b), uniform(a, b)]
v6 = [-1.0+uniform(a, b), uniform(a, b), uniform(a, b)]

v3 = [1.0+uniform(a, b), uniform(a, b), uniform(a, b)]
v5 = [1.0+uniform(a, b), uniform(a, b), uniform(a, b)]
v7 = [1.0+uniform(a, b), uniform(a, b), uniform(a, b)]

v8 = [0.6 + uniform(a, b), 0.6 + uniform(a, b), 1.0+uniform(a, b)]



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


alg = AlgRealEmbeddings('Ring8vertices', name='8Ring_rnd_TBR')
alg.findMoreEmbeddings(lengths)










