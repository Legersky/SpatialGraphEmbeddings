from algRealEmbeddings import *
from random import uniform

def dist( u, v):
    return float(np.sqrt( (u[0]-v[0])**2 + (u[1]-v[1])**2 + (u[2]-v[2])**2))

a = -0.2
b = 0.2

for i in range(0, 100):
    v1 = [0.6 + uniform(a, b), 0.6 + uniform(a, b), -3.0+ uniform(a, b)]
    v2 = [-1.0+uniform(a, b), uniform(a, b), uniform(a, b)]
    v3 = [1.0+uniform(a, b), uniform(a, b), uniform(a, b)]
    v4 = [-1.0+uniform(a, b), uniform(a, b), uniform(a, b)]
    v5 = [1.0+uniform(a, b), uniform(a, b), uniform(a, b)]
    v6 = [1.0+uniform(a, b), 1.0+uniform(a, b), uniform(a, b)]
    v7 = [0.6 + uniform(a, b), 0.6 + uniform(a, b), 1.0+uniform(a, b)]

    lengths = {(1, 2) : dist(v1,v2),
               (1, 3) : dist(v1,v3),
               (1, 4) : dist(v1,v4),
               (1, 5) : dist(v1,v5),
               (1, 6) : dist(v1,v6),
               (2, 7) : dist(v7,v2),
               ( 3, 7) : dist(v7,v3),
               ( 4, 7) : dist(v7,v4),
               ( 5, 7) : dist(v7,v5),
               ( 6, 7) : dist(v7,v6),
                (2, 3) : dist(v2,v3), 
                (3, 4) : dist(v3,v4), 
                (4, 5) : dist(v4,v5), 
                (5, 6) : dist(v5,v6), 
                (2, 6) : dist(v2,v6), 
                }
    alg = AlgRealEmbeddings('Max7vertices',  name='random7vert')
    lengths_7 = alg.findMoreEmbeddings(lengths, required_num=36)

    G = GraphEmbedding(lengths_7, 'Max7vertices')
    v1, v2, v3, v4, v5, v6, v7 = G.getEmbedding()

    v8 = [0, 0, 0]

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

    alg2 = AlgRealEmbeddings('Max8vertices', name='from7to8from_rnd')
    alg2.findMoreEmbeddings(lengths)
    del alg
    del alg2

