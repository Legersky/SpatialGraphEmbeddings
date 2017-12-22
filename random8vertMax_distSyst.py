from algRealEmbeddings import *
from graphEmbedding import *
from random import uniform

a = -0.1
b = 0.1
c = -10.0
d = 10.0

def dist( u, v):
        return float(np.sqrt( (u[0]-v[0])**2 + (u[1]-v[1])**2 + (u[2]-v[2])**2))

for i in range(0, 100):
    try:
        v1 = [uniform(a, b), uniform(a, b), uniform(a, b)]
        v8 = [uniform(a, b), uniform(a, b), uniform(a, b)]

        v2 = [1.0 + uniform(a, b), uniform(a, b), uniform(a, b)]
        v4 = [1.0 + uniform(a, b), uniform(a, b), uniform(a, b)]

        v6 = [uniform(a, b), 1.0 + uniform(a, b), uniform(a, b)]
        v7 = [uniform(a, b), 1.0 + uniform(a, b), uniform(a, b)]

        v3 = [uniform(a, b), uniform(a, b), 1.0 + uniform(a, b)]
        v5 = [uniform(a, b), uniform(a, b), 1.0 + uniform(a, b)]

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


        alg = AlgRealEmbeddings('Max8vertices_distSyst', name='8vert_random_fromTetrahedron_dist')
        alg.findMoreEmbeddings(lengths)
    except:
        pass
#   ---------------------------------
    try:
        v1 = [uniform(c, d), uniform(c, d), uniform(c, d)]
        v8 = [uniform(c, d), uniform(c, d), uniform(c, d)]

        v2 = [uniform(a, b), uniform(a, b), uniform(a, b)]
        v5 = [uniform(a, b), uniform(a, b), uniform(a, b)]

        v6 = [uniform(c, d), uniform(c, d), uniform(c, d)]
        v7 = [uniform(c, d), uniform(c, d), uniform(c, d)]

        v3 = [uniform(c, d), uniform(c, d), uniform(c, d)]
        v4 = [uniform(c, d), uniform(c, d), uniform(c, d)]

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


        alg = AlgRealEmbeddings('Max8vertices_distSyst', name='8vert_random_fromTwoDegreeCoincide_dist')
        alg.findMoreEmbeddings(lengths)
    except:
        pass
#   ---------------------------------
    try:
        v1 = [uniform(c, d), uniform(c, d), uniform(c, d)]
        v8 = [uniform(c, d), uniform(c, d), uniform(c, d)]

        v2 = [uniform(c, d), uniform(c, d), uniform(c, d)]
        v5 = [uniform(c, d), uniform(c, d), uniform(c, d)]

        v6 = [uniform(c, d), uniform(c, d), uniform(c, d)]
        v7 = [uniform(c, d), uniform(c, d), uniform(c, d)]

        v3 = [uniform(c, d), uniform(c, d), uniform(c, d)]
        v4 = [uniform(c, d), uniform(c, d), uniform(c, d)]

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


        alg = AlgRealEmbeddings('Max8vertices_distSyst', name='8vert_random_dist')
        alg.findMoreEmbeddings(lengths)
    except:
        pass







