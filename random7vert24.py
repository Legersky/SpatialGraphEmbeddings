from algRealEmbeddings import *
from graphEmbedding import *
from random import uniform

def dist( u, v):
    return float(np.sqrt( (u[0]-v[0])**2 + (u[1]-v[1])**2 + (u[2]-v[2])**2))

a = -5.0
b = 5.0

v1 = [uniform(a, b), uniform(a, b), uniform(a, b)]
v2 = [uniform(a, b), uniform(a, b), uniform(a, b)]
v3 = [uniform(a, b), uniform(a, b), uniform(a, b)]
v4 = [uniform(a, b), uniform(a, b), uniform(a, b)]
v5 = [uniform(a, b), uniform(a, b), uniform(a, b)]
v6 = [uniform(a, b), uniform(a, b), uniform(a, b)]
v7 = [uniform(a, b), uniform(a, b), uniform(a, b)]


lengths = {(1, 2) : dist(v1,v2),
        (1, 3) : dist(v1,v3),
        (1, 4) : dist(v1,v4),
        (1, 5) : dist(v1,v5),
        (2, 3) : dist(v2,v3),
        (2, 5) : dist(v2,v5),
        (2, 6) : dist(v2,v6),
        (2, 7) : dist(v2,v7),
        (3, 4) : dist(v3,v4),
        (3, 6) : dist(v3,v6),
        (3, 7) : dist(v3,v7),
        (4, 6) : dist(v4,v6),
        (4, 7) : dist(v4,v7),
        (5, 6) : dist(v5,v6),
        (5, 7) : dist(v5,v7),
            }
G = GraphEmbedding(lengths, '7vert24')
print len(G.findEmbeddings()['real'])

alg = AlgRealEmbeddings('7vert24',  name='7vert24_random')
alg.findMoreEmbeddings(lengths)










