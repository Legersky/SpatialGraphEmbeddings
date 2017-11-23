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
        (1, 6) : dist(v1,v6),
        (2, 3) : dist(v2,v3),
        (2, 4) : dist(v2,v4),
        (2, 5) : dist(v2,v5),
        (3, 4) : dist(v3,v4),
        (3, 5) : dist(v3,v5),
        (3, 6) : dist(v3,v6),
        (3, 7) : dist(v3,v7),
        (4, 7) : dist(v4,v7),
        (5, 6) : dist(v5,v6),
        (5, 7) : dist(v5,v7),
        (6, 7) : dist(v6,v7)
            }

G = GraphEmbedding(lengths, '7vert32a')
r = len(G.findEmbeddings()['real'])
print r
alg = AlgRealEmbeddings('7vert32a',  name='7vert32_a_random', num_phi=12, num_theta=12)
#    alg.findMoreEmbeddings_tree(lengths, onlyOne=True)
alg.findMoreEmbeddings(lengths)










