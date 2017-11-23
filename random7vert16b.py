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
        (3, 5) : dist(v3,v5),
        (3, 6) : dist(v3,v6),
        (3, 7) : dist(v3,v7),
        (4, 5) : dist(v4,v5),
        (4, 6) : dist(v4,v6),
        (4, 7) : dist(v4,v7),
        (6, 7) : dist(v6,v7),
            }
lengths = {(2, 7): 5.895983148228586, (1, 2): 4.6193863752177702, (4, 7): 4.461504253664427, (2, 6): 7.47123257104987, (6, 7): 3.0880879139441535, (4, 6): 7.067863188235056, (4, 5): 7.716042577713717, (1, 4): 6.514961476231782, (1, 5): 5.692420153393312, (1, 3): 3.526562220615557, (2, 3): 7.686019585902205, (3, 6): 6.432943048336337, (3, 7): 5.762187187187983, (2, 5): 9.484862380069352, (3, 5): 6.097065470256926}
for e in lengths:
    lengths[e] = round(lengths[e],2)
print lengths
G = GraphEmbedding(lengths, '7vert16b')
print len(G.findEmbeddings()['real'])

#alg = AlgRealEmbeddings('7vert16b',  name='7vert16b_random')
#alg.findMoreEmbeddings(lengths)
