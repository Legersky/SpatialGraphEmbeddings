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
        (4, 5) : dist(v4,v5),
        (4, 7) : dist(v4,v7),
        (5, 6) : dist(v5,v6),
        (6, 7) : dist(v6,v7),
            }

lengths = {(2, 7): 7.681242594724329, (1, 2): 11.063592708051631, (4, 7): 85.49235355033326, (2, 6): 7.109292101588832, (6, 7): 9.293565692425668, (4, 5): 78.53407268199916, (1, 4): 87.33452873149585, (1, 5): 21.489018832266296, (1, 3): 10.814254832812079, (2, 3): 4.469796047707818, (3, 6): 7.525153223773237, (5, 6): 22.077588664345829, (3, 7): 7.0981057752669106, (2, 5): 20.699169862228512, (3, 4): 84.16834214597186}
for e in lengths:
    lengths[e] = round(lengths[e],2)
print lengths

G = GraphEmbedding(lengths, '7vert32b')
print len(G.findEmbeddings()['real'])

#alg = AlgRealEmbeddings('7vert32b',  name='7vert32b_random')
#alg.findMoreEmbeddings(lengths)










