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
lengths = {(2, 7): 5.99642382396946, (1, 2): 11.051073596941777, (4, 7): 5.650901441338823, (2, 6): 5.701163774160325, (4, 6): 6.4865494481814689, (5, 6): 4.703787304280908, (5, 7): 5.774355494161055, (1, 4): 8.329818799504528, (1, 5): 9.395758665924081, (1, 3): 4.7651114601582645, (2, 3): 10.306127329218256, (3, 6): 8.569699108712607, (3, 7): 7.0971921318881925, (2, 5): 9.321205128818333, (3, 4): 7.642313190296686}
for e in lengths:
    lengths[e] = round(lengths[e],2)
print lengths

G = GraphEmbedding(lengths, '7vert24')
print len(G.findEmbeddings()['real'])

#alg = AlgRealEmbeddings('7vert24',  name='7vert24_random')
#alg.findMoreEmbeddings(lengths)


for e in lengths:
    print '(%d) to (%d) ' %e, 







