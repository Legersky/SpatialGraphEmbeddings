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

lengths = {(1, 2): 6.845055297116067, (4, 7): 6.45129455236443, (1, 3): 7.811743249111052, (6, 7): 3.9135056088188374, (5, 6): 6.911895195279105, 
           (5, 7): 5.525090125221949, (1, 4): 7.049817468380529, (1, 6): 7.7694583789727103, (3, 6): 1.0228847255075912, (2, 3): 7.863694779624861,
           (3, 7): 4.400592844412562, (2, 5): 8.90115337290202, (3, 4): 5.666400953582232, (2, 4): 8.241835246296734, (3, 5): 7.4217088307340635}
G = GraphEmbedding(lengths, '7vert32a')
print len(G.findEmbeddings()['real'])

alg = AlgRealEmbeddings('7vert32a',  name='7vert32_a_random')
#    alg.findMoreEmbeddings_tree(lengths, onlyOne=True)
alg.findMoreEmbeddings(lengths)










