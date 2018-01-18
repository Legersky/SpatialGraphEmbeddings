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
lengths = {(1, 2): 10.9464978614277, (4, 7): 8.7388471903030389, (1, 3): 6.27469003969716, (6, 7): 9.284894093430244, (5, 6): 9.229604229826325, (5, 7): 7.8831186688053965, (1, 4): 8.062543016998301, (1, 6): 11.559732906664415, (3, 6): 8.255957517599734, (2, 3): 8.827097289782333, (3, 7): 5.617334110160578, (2, 5): 9.741464052378111, (3, 4): 6.1089666048644204, (2, 4): 8.949700966931744, (3, 5): 5.5980614315661095}
for e in lengths:
    lengths[e] = round(lengths[e],2)
print lengths
G = GraphEmbedding(lengths, '7vert32a')
r = len(G.findEmbeddings()['real'])
print r
#alg = AlgRealEmbeddings('7vert32a',  name='7vert32_a_random', num_phi=12, num_theta=12)
##    alg.findMoreEmbeddings_tree(lengths, onlyOne=True)
#alg.findMoreEmbeddings(lengths)

for e in lengths:
    print '(%d) to (%d) ' %e, 









