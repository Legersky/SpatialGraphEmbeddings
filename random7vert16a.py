from algRealEmbeddings import *
from graphEmbedding import *
from random import uniform

def dist( u, v):
    return float(np.sqrt( (u[0]-v[0])**2 + (u[1]-v[1])**2 + (u[2]-v[2])**2))

a = -5.0
b = 5.0
#for i in range(0, 100):
#    v1 = [uniform(a, b), uniform(a, b), uniform(a, b)]
#    v2 = [uniform(a, b), uniform(a, b), uniform(a, b)]
#    v3 = [uniform(a, b), uniform(a, b), uniform(a, b)]
#    v4 = [uniform(a, b), uniform(a, b), uniform(a, b)]
#    v5 = [uniform(a, b), uniform(a, b), uniform(a, b)]
#    v6 = [uniform(a, b), uniform(a, b), uniform(a, b)]
#    v7 = [uniform(a, b), uniform(a, b), uniform(a, b)]
#
#
#    lengths = {(1, 2) : dist(v1,v2),
#            (1, 3) : dist(v1,v3),
#            (1, 6) : dist(v1,v6),
#            (1, 7) : dist(v1,v7),
#            (2, 3) : dist(v2,v3),
#            (2, 4) : dist(v2,v4),
#            (2, 5) : dist(v2,v5),
#            (3, 4) : dist(v3,v4),
#            (3, 5) : dist(v3,v5),
#            (3, 6) : dist(v3,v6),
#            (3, 7) : dist(v3,v7),
#            (4, 6) : dist(v4,v6),
#            (4, 7) : dist(v4,v7),
#            (5, 6) : dist(v5,v6),
#            (5, 7) : dist(v5,v7),
#                }
#    G = GraphEmbedding(lengths, '7vert16a')
#    r = len(G.findEmbeddings()['real'])
#    print r
#    if r==16:
#        alg = AlgRealEmbeddings('7vert16a',  name='7vert16a_random')
#        alg.findMoreEmbeddings(lengths)


lengths = {(4, 7): 7.194704049075016, (1, 3): 5.747980801274709, (5, 6): 7.9040999968533, (1, 6): 8.475244632380354, (3, 7): 5.906043147314687, (2, 5): 7.154928516304775, (3, 5): 5.085705612775616, (1, 2): 4.358218168288518, (4, 6): 8.782535437545246, (5, 7): 10.216775546892137, (3, 6): 7.058175756732775, (1, 7): 3.771401840715499, (2, 3): 3.80939380782979, (3, 4): 3.227767798526635, (2, 4): 6.052913431118577}
for e in lengths:
    lengths[e] = round(lengths[e],2)
print lengths
G = GraphEmbedding(lengths, '7vert16a')
r = len(G.findEmbeddings()['real'])
print r


#from phcpy.solver import solve
#syst = G.getEquations()
#sol = solve(syst, verbose=0, tasks=2)

for e in lengths:
    print '(%d) to (%d) ' %e, 

