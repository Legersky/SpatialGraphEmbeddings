'''
This script looks for edge lengths of G16a with 16 real embeddings by random guessing.
'''
print __doc__

import sys
sys.path.append('..')
from graphEmbeddings3D.algRealEmbeddings import AlgRealEmbeddings
from graphEmbeddings3D.graphEmbedding import getEdgeLengthsByEmbedding,  GraphEmbedding
from random import uniform



a = -5.0
b = 5.0
for i in range(0, 100):
    v1 = [uniform(a, b), uniform(a, b), uniform(a, b)]
    v2 = [uniform(a, b), uniform(a, b), uniform(a, b)]
    v3 = [uniform(a, b), uniform(a, b), uniform(a, b)]
    v4 = [uniform(a, b), uniform(a, b), uniform(a, b)]
    v5 = [uniform(a, b), uniform(a, b), uniform(a, b)]
    v6 = [uniform(a, b), uniform(a, b), uniform(a, b)]
    v7 = [uniform(a, b), uniform(a, b), uniform(a, b)]


    lengths = getEdgeLengthsByEmbedding('7vert16a', [v1, v2, v3, v4, v5, v6,v7])
    G = GraphEmbedding(lengths, '7vert16a')
    r = len(G.findEmbeddings()['real'])
    print r
    if r==16:
        alg = AlgRealEmbeddings('7vert16a',  name='7vert16a_random')
        alg.findMoreEmbeddings(lengths)
        break

#one of the obtained edge lengths:
#lengths = {(4, 7): 7.194704049075016, (1, 3): 5.747980801274709, (5, 6): 7.9040999968533, (1, 6): 8.475244632380354, 
#           (3, 7): 5.906043147314687, (2, 5): 7.154928516304775, (3, 5): 5.085705612775616, (1, 2): 4.358218168288518, 
#           (4, 6): 8.782535437545246, (5, 7): 10.216775546892137, (3, 6): 7.058175756732775, (1, 7): 3.771401840715499, 
#           (2, 3): 3.80939380782979, (3, 4): 3.227767798526635, (2, 4): 6.052913431118577}


