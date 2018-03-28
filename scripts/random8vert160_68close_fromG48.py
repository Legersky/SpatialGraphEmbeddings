'''
 This script uses the sampling method in order to find edge lengths of G160 with many real embeddings.
 Various kinds of random embeddings are used for obtaining starting edge lengths.
'''
print __doc__

import sys
sys.path.append('..')
from graphEmbeddings3D.algRealEmbeddings import AlgRealEmbeddings
from graphEmbeddings3D.graphEmbedding import getEdgeLengthsByEmbedding, GraphEmbedding
from random import uniform

a = -0.05
b = 0.05

lengths_7 = {(2, 7): 10.5361, (4, 7): 11.8471, (2, 6): 1.002, (5, 6): 4.4449, (1, 4): 5.7963, (1, 3): 1.9342, (1, 6): 2.0001, (3, 7): 10.5245, 
             (1, 2): 1.9999, (6, 7): 10.5365, (5, 7): 11.2396, (1, 5): 4.4024, (4, 5): 7.0744, (2, 3): 0.55, (3, 4): 5.4247}

G = GraphEmbedding(lengths_7, 'Max7vertices')
v1, v2, v3, v4, v5, v6, v7 = G.getEmbedding()

v8 = [v6[0] + uniform(a, b), v6[1] + uniform(a, b), v6[2] + uniform(a, b)]

lengths = getEdgeLengthsByEmbedding('Max8vertices', [v1, v2, v3, v4, v5, v6, v7, v8])

alg = AlgRealEmbeddings('Max8vertices', name='8vert_6and8close_fromG48')
alg.findMoreEmbeddings(lengths)
