'''
 This script uses the sampling method in order to find edge lengths of G160 with many real embeddings.
 Starting point G_48-48 embeddings
(6 vertex (almost) coincides with 5):

w1 = [0, 0, 0]
w2 = [0, 5.73225265711452, 0]
w3 = [0, 0.148291560665423, 3.30491899417522]
w4 = [-3.1508917850189, 0.846961251730797, 0.259525376012296]
w5 = [-2.49285516161801, 1.67204766459525, -0.558220040757827]
w6 = [-2.59381676560345, 1.54400553897266, -0.398999543683517]
w7 = [-1.69450165498139, 1.01237857741641, -0.346180158496239]
Permutation from G48 {1: 3, 2: 1, 3: 2, 4: 5, 5: 7, 6: 4, 7: 8}+{v6=v5+pert} .
 

 
'''
print __doc__

import sys
sys.path.append('..')
from graphEmbeddings3D.algRealEmbeddings import AlgRealEmbeddings
from graphEmbeddings3D.graphEmbedding import getEdgeLengthsByEmbedding,  GraphEmbedding
from random import uniform
#48 embedding for G48 
a = -0.05
b = 0.05
c = -10.0
d = 10.0
edges_G160=[(1, 2), (1, 3), (1, 4), (1, 8), (2, 3), (2, 5), (2, 6), (2, 8), (3, 4), (3, 5), (3, 7),(4, 7), (4, 8), (5, 6), (5, 7), (6, 7), (6, 8), (7, 8)]
#try:
v3 = [0, 0, 0]
v1 = [0, 5.73225265711452, 0]
v2 = [0, 0.148291560665423, 3.30491899417522]

v5 = [-3.1508917850189, 0.846961251730797, 0.259525376012296]
v7 = [-2.49285516161801, 1.67204766459525, -0.558220040757827]
v4 = [-2.59381676560345, 1.54400553897266, -0.398999543683517]
v8 = [-1.69450165498139, 1.01237857741641, -0.346180158496239]

v6 = [-3.1508917850189+uniform(a,b), 0.846961251730797+uniform(a,b), 0.259525376012296+uniform(a,b)]


lengths = getEdgeLengthsByEmbedding('edges',  [v1, v2, v3, v4, v5, v6, v7, v8], edges=edges_G160)

alg = AlgRealEmbeddings('edges', name='G160',  edges=edges_G160,  num_sols=160, allowedNumberOfMissing=4, moreSamplingSubgraphs=False)
lengths_final =alg.findMoreEmbeddings(lengths)
G = GraphEmbedding(lengths_final, 'edges', num_sols=160)
print 'There are ', len(G.findEmbeddings()['real']), 'real embeddings with the following edge lengths:'
print lengths_final
#except:
#    pass
#   ---------------------------------
