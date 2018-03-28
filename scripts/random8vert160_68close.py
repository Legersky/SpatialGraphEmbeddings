'''
 This script uses the sampling method in order to find edge lengths of G160 with many real embeddings.
 Various kinds of random embeddings are used for obtaining starting edge lengths.
'''
print __doc__

import sys
sys.path.append('..')
from graphEmbeddings3D.algRealEmbeddings import AlgRealEmbeddings
from graphEmbeddings3D.graphEmbedding import getEdgeLengthsByEmbedding
from random import uniform

a = -0.05
b = 0.05

v1 = [0.6 + uniform(a, b), 0.6 + uniform(a, b), -3.0+ uniform(a, b)]

v2 = [-1.0+uniform(a, b), uniform(a, b), uniform(a, b)]

v3 = [1.0+uniform(a, b), uniform(a, b), uniform(a, b)]

v4 = [-1.0+uniform(a, b), uniform(a, b), uniform(a, b)]

v5 = [1.0+uniform(a, b), uniform(a, b), uniform(a, b)]

v6 = [1.0+uniform(a, b), 1.0+uniform(a, b), uniform(a, b)]
v8 = [1.0+uniform(a, b), 1.0+uniform(a, b), uniform(a, b)]

v7 = [0.6 + uniform(a, b), 0.6 + uniform(a, b), 1.0+uniform(a, b)]

lengths = getEdgeLengthsByEmbedding('Max8vertices', [v1, v2, v3, v4, v5, v6, v7, v8])

alg = AlgRealEmbeddings('Max8vertices', name='8vert_random_6and8close')
alg.findMoreEmbeddings(lengths)
