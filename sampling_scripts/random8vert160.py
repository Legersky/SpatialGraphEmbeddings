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

a = -0.1
b = 0.1
c = -10.0
d = 10.0

try:
    v1 = [uniform(a, b), uniform(a, b), uniform(a, b)]
    v8 = [uniform(a, b), uniform(a, b), uniform(a, b)]

    v2 = [1.0 + uniform(a, b), uniform(a, b), uniform(a, b)]
    v4 = [1.0 + uniform(a, b), uniform(a, b), uniform(a, b)]

    v6 = [uniform(a, b), 1.0 + uniform(a, b), uniform(a, b)]
    v7 = [uniform(a, b), 1.0 + uniform(a, b), uniform(a, b)]

    v3 = [uniform(a, b), uniform(a, b), 1.0 + uniform(a, b)]
    v5 = [uniform(a, b), uniform(a, b), 1.0 + uniform(a, b)]

    lengths = getEdgeLengthsByEmbedding('Max8vertices', [v1, v2, v3, v4, v5, v6, v7, v8])

    alg = AlgRealEmbeddings('Max8vertices', name='8vert_random_fromTetrahedron',  allowedNumberOfMissing=14)
    alg.findMoreEmbeddings(lengths)
except:
    pass
#   ---------------------------------
try:
    v1 = [uniform(c, d), uniform(c, d), uniform(c, d)]
    v8 = [uniform(c, d), uniform(c, d), uniform(c, d)]

    v2 = [uniform(a, b), uniform(a, b), uniform(a, b)]
    v5 = [uniform(a, b), uniform(a, b), uniform(a, b)]

    v6 = [uniform(c, d), uniform(c, d), uniform(c, d)]
    v7 = [uniform(c, d), uniform(c, d), uniform(c, d)]

    v3 = [uniform(c, d), uniform(c, d), uniform(c, d)]
    v4 = [uniform(c, d), uniform(c, d), uniform(c, d)]

    lengths = getEdgeLengthsByEmbedding('Max8vertices', [v1, v2, v3, v4, v5, v6, v7, v8])

    alg = AlgRealEmbeddings('Max8vertices', name='8vert_random_fromTwoDegreeCoincide')
    alg.findMoreEmbeddings(lengths)
except:
    pass
#   ---------------------------------
try:
    v1 = [uniform(c, d), uniform(c, d), uniform(c, d)]
    v8 = [uniform(c, d), uniform(c, d), uniform(c, d)]

    v2 = [uniform(c, d), uniform(c, d), uniform(c, d)]
    v5 = [uniform(c, d), uniform(c, d), uniform(c, d)]

    v6 = [uniform(c, d), uniform(c, d), uniform(c, d)]
    v7 = [uniform(c, d), uniform(c, d), uniform(c, d)]

    v3 = [uniform(c, d), uniform(c, d), uniform(c, d)]
    v4 = [uniform(c, d), uniform(c, d), uniform(c, d)]

    lengths = getEdgeLengthsByEmbedding('Max8vertices', [v1, v2, v3, v4, v5, v6, v7, v8])

    alg = AlgRealEmbeddings('Max8vertices', name='8vert_random')
    alg.findMoreEmbeddings(lengths)
except:
    pass







