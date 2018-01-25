import sys
sys.path.append('..')
from graphEmbeddings3D.algRealEmbeddings import AlgRealEmbeddings
from graphEmbeddings3D.graphEmbedding import getEdgeLengthsByEmbedding
from random import uniform

print '''
 This script uses the sampling method in order to find edge lengths of G160 with many real embeddings.
 Various kinds of random embeddings are used for obtaining starting edge lengths.
 The distance system is used for speed up, but it is not so robust.
'''
a = -0.1
b = 0.1
c = -10.0
d = 10.0

for i in range(0, 100):
    try:
        v1 = [uniform(a, b), uniform(a, b), uniform(a, b)]
        v8 = [uniform(a, b), uniform(a, b), uniform(a, b)]

        v2 = [1.0 + uniform(a, b), uniform(a, b), uniform(a, b)]
        v4 = [1.0 + uniform(a, b), uniform(a, b), uniform(a, b)]

        v6 = [uniform(a, b), 1.0 + uniform(a, b), uniform(a, b)]
        v7 = [uniform(a, b), 1.0 + uniform(a, b), uniform(a, b)]

        v3 = [uniform(a, b), uniform(a, b), 1.0 + uniform(a, b)]
        v5 = [uniform(a, b), uniform(a, b), 1.0 + uniform(a, b)]

        lengths = getEdgeLengthsByEmbedding('Max8vertices_distSyst', [v1, v2, v3, v4, v5, v6, v7, v8])

        alg = AlgRealEmbeddings('Max8vertices_distSyst', name='8vert_random_fromTetrahedron_dist')
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

        lengths = getEdgeLengthsByEmbedding('Max8vertices_distSyst', [v1, v2, v3, v4, v5, v6, v7, v8])
        
        alg = AlgRealEmbeddings('Max8vertices_distSyst', name='8vert_random_fromTwoDegreeCoincide_dist')
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

        lengths = getEdgeLengthsByEmbedding('Max8vertices_distSyst', [v1, v2, v3, v4, v5, v6, v7, v8])

        alg = AlgRealEmbeddings('Max8vertices_distSyst', name='8vert_random_dist')
        alg.findMoreEmbeddings(lengths)
    except:
        pass







