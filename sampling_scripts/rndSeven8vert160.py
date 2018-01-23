from graphEmbeddings3D.algRealEmbeddings import AlgRealEmbeddings
from graphEmbeddings3D.graphEmbedding import getEdgeLengthsByEmbedding, GraphEmbedding
from random import uniform

print '''
 Here, we try to find edge elnegths of G48 with at least 36 real embeddings, and then, we use these lengths as the starting ones for G160.
'''
a = -0.2
b = 0.2

for i in range(0, 100):
    v1 = [0.6 + uniform(a, b), 0.6 + uniform(a, b), -3.0+ uniform(a, b)]
    v2 = [-1.0+uniform(a, b), uniform(a, b), uniform(a, b)]
    v3 = [1.0+uniform(a, b), uniform(a, b), uniform(a, b)]
    v4 = [-1.0+uniform(a, b), uniform(a, b), uniform(a, b)]
    v5 = [1.0+uniform(a, b), uniform(a, b), uniform(a, b)]
    v6 = [1.0+uniform(a, b), 1.0+uniform(a, b), uniform(a, b)]
    v7 = [0.6 + uniform(a, b), 0.6 + uniform(a, b), 1.0+uniform(a, b)]

    lengths = getEdgeLengthsByEmbedding('Max7vertices', [v1, v2, v3, v4, v5, v6, v7])
    
    alg = AlgRealEmbeddings('Max7vertices',  name='random7vert')
    lengths_7 = alg.findMoreEmbeddings(lengths, required_num=36)

    G = GraphEmbedding(lengths_7, 'Max7vertices')
    v1, v2, v3, v4, v5, v6, v7 = G.getEmbedding()

    v8 = [0, 0, 0]

    lengths = getEdgeLengthsByEmbedding('Max8vertices', [v1, v2, v3, v4, v5, v6, v7, v8])

    alg2 = AlgRealEmbeddings('Max8vertices', name='from_rnd7to8')
    alg2.findMoreEmbeddings(lengths)
    del alg
    del alg2

