'''
 Here, we try to find edge lengths of G48 with at least 36 real embeddings using edge lengths of g16 with many embeddings, and then, we use these lengths as the starting ones for G160.
'''
print __doc__

import sys
sys.path.append('..')
from graphEmbeddings3D.algRealEmbeddings import AlgRealEmbeddings
from graphEmbeddings3D.graphEmbedding import getEdgeLengthsByEmbedding, GraphEmbedding
from random import uniform

a = -.0
b = 1.0

for i in range(0, 100):
    v1 = [uniform(a, b), uniform(a, b), uniform(a, b)]
    v2 = [uniform(a, b), uniform(a, b), uniform(a, b)]
    v3 = [uniform(a, b), uniform(a, b), uniform(a, b)]
    v4 = [uniform(a, b), uniform(a, b), uniform(a, b)]
    v5 = [uniform(a, b), uniform(a, b), uniform(a, b)]
    v6 = [uniform(a, b), uniform(a, b), uniform(a, b)]

    lengths_6vert = getEdgeLengthsByEmbedding('Max6vertices', [v1, v2, v3, v4, v5, v6])
    alg6 = AlgRealEmbeddings('Max6vertices', num_phi=12, num_theta=12,  name='rnd_6vert')
    name6=alg6._fileNamePref
    lengths_6max = alg6.findMoreEmbeddings(lengths_6vert,  allowed_repetition=1, required_num=12)
        
    #---------------------creating 7 vertex lengths---------------------------------------------
    G = GraphEmbedding(lengths_6max, 'Max6vertices')
    v1, v2, v3, v7, v4, v5 = G.getEmbedding()
    v6 = [0, 0, 0]
    
    lengths_7vert = getEdgeLengthsByEmbedding('Max7vertices', [v1, v2, v3, v4, v5, v6, v7])

    alg = AlgRealEmbeddings('Max7vertices', name='7vert_from_6vert')
    name1st=alg._fileNamePref
    lengths_7 = alg.findMoreEmbeddings(lengths_7vert, required_num=36)

    G = GraphEmbedding(lengths_7, 'Max7vertices')
    v1, v2, v3, v4, v5, v6, v7 = G.getEmbedding()

    v8 = [0, 0, 0]

    lengths = getEdgeLengthsByEmbedding('Max8vertices', [v1, v2, v3, v4, v5, v6, v7, v8])

    alg2 = AlgRealEmbeddings('Max8vertices', name='from_rnd7to8')
    alg2.findMoreEmbeddings(lengths)
    del alg
    del alg2

