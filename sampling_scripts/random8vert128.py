'''
This script looks for edge lengths of G128 with many real embeddings by the sampling method based on coupler curves.
The starting edge lengths are chosen randomly so that top, ring and bottom edges have similar lengths, respectively.
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
v4 = [-1.0+uniform(a, b), uniform(a, b), uniform(a, b)]
v6 = [-1.0+uniform(a, b), uniform(a, b), uniform(a, b)]

v3 = [1.0+uniform(a, b), uniform(a, b), uniform(a, b)]
v5 = [1.0+uniform(a, b), uniform(a, b), uniform(a, b)]
v7 = [1.0+uniform(a, b), uniform(a, b), uniform(a, b)]

v8 = [0.6 + uniform(a, b), 0.6 + uniform(a, b), 1.0+uniform(a, b)]



lengths = getEdgeLengthsByEmbedding('Ring8vertices', [v1, v2, v3, v4, v5, v6, v7, v8])

alg = AlgRealEmbeddings('Ring8vertices', name='8Ring_rnd_TBR')
alg.findMoreEmbeddings(lengths)










