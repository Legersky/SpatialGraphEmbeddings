'''
This script looks for edge lengths of G24b with 24 real embeddings by the sampling method based on coupler curves.
'''
print __doc__

import sys
sys.path.append('..')
from graphEmbeddings3D.algRealEmbeddings import AlgRealEmbeddings
from graphEmbeddings3D.graphEmbedding import getEdgeLengthsByEmbedding
from random import uniform
 
a = -5.0
b = 5.0

v1 = [uniform(a, b), uniform(a, b), uniform(a, b)]
v2 = [uniform(a, b), uniform(a, b), uniform(a, b)]
v3 = [uniform(a, b), uniform(a, b), uniform(a, b)]
v4 = [uniform(a, b), uniform(a, b), uniform(a, b)]
v5 = [uniform(a, b), uniform(a, b), uniform(a, b)]
v6 = [uniform(a, b), uniform(a, b), uniform(a, b)]
v7 = [uniform(a, b), uniform(a, b), uniform(a, b)]

lengths = getEdgeLengthsByEmbedding('7vert24', [v1, v2, v3, v4, v5, v6,v7])

alg = AlgRealEmbeddings('7vert24',  name='7vert24_random', num_phi=12, num_theta=12)
alg.findMoreEmbeddings(lengths)

#one of the obtained edge lengths:
#lengths = {(2, 7): 5.99642382396946, (1, 2): 11.051073596941777, (4, 7): 5.650901441338823, (2, 6): 5.701163774160325, 
#           (4, 6): 6.4865494481814689, (5, 6): 4.703787304280908, (5, 7): 5.774355494161055, (1, 4): 8.329818799504528,
#           (1, 5): 9.395758665924081, (1, 3): 4.7651114601582645, (2, 3): 10.306127329218256, (3, 6): 8.569699108712607, 
#           (3, 7): 7.0971921318881925, (2, 5): 9.321205128818333, (3, 4): 7.642313190296686}






