import sys
sys.path.append('..')
from graphEmbeddings3D.algRealEmbeddings import AlgRealEmbeddings
from graphEmbeddings3D.graphEmbedding import getEdgeLengthsByEmbedding
from random import uniform

print '''
This script looks for edge lengths of G32b with 32 real embeddings by the sampling method based on coupler curves.
'''
a = -5.0
b = 5.0

v1 = [uniform(a, b), uniform(a, b), uniform(a, b)]
v2 = [uniform(a, b), uniform(a, b), uniform(a, b)]
v3 = [uniform(a, b), uniform(a, b), uniform(a, b)]
v4 = [uniform(a, b), uniform(a, b), uniform(a, b)]
v5 = [uniform(a, b), uniform(a, b), uniform(a, b)]
v6 = [uniform(a, b), uniform(a, b), uniform(a, b)]
v7 = [uniform(a, b), uniform(a, b), uniform(a, b)]

lengths = getEdgeLengthsByEmbedding('7vert32b', [v1, v2, v3, v4, v5, v6,v7])

alg = AlgRealEmbeddings('7vert32b',  name='7vert32b_random', num_phi=12, num_theta=12)
alg.findMoreEmbeddings(lengths)

#one of the obtained edge lengths:
#lengths = {(2, 7): 7.681242594724329, (1, 2): 11.063592708051631, (4, 7): 85.49235355033326, (2, 6): 7.109292101588832, 
#            (6, 7): 9.293565692425668, (4, 5): 78.53407268199916, (1, 4): 87.33452873149585, (1, 5): 21.489018832266296,
#            (1, 3): 10.814254832812079, (2, 3): 4.469796047707818, (3, 6): 7.525153223773237, (5, 6): 22.077588664345829, 
#            (3, 7): 7.0981057752669106, (2, 5): 20.699169862228512, (3, 4): 84.16834214597186}






