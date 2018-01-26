'''
This script looks for edge lengths of G16b with 16 real embeddings by the sampling method based on coupler curves.
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

lengths = getEdgeLengthsByEmbedding('7vert16b', [v1, v2, v3, v4, v5, v6,v7])

alg = AlgRealEmbeddings('7vert16b',  name='7vert16b_random')
alg.findMoreEmbeddings(lengths)

#one of the obtained edge lengths:
#lengths = {(2, 7): 5.895983148228586, (1, 2): 4.6193863752177702, (4, 7): 4.461504253664427, (2, 6): 7.47123257104987, 
#           (6, 7): 3.0880879139441535, (4, 6): 7.067863188235056, (4, 5): 7.716042577713717, (1, 4): 6.514961476231782,
#           (1, 5): 5.692420153393312, (1, 3): 3.526562220615557, (2, 3): 7.686019585902205, (3, 6): 6.432943048336337, 
#           (3, 7): 5.762187187187983, (2, 5): 9.484862380069352, (3, 5): 6.097065470256926}
