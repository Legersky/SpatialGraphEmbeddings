'''
This script looks for edge lengths of G48 with 48 real embeddings by the sampling method based on coupler curves.
The starting edge lengths are chosen randomly so that top, ring and bottom edges have similar lengths, respectively.
'''

print __doc__

import sys
sys.path.append('..')
from graphEmbeddings3D.algRealEmbeddings import AlgRealEmbeddings
from graphEmbeddings3D.graphEmbedding import getEdgeLengthsByEmbedding
from random import uniform
import time

a = -0.05
b = 0.05

v1 = [0.6 + uniform(a, b), 0.6 + uniform(a, b), -3.0+ uniform(a, b)]
v2 = [-1.0+uniform(a, b), uniform(a, b), uniform(a, b)]
v3 = [1.0+uniform(a, b), uniform(a, b), uniform(a, b)]
v4 = [-1.0+uniform(a, b), uniform(a, b), uniform(a, b)]
v5 = [1.0+uniform(a, b), uniform(a, b), uniform(a, b)]
v6 = [1.0+uniform(a, b), 1.0+uniform(a, b), uniform(a, b)]
v7 = [0.6 + uniform(a, b), 0.6 + uniform(a, b), 1.0+uniform(a, b)]


lengths = getEdgeLengthsByEmbedding('Max7vertices', [v1, v2, v3, v4, v5, v6,v7])

# or the following edge lengths yields ones with 48 real embeddings:
lengths = {(7, 3): 1.1607215105551865, (5, 6): 0.9637176545409848, (1, 3): 3.150534898739204, (4, 5): 2.052858811455116, (1, 4): 3.501511838512976,
           (7, 5): 1.2522096387483654, (1, 6): 3.065415486835149, (7, 2): 1.9742110125223091, (1, 2): 3.483745586129554, (7, 6): 1.1306270350127694,
           (2, 6): 2.1832117433322344, (7, 4): 1.9995809290124522, (1, 5): 3.062252806795317, (2, 3): 2.0335744682819774, (3, 4): 2.0640546146363996}

start = time.time()

alg = AlgRealEmbeddings('Max7vertices',  name='random7vert_tbr')
alg.findMoreEmbeddings(lengths)


print '\n\nSampling time: ',  time.time() - start


# possible result:
#lengths = {(1, 2): 3.2154982220265791, (2, 7): 1.4364348191352576, (4, 7): 17.879677295730822, (2, 6): 1.0237800945583886, (6, 7): 1.0148514465978793, 
#           (5, 6): 0.09495461781463814, (5, 7): 1.012349685743454, (1, 4): 18.137675008130483, (1, 5): 3.0323381652662267, (1, 3): 9.97279438841746, 
#           (1, 6): 3.0325698054329457, (4, 5): 17.849062036359275, (3, 7): 9.215048223037243, (3, 4): 19.765613007713743, (2, 3): 8.476438875510395}







