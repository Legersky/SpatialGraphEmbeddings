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
import time

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

# or the following edge lengths yields ones with 132 real embeddings:
lengths = {(2, 7): 10.536100000000047, (3, 2): 0.55, (2, 6): 1.002, (6, 8): 0.01656274447928313, (7, 8): 10.542352941009613, (6, 1): 2.0001, 
           (3, 1): 1.9342, (2, 8): 0.9868506797409413, (4, 7): 11.847100000000042, (2, 1): 1.9999, (5, 8): 4.4330622201743095, (4, 3): 5.4247000000000005, 
           (5, 1): 4.402399999999999, (5, 4): 7.0744, (3, 7): 10.524500000000048, (4, 1): 5.7963000000000005, (6, 5): 4.4449, (5, 7): 11.239600000000046}


start = time.time()
alg = AlgRealEmbeddings('Max8vertices', name='8vert_random_6and8close')
alg.findMoreEmbeddings(lengths)

print '\n\nSampling time: ',  time.time() - start
