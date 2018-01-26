'''
This script looks for edge lengths of G32a with 32 real embeddings by the sampling method based on coupler curves.
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

lengths = getEdgeLengthsByEmbedding('7vert32a', [v1, v2, v3, v4, v5, v6,v7])

alg = AlgRealEmbeddings('7vert32a',  name='7vert32a_random', num_phi=12, num_theta=12)
alg.findMoreEmbeddings(lengths)

#one of the obtained edge lengths:
#lengths = {(1, 2): 10.9464978614277, (4, 7): 8.7388471903030389, (1, 3): 6.27469003969716, (6, 7): 9.284894093430244, 
#           (5, 6): 9.229604229826325, (5, 7): 7.8831186688053965, (1, 4): 8.062543016998301, (1, 6): 11.559732906664415, 
#           (3, 6): 8.255957517599734, (2, 3): 8.827097289782333, (3, 7): 5.617334110160578, (2, 5): 9.741464052378111, 
#           (3, 4): 6.1089666048644204, (2, 4): 8.949700966931744, (3, 5): 5.5980614315661095}






