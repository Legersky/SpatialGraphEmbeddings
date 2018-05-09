'''
This script searches edge lengths for G48 by tree approach, namely, all combinations of subgraphs are tested.
The starting lengths were obtaind by parametric search, they have 28 real embeddings.
'''
print __doc__

import sys
sys.path.append('..')
from graphEmbeddings3D.algRealEmbeddings import AlgRealEmbeddings
import time
start = time.time()
lengths = {'12': 1.99993774567597,
         '13': 1.99476987780024,
         '14': 2.003436460984393,
         '15': 2.00289249524296,
         '16': 2.000134247468136,
         '23': 0.999614322089483,
         '26': 1.001987710974071,
         '27': 10.53609172287933,
         '34': 1.00368644488060,
         '37': 10.53631716364608,
         '45': 1.001530148504854,
         '47': 10.53572330314948,
         '56': 0.995723616535744,
         '57': 10.53627365999783,
         '67': 10.53647884635266}
alg = AlgRealEmbeddings('Max7vertices',  name='par7vert_tree_search')
alg.findMoreEmbeddings_tree(lengths, onlyOne=True)

print '\n\nSampling time: ',  time.time() - start




