'''
This script searches edge lengths for G48. 
Different sets of subgraphs are used for sampling.
The starting lengths were obtaind by parametric search, they have 28 real embeddings.
'''
print __doc__

import sys
sys.path.append('..')
from graphEmbeddings3D.algRealEmbeddings import AlgRealEmbeddings

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

alg_threeIterations = AlgRealEmbeddings('Max7vertices', num_phi=24, num_theta=24, factor_second=8, name='7vert_from_Vangelis28_threeIterations')
alg_threeIterations.findMoreEmbeddings(lengths, combinations=[[5, 6, 1, 7, 4], [4, 3, 1, 7, 5],[3, 2, 1, 7, 4]])
#result:
#{(1, 2): 1.99993774567597, (2, 7): 10.53609172287933, (4, 7): 11.847062363994983, (2, 6): 1.001987710974071, (6, 7): 10.53647884635266, 
#           (5, 6): 4.444937549609401, (5, 7): 11.239645246227646, (1, 4): 5.79633829024753, (1, 5): 4.402427328037957, (1, 3): 1.9341931539002355, 
#           (1, 6): 2.000134247468136, (4, 5): 7.07440359713745, (3, 7): 10.524474002128432, (3, 4): 5.424721797777055, (2, 3): 0.549983624664826}

alg_all_subgraphs = AlgRealEmbeddings('Max7vertices', num_phi=20, num_theta=20, factor_second=4, name='7vert_from_parametric28_all_subgraphs')
alg_all_subgraphs.findMoreEmbeddings(lengths)
#result:
#{(1, 2): 3.6547267239838073, (2, 7): 10.931902209424029, (4, 7): 10.528964771703047, (2, 6): 12.36535326089277, (6, 7): 6.614718312956418, 
#           (5, 6): 15.12587430560531, (5, 7): 13.990265768289897, (1, 4): 1.9423383670562493, (1, 5): 9.822922401885547, (1, 3): 1.9354608517198846,
#           (1, 6): 11.884171167050836, (4, 5): 9.723833100136652, (3, 7): 10.524273284504055, (3, 4): 0.5383946923429196, (2, 3): 3.283167602028824}

alg_ring_subgraphs = AlgRealEmbeddings('Max7vertices', num_phi=20, num_theta=20, factor_second=4, name='7vert_from_parametric28_ring_subgraphs')
alg_ring_subgraphs.findMoreEmbeddings(lengths,  combinations=[[2, 3, 1, 7, 6], [3, 4, 1, 7, 2], [4, 5, 1, 7, 3], [5, 6, 1, 7, 4], [6, 2, 1, 7, 5]],  allowed_repetition=2)
#result:
#{(1, 2): 4.1414653501173975, (2, 7): 11.139850811829968, (4, 7): 10.523551385333088, (2, 6): 4.1677512911570815, (6, 7): 10.542828736351707, 
#           (5, 6): 1.0309584430720742, (5, 7): 10.533375487651346, (1, 4): 1.9363529242001118, (1, 5): 1.9873340498152834, (1, 3): 5.959257111560759, 
#           (1, 6): 2.0323899508242214, (4, 5): 0.4412801229952831, (3, 7): 11.949398151358958, (3, 4): 5.6557115598418122, (2, 3): 7.295330049669027}













