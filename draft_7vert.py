from algRealEmbeddings import *

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

#alg = AlgRealEmbeddings('Max7vertices', num_phi=20, num_theta=20, factor_second=4, name='7vert_from_Vangelis28_all_subgraphs')
#alg.findMoreEmbeddings(lengths)

alg = AlgRealEmbeddings('Max7vertices', num_phi=20, num_theta=20, factor_second=4, name='7vert_from_Vangelis28_only_GUI_subgraphs')
alg.findMoreEmbeddings(lengths,  combinations=[[2, 3, 1, 7, 6], [3, 4, 1, 7, 2], [4, 5, 1, 7, 3], [5, 6, 1, 7, 4], [6, 2, 1, 7, 5]],  allowed_repetition=2)


#alg = AlgRealEmbeddings('Max7vertices', name='7vert_from_Vangelis28_MoreDirect')
#alg.findMoreEmbeddings(lengths, num_phi=24, num_theta=24, factor_second=8, combinations=[[5, 6, 1, 7, 4], [4, 3, 1, 7, 5],[3, 2, 1, 7, 4]])















