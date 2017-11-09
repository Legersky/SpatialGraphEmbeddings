from algRealEmbeddings import *

#lengths_7 = {(1, 2): 5.7322526571145165, (2, 7): 5.026767140924392, (4, 7): 1.5859742422378464, (2, 6): 4.942519597534114, (6, 7): 1.0460328780197365, 
#             (5, 6): 0.22790172904533793, (5, 7): 1.0571151960916287, (1, 4): 3.273043816331718, (1, 5): 3.0531426208724093, (1, 3): 3.308244239022385, 
#             (1, 6): 3.044838116183432, (4, 5): 1.3350982559642988, (3, 7): 4.116856477203278, (3, 4): 4.437418243490294, (2, 3): 6.488691014736114}
##             /home/jan/Research/MaximumReal3dEmbeddings/CouplerCurve/important/7vert/48_embd_from_16_random_TBR.txt
lengths_7 = {(1, 2): 6.274873914397456, (2, 7): 12.093170928755793, (4, 7): 10.52881316000657, (2, 6): 10.479941653322763, (6, 7): 8.11321631042846, 
             (5, 6): 9.459324566316262, (5, 7): 10.73159467196586, (1, 4): 1.9761239449935557, (1, 5): 2.8483166740669206, (1, 3): 1.9564840951854308, 
             (1, 6): 9.4046159048512159, (4, 5): 2.313950594269001, (3, 7): 10.525064001461296, (3, 4): 0.47038007668728976, (2, 3): 5.952813550679546}
#/home/jan/Research/MaximumReal3dEmbeddings/CouplerCurve/important/7vert/48_embd_from_28_7vert_53129440.txt
G = GraphEmbedding(lengths_7, 'Max7vertices')
v1, v2, v3, v4, v5, v6, v7 = G.getEmbedding()

def dist( u, v):
    return float(np.sqrt( (u[0]-v[0])**2 + (u[1]-v[1])**2 + (u[2]-v[2])**2))

v8 = [0, 0, 0]

lengths = {
           (2, 1) : dist(v2,v1),
            (2, 7) : dist(v2,v7),
            (2, 6) : dist(v2,v6),
            (3, 1) : dist(v3,v1),
            (3, 7) : dist(v3,v7),
            (3, 2) : dist(v3,v2),
            (4, 1) : dist(v4,v1),
            (4, 7) : dist(v4,v7),
            (4, 3) : dist(v4,v3),
            (5, 1) : dist(v5,v1),
            (5, 7) : dist(v5,v7),
            (5, 4) : dist(v5,v4),
            (6, 1) : dist(v6,v1),
            (6, 5) : dist(v6,v5),
            (5, 8) : dist(v5,v8),
            (6, 8) : dist(v6,v8),
            (7, 8) : dist(v7,v8),
            (2, 8) : dist(v2,v8),
           }

all_comb = [[3, 7, 2, 4, 1], [6, 8, 2, 5, 1], [4, 5, 1, 7, 3], [8, 7, 2, 5, 6], [4, 3, 1, 7, 5], [6, 2, 1, 8, 5], [3, 2, 1, 7, 4], [8, 5, 7, 6, 2], 
                                  [4, 1, 3, 5, 7], [6, 5, 1, 8, 2], [3, 1, 2, 4, 7], [8, 2, 7, 6, 5], [4, 7, 3, 5, 1], [6, 1, 2, 5, 8], [3, 4, 1, 7, 2], [8, 6, 2, 5, 7]]
alg = AlgRealEmbeddings('Max8vertices', name='8vert_from_7vert')
#alg.findMoreEmbeddings_tree(lengths,  combinations=[comb for comb in all_comb if comb[0]==6])
alg.findMoreEmbeddings(lengths, combinations=[comb for comb in all_comb if comb[0]==8])
