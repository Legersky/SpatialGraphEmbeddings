from algRealEmbeddings import *
from graphEmbedding import *
from random import uniform

def dist( u, v):
    return float(np.sqrt( (u[0]-v[0])**2 + (u[1]-v[1])**2 + (u[2]-v[2])**2))
#---------------------random 6 vertex---------------------------------------------
a = -10
b = 10
choice = 'closestToAverageLength'
for i in range(0, 10):
    v1 = [uniform(a, b), uniform(a, b), uniform(a, b)]
    v2 = [uniform(a, b), uniform(a, b), uniform(a, b)]
    v3 = [uniform(a, b), uniform(a, b), uniform(a, b)]
    v4 = [uniform(a, b), uniform(a, b), uniform(a, b)]
    v5 = [uniform(a, b), uniform(a, b), uniform(a, b)]
    v6 = [uniform(a, b), uniform(a, b), uniform(a, b)]
    lengths_6vert = {(1, 2) : dist(v1,v2), 
                (2, 3) : dist(v2,v3), 
                (3, 4) : dist(v3,v4), 
                (4, 5) : dist(v4,v5), 
                (5, 6) : dist(v5,v6), 
                (1, 6) : dist(v1,v6), 
                (1, 3) : dist(v1,v3), 
                (2, 4) : dist(v2,v4), 
                (3, 5) : dist(v3,v5), 
                (4, 6) : dist(v4,v6), 
                (1, 5) : dist(v1,v5), 
                (2, 6) : dist(v2,v6)}

    alg6 = AlgRealEmbeddings('Max6vertices', num_phi=12, num_theta=12, choice_from_clusters=choice,  name='rnd_6vert_for_7')
    name6=alg6._fileNamePref
    lengths_6max = alg6.findMoreEmbeddings(lengths_6vert,  allowed_repetition=1)
        
    #---------------------creating 7 vertex lengths---------------------------------------------
    G = GraphEmbedding(lengths_6max, 'Max6vertices')
    v1, v2, v3, v7, v4, v5 = G.getEmbedding()

    v6 = [0, 0, 0]

    lengths_1st_phase = {(1, 2) : dist(v1,v2),
               (1, 3) : dist(v1,v3),
               (1, 4) : dist(v1,v4),
               (1, 5) : dist(v1,v5),
               (1, 6) : dist(v1,v6),
               (7, 2) : dist(v7,v2),
               (7, 3) : dist(v7,v3),
               (7, 4) : dist(v7,v4),
               (7, 5) : dist(v7,v5),
               (7, 6) : dist(v7,v6),
                (2, 3) : dist(v2,v3), 
                (3, 4) : dist(v3,v4), 
                (4, 5) : dist(v4,v5), 
                (5, 6) : dist(v5,v6), 
                (2, 6) : dist(v2,v6), 
                }

    all_comb = [[2, 3, 1, 7, 6], [3, 7, 2, 4, 1], [5, 6, 1, 7, 4], [6, 1, 5, 2, 7], [3, 4, 1, 7, 2], [6, 7, 5, 2, 1], [4, 3, 1, 7, 5],
                                      [5, 1, 4, 6, 7], [2, 6, 1, 7, 3], [5, 7, 4, 6, 1], [3, 2, 1, 7, 4], [4, 1, 3, 5, 7], [6, 5, 1, 7, 2], [2, 7, 6, 3, 1],
                                      [3, 1, 2, 4, 7], [5, 4, 1, 7, 6], [4, 7, 3, 5, 1], [6, 2, 1, 7, 5], [2, 1, 6, 3, 7], [4, 5, 1, 7, 3]]
    alg = AlgRealEmbeddings('Max7vertices', choice_from_clusters=choice, name='7vert_from_6vert_1st_phase_'+name6)
    name1st=alg._fileNamePref
    lengths_2nd_phase = alg.findMoreEmbeddings(lengths_1st_phase, combinations=[comb for comb in all_comb if comb[0]==6], required_num=32)

    print '*********************2nd phase starts*****************************************'

    alg2 = AlgRealEmbeddings('Max7vertices', choice_from_clusters=choice, name='2nd_phase_'+name1st)
    alg2.findMoreEmbeddings(lengths_2nd_phase)
