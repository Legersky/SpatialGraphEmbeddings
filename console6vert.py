from algRealEmbeddings import *
from graphEmbedding import *
from phcpy.solver import solve, mixed_volume
from phcpy.trackers import track, gpu_double_track
from phcpy.solutions import strsol2dict, is_real

v1, v2, v3, v4, v5, v6 =[
                        [0, 1.50000000000000, 0],
                        [0, 1.00000000000000, 0],
                        [-0.215583789828996, -0.994162101710040, 0.621965547009129],
                        [-1.68361243153425, 1.34617354563865, 0.166201392119758],
                        [0.700033202689246, 0.110676360238974, 3.34383891004353],
                        [1.00000000000000, 1.10000000000000, 0]
                         ]

def dist( u, v):
    return float(np.sqrt( (u[0]-v[0])**2 + (u[1]-v[1])**2 + (u[2]-v[2])**2))
    
lengths = {(1, 2) : dist(v1,v2), 
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

#G = GraphEmbedding(lengths, 'Max6vertices')
##sols = G.findEmbeddings()
##print len(sols['real'])
#
#
#syst = G.getEquations()
#
#start = time.time()
#
#sols = solve(syst,tasks=2)
#result_real = []
#result_complex = []
#
#for sol in sols:
#    soldic = strsol2dict(sol)
#    if is_real(sol, 1.0e-8):
#        result_real.append(soldic)
#    else:
#        result_complex.append(soldic)
#            
#end = time.time()
#print len(result_real)
#
#print 'time: ',(end - start)


possibleParametrizedVertices = {'[2, 1, 6, 3, 4]': [2, 1, 6, 3, 4],
                             '[2, 3, 1, 4, 6]': [2, 3, 1, 4, 6],
                             '[2, 4, 6, 3, 1]': [2, 4, 6, 3, 1],
                             '[2, 6, 1, 4, 3]': [2, 6, 1, 4, 3],
                             '[3, 1, 2, 5, 4]': [3, 1, 2, 5, 4],
                             '[3, 2, 1, 4, 5]': [3, 2, 1, 4, 5],
                             '[3, 4, 2, 5, 1]': [3, 4, 2, 5, 1],
                             '[3, 5, 1, 4, 2]': [3, 5, 1, 4, 2],
                             '[5, 1, 3, 6, 4]': [5, 1, 3, 6, 4],
                             '[5, 3, 1, 4, 6]': [5, 3, 1, 4, 6],
                             '[5, 4, 3, 6, 1]': [5, 4, 3, 6, 1],
                             '[5, 6, 1, 4, 3]': [5, 6, 1, 4, 3],
                             '[6, 1, 5, 2, 4]': [6, 1, 5, 2, 4],
                             '[6, 2, 1, 4, 5]': [6, 2, 1, 4, 5],
                             '[6, 4, 5, 2, 1]': [6, 4, 5, 2, 1],
                             '[6, 5, 1, 4, 2]': [6, 5, 1, 4, 2]}




alg = AlgRealEmbeddings(lengths, 'Max6vertices')
alg.findMoreEmbeddings(lengths, 
                       20, 20, 
                       possibleParametrizedVertices.values(), 
                       0, 
                       name='6vert')





















