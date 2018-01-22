from algRealEmbeddings import *
from graphEmbedding import *
#from random import uniform

lengths = {(1, 2): 14.675900223725295, (4, 7): 38.76311928178239, (2, 6): 1.6158929752494409, (4, 5): 70.34167136306003, (2, 8): 1.6098988097760762,
           (1, 4): 39.647349233450726, (1, 3): 13.926417286907368, (1, 6): 14.590181277659985, (3, 7): 8.445778625611991, (5, 8): 57.10907857052827,
           (2, 7): 10.697280470178171, (6, 8): 0.3644390005351665, (1, 5): 58.63814934691712, (7, 8): 10.594204813377264, (5, 7): 57.916102578212076, 
           (5, 6): 57.11505459632309, (2, 3): 10.46486347500465, (3, 4): 38.981201549967366} #128 embeddings


for e in lengths:
    lengths[e] = round(lengths[e],2)
print lengths

G = GraphEmbedding(lengths, 'Max8vertices')
r = len(G.findEmbeddings(tolerance=1.0e-80)['real'])
print r


#
#print '\nDistance system'
#G = GraphEmbedding(lengths, 'Max8vertices_distSyst')
#sols = G.findEmbeddings()
#print '# embeddable (dist syst):'
#print len(sols['real'])


#from phcpy.solver import solve
#syst = G.getEquations()
#sol = solve(syst, verbose=0, tasks=2)

for e in lengths:
    print '(%d) to (%d) ' %e, 
    
    
alg = AlgRealEmbeddings('Max8vertices', name='8vert_from_128_')
alg.findMoreEmbeddings(lengths)
