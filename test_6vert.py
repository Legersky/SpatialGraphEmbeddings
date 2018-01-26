'''
 This script uses the sampling method in order to find edge lengths of the only H2 6-vertex graph (cyclohexane) with many real embeddings.
 Random embedding is used for obtaining starting edge lengths.
 It is likely that edge lengths with 12 or 16 are found.
'''

print __doc__


from graphEmbeddings3D.algRealEmbeddings import AlgRealEmbeddings
from graphEmbeddings3D.graphEmbedding import getEdgeLengthsByEmbedding,  GraphEmbedding
from random import uniform

a = -1.0
b = 1.0
# a random embedding is taken:
v1 = [uniform(a, b), uniform(a, b), uniform(a, b)]
v2 = [uniform(a, b), uniform(a, b), uniform(a, b)]
v3 = [uniform(a, b), uniform(a, b), uniform(a, b)]
v4 = [uniform(a, b), uniform(a, b), uniform(a, b)]
v5 = [uniform(a, b), uniform(a, b), uniform(a, b)]
v6 = [uniform(a, b), uniform(a, b), uniform(a, b)]

# edge lengths are given as the distances in the embedding
lengths = getEdgeLengthsByEmbedding('Max6vertices', [v1, v2, v3, v4, v5, v6])

# we run the sampling method based on coupler curves
alg = AlgRealEmbeddings('Max6vertices', num_phi=5, num_theta=5,  name='random6vert')
lengths_final = alg.findMoreEmbeddings(lengths, allowed_repetition=2)

print '\n\nChecking and summary '
G = GraphEmbedding(lengths_final, 'Max6vertices')
print 'There are ', len(G.findEmbeddings()['real']), 'real embeddings of the cyclohexane with the following edge lengths:'
print lengths_final
