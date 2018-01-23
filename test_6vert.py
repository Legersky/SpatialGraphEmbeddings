
# This script uses the sampling method in order to find edge lengths of the only H2 6-vertex graph (cyclohexane) with many real embeddings
# Random embedding is used for obtaining starting edge lengths
# It is likely that edge lengths with 12 or 16 are found

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
#lengths_final = alg.findMoreEmbeddings(lengths, allowed_repetition=2)
lengths_final = {(1, 2): 1.2895877366356432, (2, 6): 0.9400153833621129, (4, 6): 0.8014773302862134, (5, 6): 2.7367590363070198, (1, 5): 2.911314881424429, (1, 3): 0.68614667990462808, (2, 3): 1.5023140682372356, (4, 5): 2.750838716246279, (1, 6): 1.096677140268674, (3, 4): 1.626487374448872, (2, 4): 0.9629467154119479, (3, 5): 2.905105239383198}

print '\n\nChecking and summary '
G = GraphEmbedding(lengths_final, 'Max6vertices')
print 'There are ', len(G.findEmbeddings()['real']), 'real embeddings of the cyclohexane with the following edge lengths:'
print lengths_final
