'''
 This script uses the sampling method in order to find edge lengths of the only H2 6-vertex graph (G16, cyclohexane) with many real embeddings.
 Random embedding is used for obtaining starting edge lengths.
 It is likely that edge lengths with 12 or 16 are found.
'''
if __name__ == '__main__':
    print __doc__


    from graphEmbeddings3D.algRealEmbeddings import AlgRealEmbeddings
    from graphEmbeddings3D.graphEmbedding import getEdgeLengthsByEmbedding,  GraphEmbedding
    from random import uniform
    import time

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
    # or the following edge lengths yields ones with 16 real embeddings:
    lengths = {(1, 3): 0.5308353556228526, (5, 6): 1.8056034169357076, (2, 6): 2.4456570316770256, (2, 3): 0.8512597859002053,
               (3, 5): 1.2853151143150572, (1, 2): 1.1757548283222, (4, 6): 0.9879010951471968, (1, 5): 0.7810003670625068,
               (4, 5): 1.6175092476119708, (1, 6): 1.4488010671291325, (3, 4): 1.1784501009786457, (2, 4): 1.9763357382306461}
    
    start = time.time()
    # we run the sampling method based on coupler curves
    alg = AlgRealEmbeddings('Max6vertices', num_phi=5, num_theta=5,  name='random6vert')
    lengths_final = alg.findMoreEmbeddings(lengths, allowed_repetition=2)
    
    print '\n\nSampling time: ',  time.time() - start

    print '\n\nChecking and summary '
    G = GraphEmbedding(lengths_final, 'Max6vertices')
    print 'There are ', len(G.findEmbeddings()['real']), 'real embeddings of the cyclohexane with the following edge lengths:'
    print lengths_final
