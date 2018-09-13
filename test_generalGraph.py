'''
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
    v7 = [uniform(a, b), uniform(a, b), uniform(a, b)]
    v8 = [uniform(a, b), uniform(a, b), uniform(a, b)]

    # edge lengths are given as the distances in the embedding
    graphs = [[(1, 2), (1, 3), (1, 5), (1, 8), (2, 3), (2, 4), (2, 8), (3, 4), (3, 5), 
                (3, 6), (3, 7), (4, 7), (4, 8), (5, 6), (5, 8), (6, 7), (6, 8), (7, 8)],
                [(1, 2), (1, 3), (1, 4), (1, 8), (2, 3), (2, 5), (2, 6), (2, 8), (3, 4),
                (3, 5), (3, 7), (4, 7), (4, 8), (5, 6), (5, 7), (6, 7), (6, 8), (7, 8)], 
                [(1, 2), (1, 3), (1, 4), (1, 5), (1, 6), (2, 3), (2, 4), (2, 5), (2, 7),
                (3, 6), (3, 7), (3, 8), (4, 7), (4, 8), (5, 6), (5, 8), (6, 8), (7, 8)], 
                [(1, 2), (1, 3), (1, 7), (1, 8), (2, 3), (2, 5), (2, 6), (2, 7), (3, 4), 
                (3, 6), (3, 8), (4, 5), (4, 7), (4, 8), (5, 7), (5, 8), (6, 7), (6, 8)],
                [(1, 2), (1, 3), (1, 4), (1, 8), (2, 3), (2, 5), (2, 6), (2, 7), (3, 4),
                (3, 5), (3, 7), (3, 8), (4, 5), (4, 8), (5, 6), (6, 7), (6, 8), (7, 8)],
                [(1, 2), (1, 3), (1, 4), (1, 6), (1, 7), (2, 3), (2, 5), (2, 7), (2, 8),
                (3, 4), (3, 6), (3, 8), (4, 5), (4, 7), (5, 6), (5, 8), (6, 8), (7, 8)], 
                [(1, 2), (1, 3), (1, 4), (1, 6), (1, 7), (2, 3), (2, 5), (2, 6), (2, 8), 
                (3, 4), (3, 7), (3, 8), (4, 5), (4, 7), (5, 6), (5, 8), (6, 8), (7, 8)], 
                [(1, 2), (1, 3), (1, 4), (1, 6), (1, 7), (2, 3), (2, 4), (2, 5), (2, 8), 
                (3, 6), (3, 7), (3, 8), (4, 5), (4, 7), (5, 6), (5, 8), (6, 8), (7, 8)], 
                [(1, 2), (1, 3), (1, 4), (1, 8), (2, 3), (2, 6), (2, 7), (2, 8), (3, 4), 
                (3, 5), (3, 6), (3, 7), (4, 5), (4, 7), (5, 6), (5, 8), (6, 8), (7, 8)], 
                [(1, 2), (1, 3), (1, 4), (1, 8), (2, 3), (2, 4), (2, 6), (2, 7), (3, 5), 
                (3, 6), (3, 7), (3, 8), (4, 5), (4, 7), (5, 6), (5, 8), (6, 8), (7, 8)], 
                [(1, 2), (1, 3), (1, 7), (1, 8), (2, 3), (2, 5), (2, 6), (2, 8), (3, 4), 
                (3, 5), (3, 7), (4, 5), (4, 6), (4, 8), (5, 7), (6, 7), (6, 8), (7, 8)], 
                [(1, 2), (1, 3), (1, 7), (1, 8), (2, 3), (2, 4), (2, 6), (2, 7), (3, 4), 
                (3, 5), (3, 8), (4, 5), (4, 6), (5, 7), (5, 8), (6, 7), (6, 8), (7, 8)], 
                [(1, 2), (1, 3), (1, 5), (1, 8), (2, 3), (2, 4), (2, 7), (2, 8), (3, 4), 
                (3, 5), (3, 6), (3, 7), (4, 5), (4, 6), (5, 8), (6, 7), (6, 8), (7, 8)], 
                [(1, 2), (1, 3), (1, 6), (1, 8), (2, 3), (2, 6), (2, 7), (2, 8), (3, 4), 
                (3, 5), (3, 6), (3, 7), (4, 5), (4, 6), (4, 8), (5, 7), (5, 8), (7, 8)], 
                [(1, 2), (1, 3), (1, 5), (1, 6), (1, 7), (2, 3), (2, 4), (2, 5), (2, 7), 
                (3, 6), (3, 7), (3, 8), (4, 5), (4, 6), (4, 8), (5, 8), (6, 8), (7, 8)], 
                [(1, 2), (1, 3), (1, 4), (1, 6), (1, 7), (2, 3), (2, 5), (2, 6), (2, 8), 
                (3, 5), (3, 7), (3, 8), (4, 5), (4, 6), (4, 8), (5, 7), (6, 8), (7, 8)], 
                [(1, 2), (1, 3), (1, 4), (1, 6), (1, 7), (2, 3), (2, 5), (2, 7), (2, 8), 
                (3, 6), (3, 7), (3, 8), (4, 5), (4, 6), (4, 8), (5, 7), (5, 8), (6, 8)], 
                [(1, 2), (1, 3), (1, 4), (1, 5), (1, 6), (2, 3), (2, 6), (2, 7), (2, 8), 
                (3, 5), (3, 7), (3, 8), (4, 5), (4, 6), (4, 7), (5, 8), (6, 8), (7, 8)]]
    
    numbersOfComplexSolutions = [128, 160, 128, 112, 96, 96, 96, 96, 80, 80, 80, 128, 96, 96, 96, 80, 96, 96]
    names = ['G128a', 'G160', 'G128b', 'G112', 'G96a', 'G96b', 'G96c', 'G96d', 'G80a', 'G80b', 'G80c', 'G128c', 'G96', 'G96e', 'G96f', 'G80d', 'G96g', 'G96h']
    for N in range(2, len(graphs)):
        lengths = getEdgeLengthsByEmbedding('edges',  [v1, v2, v3, v4, v5, v6, v7, v8], edges=graphs[N])
        print graphs[N]
        start = time.time()
        # we run the sampling method based on coupler curves
        alg = AlgRealEmbeddings('edges', name=names[N],  edges=graphs[N],  num_sols=numbersOfComplexSolutions[N],  allowedNumberOfMissing=10)
        lengths_final = alg.findMoreEmbeddings(lengths)
        
        print '\n\nSampling time: ',  time.time() - start

        print '\n\nChecking and summary '
        G = GraphEmbedding(lengths_final, 'edges', num_sols=numbersOfComplexSolutions[N])
        print 'There are ', len(G.findEmbeddings()['real']), 'real embeddings with the following edge lengths:'
        print lengths_final
