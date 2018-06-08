'''
This script tries to find edge lengths of G160 with many real embeddings stochastically.
'''



if __name__ == '__main__':
    print __doc__
    from random import random
    
    from graphEmbeddings3D.graphEmbedding import GraphEmbedding

    lengths_G160 = {(1, 2): 1.999, (4, 7): 11.993, (2, 6): 0.879, (4, 5): 7.278, (2, 8): 0.847, (1, 4):6.611,
                    (1, 3): 1.568, (1, 6): 1.994, (3, 7): 10.447, (5, 8):4.279, (2, 7): 10.536, (6, 8):0.398,
                    (1, 5): 4.402, (3, 4): 6.494, (5, 7): 11.239, (5, 6): 4.321, (2, 3): 1.426, (7, 8): 10.474}
    G = GraphEmbedding(lengths_G160, 'Max8vertices')
    G.findEmbeddings()
    edges = [(1, 2), (2, 7), (4, 7), (2, 6), (6, 8), (4, 5), (2, 8), (5, 7), (7, 8), (1, 4), (1, 5), (1, 3), (1, 6), (5, 6), (3, 7), (3, 4), (2, 3), (5, 8)]

    for _ in range(0,10):
        for e in edges:
            lengths_G160[e] = 1+0.01*random()
        G.setLengths(lengths_G160)
        real_sol = G.findEmbeddings(tolerance=1.0e-50)['real']
        print len(real_sol)
        if len(real_sol)>=132:
            with open('stoch_res.txt', 'a') as f:
                f.write(str(len(real_sol))+'\n')
                f.write(str(lengths_G160)+'\n\n')


