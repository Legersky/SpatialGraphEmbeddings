'''
This script verifies that all 7-vertex graphs have edge lengths such that all embeddings are real
The edge lengths were found by the sampling method except for G16a where they were obtained just by random guess.
'''
if __name__ == '__main__':
    print __doc__


    from graphEmbeddings3D.graphEmbedding import GraphEmbedding

    print 'Graph G16a'
    lengths_G16a = {(4, 7): 7.19, (1, 3): 5.75, (5, 6): 7.9, (1, 6): 8.48, (3, 7): 5.91, (2, 5): 7.15, (3, 5): 5.09, 
                    (1, 2): 4.36, (4, 6): 8.78, (5, 7): 10.22, (3, 6): 7.06, (1, 7): 3.77, (2, 3): 3.81, (3, 4): 3.23, (2, 4): 6.05}
    G = GraphEmbedding(lengths_G16a, '7vert16a')
    print 'There are ', len(G.findEmbeddings(tolerance=1.0e-50)['real']), 'real embeddings of G16a with the following edge lengths:'
    print lengths_G16a,  '\n'


    print 'Graph G16b'
    lengths_G16b = {(1, 2): 4.62, (4, 7): 4.46, (2, 6): 7.47, (4, 5): 7.72, (1, 4): 6.51, (1, 3): 3.53, (2, 3): 7.69, (3, 7): 5.76,
                    (2, 5): 9.48, (3, 5): 6.1, (2, 7): 5.9, (6, 7): 3.09, (4, 6): 7.07, (1, 5): 5.69, (3, 6): 6.43}
    G = GraphEmbedding(lengths_G16b, '7vert16b')
    print 'There are ', len(G.findEmbeddings(tolerance=1.0e-50)['real']), 'real embeddings of G16b with the following edge lengths:'
    print lengths_G16b,  '\n'


    print 'Graph G24'
    lengths_G24 = {(1, 2): 11.05, (4, 7): 5.65, (2, 6): 5.7, (5, 6): 4.7, (1, 4): 8.33, (1, 3): 4.77, (2, 3): 10.31, (3, 7): 7.1, 
                   (2, 5): 9.32, (2, 7): 6.0, (4, 6): 6.49, (5, 7): 5.77, (1, 5): 9.4, (3, 6): 8.57, (3, 4): 7.64}
    G = GraphEmbedding(lengths_G24, '7vert24')
    print 'There are ', len(G.findEmbeddings(tolerance=1.0e-50)['real']), 'real embeddings of G24 with the following edge lengths:'
    print lengths_G24,  '\n'


    print 'Graph G32a'
    lengths_G32a = {(4, 7): 8.74, (1, 3): 6.27, (5, 6): 9.23, (1, 4): 8.06, (2, 3): 8.83, (3, 7): 5.62, (2, 5): 9.74, (3, 5): 5.6, 
                    (1, 2): 10.95, (6, 7): 9.28, (5, 7): 7.88, (3, 6): 8.26, (1, 6): 11.56, (3, 4): 6.11, (2, 4): 8.95}
    G = GraphEmbedding(lengths_G32a, '7vert32a')
    print 'There are ', len(G.findEmbeddings(tolerance=1.0e-50)['real']), 'real embeddings of G32a with the following edge lengths:'
    print lengths_G32a,  '\n'


    print 'Graph G32b'
    lengths_G32b = {(1, 2): 11.06, (4, 7): 85.49, (2, 6): 7.11, (4, 5): 78.53, (1, 4): 87.33, (5, 6): 22.08, (1, 3): 10.81, 
                    (2, 3): 4.47, (3, 7): 7.1, (2, 5): 20.7, (2, 7): 7.68, (6, 7): 9.29, (1, 5): 21.49, (3, 6): 7.53, (3, 4): 84.17}
    G = GraphEmbedding(lengths_G32b, '7vert32b')
    print 'There are ', len(G.findEmbeddings(tolerance=1.0e-50)['real']), 'real embeddings of G32b with the following edge lengths:'
    print lengths_G32b,  '\n'


    print 'Graph G48'
    lengths_G48 = {(2, 7): 10.5361, (4, 7): 11.8471, (2, 6): 1.002, (5, 6): 4.4449, (1, 4): 5.7963, (1, 3): 1.9342, (1, 6): 2.0001, (3, 7): 10.5245, 
                   (1, 2): 1.9999, (6, 7): 10.5365, (5, 7): 11.2396, (1, 5): 4.4024, (4, 5): 7.0744, (2, 3): 0.55, (3, 4): 5.4247}
    G = GraphEmbedding(lengths_G48, 'Max7vertices')
    print 'There are ', len(G.findEmbeddings(tolerance=1.0e-50)['real']), 'real embeddings of G48 with the following edge lengths:'
    print lengths_G48,  '\n'

