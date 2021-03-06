'''
This script verifies that the 8-vertex graph G128, resp. G160, has edge lengths such that 128, resp. 132, embeddings are real.
The edge lengths were found by the sampling method.
'''

if __name__ == '__main__':
    print __doc__

    from graphEmbeddings3D.graphEmbedding import GraphEmbedding

#    print 'Graph G128'
#    lengths_G128 = {(2, 7): 9.2728, (1, 2): 8.7093, (1, 3): 10.3433, (6, 7): 1.5449, (6, 8): 10.5532, (4, 8): 10.5237, 
#                    (5, 6): 0.7536, (2, 8): 13.5773, (3, 4): 10.1636, (1, 4): 1.9373, (3, 8): 14.6173, (1, 5): 1.9379, 
#                    (1, 6): 2.0691, (1, 7): 2.1185, (4, 5): 0.0634, (2, 3): 13.5267, (7, 8): 10.5509, (5, 8): 10.5237}
#    G = GraphEmbedding(lengths_G128, 'Ring8vertices')
#    print 'There are ', len(G.findEmbeddings(tolerance=1.0e-50)['real']), 'real embeddings of G128 with the following edge lengths:'
#    print lengths_G128,  '\n'

    print 'Graph G160'
    lengths_G160 = {(1, 2): 1.999, (4, 7): 11.993, (2, 6): 0.879, (4, 5): 7.278, (2, 8): 0.847, (1, 4):6.611,
                    (1, 3): 1.568, (1, 6): 1.994, (3, 7): 10.447, (5, 8):4.279, (2, 7): 10.536, (6, 8):0.398,
                    (1, 5): 4.402, (3, 4): 6.494, (5, 7): 11.239, (5, 6): 4.321, (2, 3): 1.426, (7, 8): 10.474}
    G = GraphEmbedding(lengths_G160, 'Max8vertices')
    print 'There are ', len(G.findEmbeddings(tolerance=1.0e-50)['real']), 'real embeddings of G160 with the following edge lengths:'
    print lengths_G160,  '\n'
