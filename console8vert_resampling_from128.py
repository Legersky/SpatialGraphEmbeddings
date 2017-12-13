from algRealEmbeddings import *


L = [
#     {(2, 7): 10.667837135373238, (4, 7): 44.550486792918335, (2, 6): 1.6141574401570116, (4, 5): 72.35604118410009, (2, 8): 1.6175361079301909,
#       (1, 4): 44.81372216057459, (5, 6): 57.11467444895333, (1, 3): 13.722832335192432, (1, 6): 14.58705569793232, (3, 7): 8.570089214958891, 
#       (5, 8): 57.10811293639709, (1, 2): 14.655598375989243, (6, 8): 0.3474654025573003, (3, 4): 44.95202906796581, (5, 7): 58.21903364628109, 
#       (1, 5): 58.926778412561376, (2, 3): 10.149107373916825, (7, 8): 10.628404786753377}, 
    {(2, 7): 10.756167139788014, (4, 7): 44.04766504173981, (2, 6): 1.6141574401570116, (4, 5): 78.08873019522255, (2, 8): 1.6175361079301909, 
    (1, 4): 45.64217524019533, (5, 6): 57.11467444895333, (1, 3): 14.525549875177061, (1, 6): 14.645942367542169, (3, 7): 7.778046681061985, 
    (5, 8): 57.10811293639709, (1, 2): 14.734786893335137, (6, 8): 0.47609964944060706, (3, 4): 45.96689412631051, (5, 7): 57.125075814116045, 
    (1, 5): 57.77282917982447, (2, 3): 9.191326751698375, (7, 8): 10.540074782338603}, 
    {(2, 7): 10.579507130958465, (4, 7): 45.05330854409685, (2, 6): 1.6141574401570116, (4, 5): 66.01696292520337, (2, 8): 1.6175361079301909, 
    (1, 4): 43.99334849111959, (5, 6): 57.11467444895333, (1, 3): 12.935609795711123, (1, 6): 14.52816902832247, (3, 7): 9.367986760527923, 
    (5, 8): 57.10811293639709, (1, 2): 14.558126884505588, (6, 8): 0.2405529710012077, (3, 4): 43.99334849111959, (5, 7): 59.32938264884837, 
    (1, 5): 60.06940929460853, (2, 3): 11.096653871614114, (7, 8): 10.716734791168152}
     ]



for lengths in L:
    alg = AlgRealEmbeddings('Max8vertices', name='8vert_afterResamplingFrom128')
    alg.findMoreEmbeddings(lengths)









