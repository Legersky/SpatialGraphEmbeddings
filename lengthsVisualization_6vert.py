import ast
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


edges = [(1, 2), (2, 7), (4, 7), (2, 6), (6, 7), (4, 5), (5, 7), (1, 4), (1, 5), 
         (1, 3), (1, 6), (3, 7), (5, 6), (2, 3), (3, 4)]
edge2ind = {
    '27':0,'47':1,'26':2,'67':3,'45':4,'57':5,'14':6,'15':7,
    '13':8,'16':9,'37':10,'56':11,'23':12,'34':13
}
list_of_lengths = []
with open('../CouplerCurve/res/important/generated_48_from_first_in_7vert_73385d51.txt','r') as file:
    i = 0
    for L in file:
        lengths = ast.literal_eval(L)
        len_list = []
        for e in edges[1:]:
            len_list.append(lengths[e]/float(lengths[edges[0]]))
        list_of_lengths.append(len_list)
        i +=1
#         if i>=10:
#             break


print 'Number of lengths:'
print len(list_of_lengths)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

e1 = '23'
e2 = '14'
e3 = '45'
ax.scatter([L[edge2ind[e1]] for L in list_of_lengths],  
            [L[edge2ind[e2]] for L in list_of_lengths], 
            [L[edge2ind[e3]] for L in list_of_lengths])

ax2 = fig.add_subplot(211, projection='3d')

e1 = '27'
e2 = '16'
e3 = '14'
ax2.scatter([L[edge2ind[e1]] for L in list_of_lengths],  
            [L[edge2ind[e2]] for L in list_of_lengths], 
            [L[edge2ind[e3]] for L in list_of_lengths])
            
plt.show()


