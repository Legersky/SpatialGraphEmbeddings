import sys
import ast
import os
from phcpy.solver import solve
try:
    from phcpy import _number_of_cores_ as tasks_num
except:
    tasks_num = 2

syst = ast.literal_eval(sys.argv[1])
pref = sys.argv[2]

with open(os.path.dirname(os.path.realpath(__file__))+'/../tmp/'+pref+'solve.txt','w') as file:
    sol = solve(syst, verbose=0, tasks=tasks_num)

    file.write(str(sol))