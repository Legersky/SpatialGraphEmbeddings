import sys
import ast
from phcpy.solver import solve


syst = ast.literal_eval(sys.argv[1])

file = open('tmp/solve.txt','w') 
sol = solve(syst, verbose=0, tasks=2)

file.write(str(sol))
