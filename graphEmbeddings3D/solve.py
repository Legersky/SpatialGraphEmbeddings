#Copyright (C) 2018 Jan Legersky
#
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <https://www.gnu.org/licenses/>.



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
