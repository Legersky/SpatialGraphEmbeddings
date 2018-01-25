#Copyright (C) 2018 Jan Legerský
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
from phcpy.trackers import track
try:
    from phcpy import _number_of_cores_ as tasks_num
except:
    tasks_num = 2

syst2 = ast.literal_eval(sys.argv[1])
#syst = ast.literal_eval(sys.argv[2])
#sols = ast.literal_eval(sys.argv[3])
pref = sys.argv[3]

filePrev = sys.argv[2]

with open(filePrev, 'r') as fPrev:
    syst, sols = [ast.literal_eval(line) for line in fPrev]

file = open(os.path.dirname(os.path.realpath(__file__))+'/../tmp/'+pref+'track.txt','w') 
sol = track(syst2,syst,sols, tasks=tasks_num)

file.write(str(sol))
