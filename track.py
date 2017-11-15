import sys
import ast
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

file = open('tmp/'+pref+'track.txt','w') 
sol = track(syst2,syst,sols, tasks=tasks_num)

file.write(str(sol))
