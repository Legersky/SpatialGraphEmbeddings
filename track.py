import sys
import ast
from phcpy.trackers import track


syst2 = ast.literal_eval(sys.argv[1])
syst = ast.literal_eval(sys.argv[2])
sols = ast.literal_eval(sys.argv[3])
pref = sys.argv[4]


file = open('tmp/'+pref+'track.txt','w') 
sol = track(syst2,syst,sols, tasks=2)

file.write(str(sol))
