import sys
import ast
from phcpy.trackers import track
from phcpy.solver import solve
from phcpy.solutions import strsol2dict, is_real

try:
    from phcpy import _number_of_cores_ as tasks_num
except:
    tasks_num = 2

fileNamePref = sys.argv[1]
filePrev = sys.argv[2]
intervals = ast.literal_eval(sys.argv[3])
#prevSystem = ast.literal_eval(sys.argv[2])
#prevSolutions = ast.literal_eval(sys.argv[3])
#numAll = ast.literal_eval(sys.argv[4])
with open(filePrev, 'r') as fPrev:
    prevSystem, prevSolutions, numAll = [ast.literal_eval(line) for line in fPrev]
    
usePrev = True
tolerance = 1.0e-8
 
def findEmbeddings(syst):
    global prevSystem
    global usePrev
    global prevSolutions
    i = 0
    while True:
        if prevSystem and usePrev:
            sols = track(syst, prevSystem, prevSolutions, tasks=tasks_num)
        else:
            sols = solve(syst, verbose=0, tasks=tasks_num)

        result_real = []
        for sol in sols:
            soldic = strsol2dict(sol)
            if is_real(sol, tolerance):
                result_real.append(soldic)
        
        num_real = len(result_real)
        
        if num_real%4==0:    # and len(sols)==numAll:
            prevSystem = syst
            prevSolutions = sols
            return num_real
        else:
            usePrev = False
            i += 1
#            print 'PHC failed, trying again: '+str(i)
        if i>=3:
            print 'PHC failed 3 times'
            return -1

def findEmbeddings_dist(syst,  interval):
    global prevSystem
    global usePrev
    global prevSolutions
    i = 0
    y1_left,  y1_right, y4_left,  y4_right = interval
#    print fileNamePref
    while True:
        if prevSystem and usePrev:
            sols = track(syst, prevSystem, prevSolutions, tasks=tasks_num)
        else:
            sols = solve(syst, verbose=0, tasks=tasks_num)

        result_real = []
        for sol in sols:
            soldic = strsol2dict(sol)
            if is_real(sol, tolerance) and soldic['y1'].real>=y1_left and soldic['y1'].real<=y1_right and soldic['y4'].real>=y4_left and soldic['y4'].real<=y4_right:
                result_real.append(soldic)

        num_real = len(result_real)
        
        if num_real%2==0:    # and len(sols)==numAll:
            prevSystem = syst
            prevSolutions = sols
            return num_real
        else:
            usePrev = False
            i += 1
#            print 'PHC failed, trying again: '+str(i)
        if i>=3:
            print 'PHC failed 3 times'
            return -1


res = []
with open(fileNamePref+'eqs.txt', 'r') as fileEqs:
    if intervals:
        for eq, interval in zip(fileEqs, intervals):
            res.append(findEmbeddings_dist(ast.literal_eval(eq), interval))
    else:
        for eq in fileEqs:
            res.append(findEmbeddings(ast.literal_eval(eq)))


with open(fileNamePref+'numReal.txt','w') as file:
    file.write(str(res))










