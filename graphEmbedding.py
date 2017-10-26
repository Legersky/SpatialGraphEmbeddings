import numpy as np
#import math
from phcpy.solver import solve
from phcpy.trackers import track
from phcpy.solutions import strsol2dict, is_real
#from sympy import symbols
#import time
#import copy
#from sklearn.cluster import DBSCAN


class GraphEmbedding(object):
    def __init__(self, lengths, equationsConstructor, window=None):
        self._window=window
        self.setLengths(lengths)
        
        self._equationsConstructor = equationsConstructor

        self._prevSystem = None
        self._prevSolutions = None

    def getEdgeLength(self, edge):
        return float(self._lengths['L'+edge])

    def setLengths(self, lengths):
        self._lengths = {}
        try:
            for el in lengths.keys():
                self._lengths[el] = float(lengths[el])
                self.printLog(str(el)+': '+str(self._lengths[el]), verbose=2)
        except KeyError as e:
            self.printLog('Problem with setting lengths: '+str(e))

    def printLog(self, s, verbose=0):
        if self._window:
            self._window.printLog(s, verbose=verbose)
        else:
            print s

    def getEquations(self):
        return self._equationsConstructor()
    
    def dist(self, u, v):
        return float(np.sqrt( (u[0]-v[0])**2 + (u[1]-v[1])**2 + (u[2]-v[2])**2))

    def findEmbeddings(self, tolerance=1.0e-8,  errorMsg=True, usePrev=True):
        syst = self.getEquations()
        i = 0
        while True:
            if self._prevSystem and usePrev:
                sols = track(syst, self._prevSystem, self._prevSolutions, tasks=2)
            else:
                sols = solve(syst, verbose=0, tasks=2)
        
            result_real = []
            result_complex = []
            for sol in sols:
                soldic = strsol2dict(sol)
                if is_real(sol, tolerance):
                    result_real.append(soldic)
                else:
                    result_complex.append(soldic)
            
            num_real = len(result_real)
            num_all = len(sols)
            
            if num_real % 4 == 0 and num_all==48:
                self._prevSystem = syst
                self._prevSolutions = sols
                return {'real':result_real, 'complex':result_complex}
            else:
                usePrev = False
                i += 1
                self.printLog('PHC failed, trying again: '+str(i), verbose=1)
            if i>=10:
                self._prevSystem = []
                if errorMsg:
                    self._window.showError('PHC failed 10 times')
                print 'PHC failed 10 times'
                return {'real':[], 'complex':[]}
