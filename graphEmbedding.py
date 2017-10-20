import numpy as np
import math
from phcpy.solver import solve
from phcpy.solutions import strsol2dict, is_real
from sympy import symbols
import copy


class GraphEmbedding(object):
    def __init__(self,lengths={}, r26=1.0,   window=None):
        self._window=window
        if not lengths:
            lengths = {
                       'L12': 10.4093719791350,
                       'L13': 10.3120185827024,
                       'L14': 10.3156349780321,
                       'L15': 10.3266511996871,
                       'L16': 10.3123724234533,
                       'L27': 2.44436085715673,
                       'L37': 1.99846497342335,
                       'L47': 1.98932149236869,
                       'L57': 2.05635810311561,
                       'L67': 1.99371512508683,
                       'L23': 1.01300160414483,
                       'L34': 1.00554077490672,
                       'L45': 0.996142560078625,
                       'L56': 1.07000000000000,
                       } 
        self.setLengths(lengths)
        self.setR26(r26)
        self._branches = {}
        self._branches_mirror = {}
        self._samples = 0

    def setRequiresRecomputing(self, propagate=True):
        self._branches_updated = False
        try:
            if self._window and propagate:
                self._window.setRequiresRecomputing()
        except AttributeError:
            pass

    def setComputed(self, propagate=True):
        self._branches_updated = True
        if self._window and propagate:
            self._window.setComputed()
    
    def isComputed(self):
        return self._branches_updated

    def getEdgeLength(self, edge):
        return float(self._lengths['L'+edge])

    def getSamples(self):
        return self._samples

    def setR26(self, r26):
        self._R26 = float(r26)
    
    def getR26(self):
        return self._R26

    def getyV2(self):
        return self._v2[1]
    

    def setLengths(self, lengths):
        self._lengths = {}
        try:
            for el in ['L67', 'L47', 'L56', 'L57', 'L14', 'L15', 'L16', 'L12', 'L13', 'L45', 'L34', 'L37', 'L23', 'L27']:
                self._lengths[el] = float(lengths[el])
                self.printLog(str(el)+': '+str(self._lengths[el]))
        except KeyError as e:
            self.printLog('Problem with setting lengths: '+str(e))
        
        l12 = self._lengths['L12']
        l13 = self._lengths['L13']
        l27 = self._lengths['L27']
        l37 = self._lengths['L37']
        l23 = self._lengths['L23']
        
        theta1 = math.acos((-l13**2+l12**2+l23**2)/(2*l12*l23))
        theta7 = math.acos((-l37**2+l27**2+l23**2)/(2*l27*l23))

        self._max_x7 = math.sin(theta7)*l27
        x1 = math.sin(theta1)*l12

        y7 = math.cos(theta7)*l27
        y1 = math.cos(theta1)*l12

        y2 = -y7
        y3 = -y7+l23
        y1 = y1-y7
        
        self._v1 = [x1, y1, 0]
        self._v2 = [0, y2, 0]
        self._v3 = [0, y3, 0]
        self.setRequiresRecomputing()
    
    def setyV2(self, y2):
        self._v2 = [0, y2, 0]
        
        self._lengths['L12'] = dist(self._v1, self._v2)
        self._lengths['L27'] = dist([self._max_x7, 0, 0], self._v2)
        self._lengths['L23'] = dist(self._v3, self._v2)
    
    def printLog(self, s):
        if self._window:
            self._window.printLog(s)
        else:
            print s
            
    def getCenterOfGravity(self):
        try:
            return self._centerOfGravity
        except:
            return [0, 0, 0]

    def getBranch(self, color):
        if self.isComputed():
            return self._branches[color]
        else:
            self.printLog('Coupler curve must be computed first.')
            return []

    def getMirrorBranch(self, color):
        if self.isComputed():
            return self._branches_mirror[color]
        else:
            self.printLog('Coupler curve must be computed first.')
            return []

    def computeCouplerCurve(self, N):
        self._samples = N
        
        l14 = self._lengths['L14']
        l15 = self._lengths['L15']
        l16 = self._lengths['L16']
        l47 = self._lengths['L47']
        l57 = self._lengths['L57']
        l67 = self._lengths['L67']
        l34 = self._lengths['L34']
        l45 = self._lengths['L45']
        l56 = self._lengths['L56']


        v3 = self._v3
        v1 = self._v1
        
        self.minbound = 0
        self.maxbound = 0
        
#        def getTraceV6(solv4, solv5, solv6):
#            trace_v6 = [[]]
#            for i in range(0,N):
#                v7_x = float(self._max_x7*math.sin(2*math.pi*i/N))
#                v7_z = float(self._max_x7*math.cos(2*math.pi*i/N))
#                v7 = [v7_x, 0, v7_z]
#        #         trace_v7.append(v7)
#                intersection_v4 = getIntersectionOfThreeSpheres(v3,v1,v7,l34,l14,l47)
#                if intersection_v4 != None:
#                    intersection_v5 = getIntersectionOfThreeSpheres(intersection_v4[solv4],v1,v7,l45,l15,l57)
#                    if intersection_v5 != None:
#                        intersection_v6 = getIntersectionOfThreeSpheres(intersection_v5[solv5],v1,v7,l56,l16,l67)
#                        if intersection_v6 != None:
#                            v6 = intersection_v6[solv6]
#                            trace_v6[-1].append(v6)
#                            self.minbound = min(self.minbound, v6[0], v6[1], v6[2])
#                            self.maxbound = max(self.maxbound, v6[0], v6[1], v6[2])
#                        else:
#                            trace_v6.append([])
#                    else:
#                        trace_v6.append([])
#                else:
#                    trace_v6.append([])
#            return trace_v6
        
        self._trace_v7 = []
        self._trace_v4 = [[], []]
        self._trace_v5 = [[[], []], [[], []]]
        self._trace_v6 = [[[[], []], [[], []]], [[[], []], [[], []]]]
        
        for i in range(0,N):
            v7_x = float(self._max_x7*math.sin(2*math.pi*i/N))
            v7_z = float(self._max_x7*math.cos(2*math.pi*i/N))
            v7 = [v7_x, 0, v7_z]
            self._trace_v7.append(v7)
            
            intersection_v4 = getIntersectionOfThreeSpheres(v3,v1,v7,l34,l14,l47)
            v4 = intersection_v4[0]
            self._trace_v4[0].append(v4)
            if True:
                intersection_v5 = getIntersectionOfThreeSpheres(v4,v1,v7,l45,l15,l57)
                v5 = intersection_v5[0]
                self._trace_v5[0][0].append(v5)
                if True:
                    intersection_v6 = getIntersectionOfThreeSpheres(v5,v1,v7,l56,l16,l67)
                    v6 = intersection_v6[0]
                    self._trace_v6[0][0][0].append(v6)
                    if v6 != None:
                        self.minbound = min(self.minbound, v6[0], v6[1], v6[2])
                        self.maxbound = max(self.maxbound, v6[0], v6[1], v6[2])
                    
                    v6 = intersection_v6[1]
                    self._trace_v6[0][0][1].append(v6)
                    if v6 != None:
                        self.minbound = min(self.minbound, v6[0], v6[1], v6[2])
                        self.maxbound = max(self.maxbound, v6[0], v6[1], v6[2])
                
                v5 = intersection_v5[1]
                self._trace_v5[0][1].append(v5)
                if True:
                    intersection_v6 = getIntersectionOfThreeSpheres(v5,v1,v7,l56,l16,l67)
                    v6 = intersection_v6[0]
                    self._trace_v6[0][1][0].append(v6)
                    if v6 != None:
                        self.minbound = min(self.minbound, v6[0], v6[1], v6[2])
                        self.maxbound = max(self.maxbound, v6[0], v6[1], v6[2])
                    
                    v6 = intersection_v6[1]
                    self._trace_v6[0][1][1].append(v6)
                    if v6 != None:
                        self.minbound = min(self.minbound, v6[0], v6[1], v6[2])
                        self.maxbound = max(self.maxbound, v6[0], v6[1], v6[2])
            
            v4 = intersection_v4[1]
            self._trace_v4[1].append(v4)
            if True:
                intersection_v5 = getIntersectionOfThreeSpheres(v4,v1,v7,l45,l15,l57)
                v5 = intersection_v5[0]
                self._trace_v5[1][0].append(v5)
                if True:
                    intersection_v6 = getIntersectionOfThreeSpheres(v5,v1,v7,l56,l16,l67)
                    v6 = intersection_v6[0]
                    self._trace_v6[1][0][0].append(v6)
                    if v6 != None:
                        self.minbound = min(self.minbound, v6[0], v6[1], v6[2])
                        self.maxbound = max(self.maxbound, v6[0], v6[1], v6[2])
                    
                    v6 = intersection_v6[1]
                    self._trace_v6[1][0][1].append(v6)
                    if v6 != None:
                        self.minbound = min(self.minbound, v6[0], v6[1], v6[2])
                        self.maxbound = max(self.maxbound, v6[0], v6[1], v6[2])
                
                v5 = intersection_v5[1]
                self._trace_v5[1][1].append(v5)
                if True:
                    intersection_v6 = getIntersectionOfThreeSpheres(v5,v1,v7,l56,l16,l67)
                    v6 = intersection_v6[0]
                    self._trace_v6[1][1][0].append(v6)
                    if v6 != None:
                        self.minbound = min(self.minbound, v6[0], v6[1], v6[2])
                        self.maxbound = max(self.maxbound, v6[0], v6[1], v6[2])
                    
                    v6 = intersection_v6[1]
                    self._trace_v6[1][1][1].append(v6)
                    if v6 != None:
                        self.minbound = min(self.minbound, v6[0], v6[1], v6[2])
                        self.maxbound = max(self.maxbound, v6[0], v6[1], v6[2])
            

        def getTraceV6(a, b, c):
            res = [[]]
            for v in self._trace_v6[a][b][c]:
                if v != None:
                    res[-1].append(v)
                else:
                    if res[-1]:
                        res.append([])
            return res
            

        self._branches['orange'] = getTraceV6(0, 0, 0)
        self._branches['red'] = getTraceV6(0, 1, 0)
        self._branches['green'] = getTraceV6(0, 0, 1)
        self._branches['blue'] = getTraceV6(0, 1, 1)
        
        self._branches_mirror['orange'] = getTraceV6(1, 1, 1)
        self._branches_mirror['red'] = getTraceV6(1, 0, 1)
        self._branches_mirror['green'] = getTraceV6(1, 1, 0)
        self._branches_mirror['blue'] = getTraceV6(1, 0, 0)
        
        s_x, s_y, s_z = 0, 0, 0
        for color in ['orange', 'red', 'green', 'blue']:
            for part in self._branches_mirror[color]+self._branches[color]:
                for x, y, z in part:
                    s_x += x
                    s_y += y
                    s_z += z
        
        self._centerOfGravity = [s_x/float(8*N), s_y/float(8*N), s_z/float(8*N)]
        
        self.computeIntersections(nocheck=True)
        
        self.setComputed()

    def computeIntersections(self, nocheck=False):
        if self.isComputed() or nocheck:
            def getIntersectionWithSphere(A, B):
                x_A, y_A, z_A = A
                x_B, y_B, z_B = B
                x_C, y_C, z_C = self._v2
                r = self._R26
                res = []
                nom_a = x_A**2 - x_A*x_B - (x_A - x_B)*x_C + y_A**2 - y_A*y_B - (y_A - y_B)*y_C + z_A**2 - z_A*z_B - (z_A - z_B)*z_C
                D = (r**2*x_A**2 - 2*r**2*x_A*x_B + r**2*x_B**2 + (r**2 - x_B**2 + 2*x_B*x_C - x_C**2)*y_A**2 
                     - 2*(r**2 - x_A*x_B + (x_A + x_B)*x_C - x_C**2)*y_A*y_B + (r**2 - x_A**2 + 2*x_A*x_C - x_C**2)*y_B**2 
                     - (x_A**2 - 2*x_A*x_B + x_B**2)*y_C**2 + (r**2 - x_B**2 + 2*x_B*x_C - x_C**2 - y_B**2 + 2*y_B*y_C - y_C**2)*z_A**2 
                     - 2*(r**2 - x_A*x_B + (x_A + x_B)*x_C - x_C**2 - y_A*y_B + (y_A + y_B)*y_C - y_C**2)*z_A*z_B 
                     + (r**2 - x_A**2 + 2*x_A*x_C - x_C**2 - y_A**2 + 2*y_A*y_C - y_C**2)*z_B**2 
                     - (x_A**2 - 2*x_A*x_B + x_B**2 + y_A**2 - 2*y_A*y_B + y_B**2)*z_C**2 - 2*((x_A*x_B - x_B**2 - (x_A - x_B)*x_C)*y_A 
                     - (x_A**2 - x_A*x_B - (x_A - x_B)*x_C)*y_B)*y_C 
                     - 2*((x_A*x_B - x_B**2 - (x_A - x_B)*x_C + y_A*y_B - y_B**2 - (y_A - y_B)*y_C)*z_A - (x_A**2 - x_A*x_B - (x_A - x_B)*x_C + y_A**2 - y_A*y_B - (y_A - y_B)*y_C)*z_B)*z_C)
                denom = (x_A**2 - 2*x_A*x_B + x_B**2 + y_A**2 - 2*y_A*y_B + y_B**2 + z_A**2 - 2*z_A*z_B + z_B**2)
                if D>=0:
                    t1 = (nom_a - np.sqrt(D))/float(denom)
                    t2 = (nom_a + np.sqrt(D))/float(denom)
                    if t1>=0 and t1<=1:
                        res.append([x_A +t1*(x_B-x_A), y_A +t1*(y_B-y_A), z_A +t1*(z_B-z_A)])
                    if t2>=0 and t2<=1:
                        res.append([x_A +t2*(x_B-x_A), y_A +t2*(y_B-y_A), z_A +t2*(z_B-z_A)])
                return res
            
            self.intersections = []
            self.intersections_mirror = []
            for color in ['orange', 'red', 'green', 'blue']:
                for part in self._branches[color]:
                    if len(part)>1:
                        for i in range(0, len(part)-1):
                            self.intersections += getIntersectionWithSphere(part[i], part[i+1])
                for part in self._branches_mirror[color]:
                    if len(part)>1:
                        for i in range(0, len(part)-1):
                            self.intersections_mirror += getIntersectionWithSphere(part[i], part[i+1])
            
            self.printLog('Intersections:')
            self.printLog(str(len(self.intersections)+len(self.intersections_mirror)))
        else:
            self.printLog('Coupler curve must be computed first.')
        

    def getPositionsOfVertices(self, i, solv4, solv5, solv6):
        if self.isComputed():
            v4 = self._trace_v4[solv4][i]
            v5 = self._trace_v5[solv4][solv5][i]
            v6 = self._trace_v6[solv4][solv5][solv6][i]
            v7 = self._trace_v7[i]
            if v4==None or v5==None or v6==None:
                return None
            else:
                return [self._v1, self._v2, self._v3, v4, v5, v6,  v7]
        else:
            self.printLog('Coupler curve must be computed first.')
    
    def getSphere(self,  N):
        u = np.linspace(0, 2 * np.pi, N)
        v = np.linspace(0, np.pi, N)
        x = self.getEdgeLength('16') * np.outer(np.cos(u), np.sin(v))
        y = self.getEdgeLength('16') * np.outer(np.sin(u), np.sin(v))
        z = self.getEdgeLength('16') * np.outer(np.ones(np.size(u)), np.cos(v))
        x = np.add(x, self._v1[0]* np.ones(np.size(u)))
        y = np.add(y, self._v1[1]* np.ones(np.size(u)))
        z = np.add(z, self._v1[2]* np.ones(np.size(u)))
        return (x, y, z)
        
        
    def getEquations(self):
        '''system with correct mixed volume'''
        x4, y4, z4 = symbols('x4 y4 z4')
        x5, y5, z5 = symbols('x5 y5 z5')
        x6, y6, z6 = symbols('x6 y6 z6')
        x7, y7, z7 = symbols('x7 y7 z7')
#        L12 = self.getEdgeLength('12')
#        L13 = self.getEdgeLength('13')
        L14 = self.getEdgeLength('14')
        L15 = self.getEdgeLength('15')
        L16 = self.getEdgeLength('16')
        
#        L27 = self.getEdgeLength('27')
        L37 = self.getEdgeLength('37')
        L47 = self.getEdgeLength('47')
        L57 = self.getEdgeLength('57')
        L67 = self.getEdgeLength('67')
        
#        L23 = self.getEdgeLength('23')
        L34 = self.getEdgeLength('34')
        L45 = self.getEdgeLength('45')
        L56 = self.getEdgeLength('56')
        L26 = self.getR26()
        
        X1, Y1, _ = self._v1
        X2, Y2, _ = self._v2
        X3, Y3, _ = self._v3
        eqs = [
            L26**2 - (Y2 - y6)**2 - x6**2 - z6**2 ,
            L37**2 - Y3**2 - x7**2 - z7**2 ,
            L34**2 - (Y3 - y4)**2 - x4**2 - z4**2 ,
            L14**2 - L34**2 - X1**2 - Y1**2 + Y3**2 + 2*X1*x4 + 2*Y1*y4 - 2*Y3*y4 ,
            L15**2 - (X1 - x5)**2 - (Y1 - y5)**2 - z5**2 ,
            L16**2 - L26**2 - X1**2 - Y1**2 + Y2**2 + 2*X1*x6 + 2*Y1*y6 - 2*Y2*y6 ,
            -L34**2 - L37**2 + L47**2 + 2*Y3**2 + 2*x4*x7 - 2*Y3*y4 + 2*z4*z7 ,
            -L15**2 - L37**2 + L57**2 + X1**2 + Y1**2 + Y3**2 - 2*X1*x5 + 2*x5*x7 - 2*Y1*y5 + 2*z5*z7 ,
            -L15**2 - L34**2 + L45**2 + X1**2 + Y1**2 + Y3**2 - 2*X1*x5 + 2*x4*x5 - 2*Y3*y4 - 2*Y1*y5 + 2*y4*y5 + 2*z4*z5 ,
            -L26**2 - L37**2 + L67**2 + Y2**2 + Y3**2 + 2*x6*x7 - 2*Y2*y6 + 2*z6*z7 ,
            -L15**2 - L26**2 + L56**2 + X1**2 + Y1**2 + Y2**2 - 2*X1*x5 + 2*x5*x6 - 2*Y1*y5 - 2*Y2*y6 + 2*y5*y6 + 2*z5*z6 ,
        ]
        res = []
        for eq in eqs:
            res.append(str(eq)+';')
        return res

    def findEmbeddings(self, tolerance=1.0e-8,  verbose=True):
        syst = self.getEquations()
        if verbose:
            print 'solving system:'
            for pol in syst:
                print pol
        sols = solve(syst,tasks=2)
        result_real = []
        result_complex = []
        if verbose:
            print 'system has', len(sols), 'solutions'
        if len(sols)<48:
            self._window.showError('PHC found only '+str(len(sols))+' solutions!')
        for sol in sols:
            soldic = strsol2dict(sol)
            if is_real(sol, tolerance):
                result_real.append(soldic)
            else:
                result_complex.append(soldic)
        return {'real':result_real, 'complex':result_complex}
        
#    def transformSolution(self, sol):
#        res = copy.deepcopy(sol)
##        res['x6'] = self._x6_V2atOrigin
#        self.printLog(str(res))
#        return res
    
    def getSolutionsForV6(self):    
        real_sol = self.findEmbeddings()['real']
        coordinates = []
        for sol in real_sol:
#            sol_transformed = self.transformSolution(sol)
            coordinates.append([sol['x6'].real, 
                               sol['y6'].real, 
                               sol['z6'].real])
        return coordinates



def getIntersectionOfThreeSpheres(c1,c2,c3,r1,r2,r3):
    if c1==None or c2==None or c3==None:
        return [None, None]
    x2 = float(c2[0]) - float(c1[0])
    y2 = float(c2[1]) - float(c1[1])
    z2 = float(c2[2]) - float(c1[2])
    
    x3 = float(c3[0]) - float(c1[0])
    y3 = float(c3[1]) - float(c1[1])
    z3 = float(c3[2]) - float(c1[2])
    
    d31 = -r1**2 + r3**2 - (x3**2 + y3**2 + z3**2)
    d21 = -r1**2 + r2**2 - (x2**2 + y2**2 + z2**2)
    try:
        z_sqrt = math.sqrt(-d31**2*x2**2 + 2*d21*d31*x2*x3 - d21**2*x3**2 + (4*r1**2*x3**2 - d31**2)*y2**2 
                          - 2*(4*r1**2*x2*x3 - d21*d31)*y2*y3 + (4*r1**2*x2**2 - d21**2)*y3**2 
                          + (4*r1**2*x3**2 + 4*r1**2*y3**2 - d31**2)*z2**2 - 2*(4*r1**2*x2*x3 + 4*r1**2*y2*y3 - d21*d31)*z2*z3 
                          + (4*r1**2*x2**2 + 4*r1**2*y2**2 - d21**2)*z3**2)
        z_a = (1/2.0)*((d31*x2*x3 - d21*x3**2 + d31*y2*y3 - d21*y3**2)*z2 - (d31*x2**2 - d21*x2*x3 + d31*y2**2 - d21*y2*y3)*z3 
                   + z_sqrt*(x3*y2 - x2*y3)
                  )/float(x3**2*y2**2 - 2*x2*x3*y2*y3 + x2**2*y3**2 + (x3**2 + y3**2)*z2**2 
                     - 2*(x2*x3 + y2*y3)*z2*z3 + (x2**2 + y2**2)*z3**2)
        y_a = -1/2.0*(2*x3*z_a*z2 - 2*x2*z_a*z3 - d31*x2 + d21*x3)/(x3*y2 - x2*y3)
        x_a = -1/2.0*(2*y_a*y2 + 2*z_a*z2 + d21)/x2
        
        z_b = 1/2.0*((d31*x2*x3 - d21*x3**2 + d31*y2*y3 - d21*y3**2)*z2 - (d31*x2**2 - d21*x2*x3 + d31*y2**2 - d21*y2*y3)*z3 
                   - z_sqrt*(x3*y2 - x2*y3)
                  )/(x3**2*y2**2 - 2*x2*x3*y2*y3 + x2**2*y3**2 + (x3**2 + y3**2)*z2**2 
                     - 2*(x2*x3 + y2*y3)*z2*z3 + (x2**2 + y2**2)*z3**2)
        y_b = -1/2.0*(2*x3*z_b*z2 - 2*x2*z_b*z3 - d31*x2 + d21*x3)/(x3*y2 - x2*y3)
        x_b = -1/2.0*(2*y_b*y2 + 2*z_b*z2 + d21)/x2
        return  [[x_a + c1[0], y_a + c1[1], z_a + c1[2]],
                [x_b + float(c1[0]), y_b + float(c1[1]), z_b + float(c1[2])]
               ]
    except Exception as e:
        if type(e)==ValueError:
            return [None, None]
 
def dist(u,v):
    return float(np.sqrt( (u[0]-v[0])**2 + (u[1]-v[1])**2 + (u[2]-v[2])**2))
