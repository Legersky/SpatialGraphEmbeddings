import numpy as np
import math
from phcpy.solver import solve
from phcpy.trackers import track
from phcpy.solutions import strsol2dict, is_real
from sympy import symbols
import time
import copy
from sklearn.cluster import DBSCAN

from graphEmbedding import *


class GraphEmbeddingVangelis(GraphEmbedding):
    def __init__(self,lengths={}, r26=0.550000000000000,   window=None):

        if not lengths:
            lengths = {'L12': 0.956581258335196,
                     'L13': 2.651661708696340,
                     'L14': 3.74969777489766,
                     'L15': 3.076648111148685,
                     'L16': 1.101294935802410,
                     'L23': 2.800000000000000,
                     'L27': 1.526854929507074,
                     'L34': 4.49999748323629,
                     'L37': 3.12265521649798,
                     'L45': 6.45000000000000,
                     'L47': 2.705062963547708,
                     'L56': 3.407714967884632,
                     'L57': 4.13929267237052,
                     'L67': 1.218419216365267}

        if not 'L26' in lengths.keys():
            lengths['L26'] = r26
        
        super(GraphEmbeddingVangelis, self).__init__(lengths, self.constructEquations, window=window)
        
        self.updateFixedTriangle()
        
        self._branches = {}
        self._branches_mirror = {}
        self._samples = 0
    
    def setLengthsAndUpdateFixedTriangle(self, lengths):
        self.setLengths(lengths)
        self.updateFixedTriangle()

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

    def getSamples(self):
        return self._samples

    def setR26(self, r26):
        self._lengths['L26'] = float(r26)
    
    def getR26(self):
        return self.getEdgeLength('26')

    def getyV2(self):
        return self._v2[1]
    
    def getPhiRadian(self):
#        try:
            return math.asin((self._v1[1]-self._v2[1])/float(self.getEdgeLength('12')))
#        except:
#            self.printLog('Math error in Phi')
#            return 0
    
    def getPhiDegree(self):
        return self.getPhiRadian()/math.pi*180.0
    
    def getThetaRadian(self):
#        try:
            r26 = self.getR26()
            l12 = self.getEdgeLength('12')
            l16 = self.getEdgeLength('16')
            return math.acos((-r26**2+l12**2+l16**2)/float(2*l12*l16))
#        except:
#            self.printLog('Math error in Theta ')
#            return 0

    def getThetaDegree(self):
        return self.getThetaRadian()/math.pi*180.0

    def updateFixedTriangle(self):       
        l12 = self.getEdgeLength('12')
        l13 = self.getEdgeLength('13')
        l27 = self.getEdgeLength('27')
        l37 = self.getEdgeLength('37')
        l23 = self.getEdgeLength('23')
        
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
    
    def setPhiRadian(self, phi):
        y2 = self._v1[1] - math.tan(phi)*self._v1[0]
        self.setyV2(y2)
    
    def setPhiDegree(self, phi):
        self.setPhiRadian(phi*math.pi/180.0)

    def setThetaRadian(self, theta):
        l12 = self.getEdgeLength('12')
        l16 = self.getEdgeLength('16')
        r26 = np.sqrt(l12**2+l16**2-2*l12*l16*math.cos(theta))
        self.setR26(r26)
    
    def setThetaDegree(self, phi):
        self.setThetaRadian(phi*math.pi/180.0)
    
    def setyV2(self, y2):
        self._v2 = [0, y2, 0]
        
        self._lengths['L12'] = self.dist(self._v1, self._v2)
        self._lengths['L27'] = self.dist([self._max_x7, 0, 0], self._v2)
        self._lengths['L23'] = self.dist(self._v3, self._v2)

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
        
        self._trace_v7 = []
        self._trace_v4 = [[], []]
        self._trace_v5 = [[[], []], [[], []]]
        self._trace_v6 = [[[[], []], [[], []]], [[[], []], [[], []]]]
        
        for i in range(0,N):
            v7_x = float(self._max_x7*math.sin(2*math.pi*i/N))
            v7_z = float(self._max_x7*math.cos(2*math.pi*i/N))
            v7 = [v7_x, 0, v7_z]
            self._trace_v7.append(v7)
            
            intersection_v4 = self.getIntersectionOfThreeSpheres(v3,v1,v7,l34,l14,l47)
            v4 = intersection_v4[0]
            self._trace_v4[0].append(v4)
            if True:
                intersection_v5 = self.getIntersectionOfThreeSpheres(v4,v1,v7,l45,l15,l57)
                v5 = intersection_v5[0]
                self._trace_v5[0][0].append(v5)
                if True:
                    intersection_v6 = self.getIntersectionOfThreeSpheres(v5,v1,v7,l56,l16,l67)
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
                    intersection_v6 = self.getIntersectionOfThreeSpheres(v5,v1,v7,l56,l16,l67)
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
                intersection_v5 = self.getIntersectionOfThreeSpheres(v4,v1,v7,l45,l15,l57)
                v5 = intersection_v5[0]
                self._trace_v5[1][0].append(v5)
                if True:
                    intersection_v6 = self.getIntersectionOfThreeSpheres(v5,v1,v7,l56,l16,l67)
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
                    intersection_v6 = self.getIntersectionOfThreeSpheres(v5,v1,v7,l56,l16,l67)
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
                r = self.getR26()
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
    
#    def getSphere(self,  N):
#        u = np.linspace(0, 2 * np.pi, N)
#        v = np.linspace(0, np.pi, N)
#        x = self.getEdgeLength('16') * np.outer(np.cos(u), np.sin(v))
#        y = self.getEdgeLength('16') * np.outer(np.sin(u), np.sin(v))
#        z = self.getEdgeLength('16') * np.outer(np.ones(np.size(u)), np.cos(v))
#        x = np.add(x, self._v1[0]* np.ones(np.size(u)))
#        y = np.add(y, self._v1[1]* np.ones(np.size(u)))
#        z = np.add(z, self._v1[2]* np.ones(np.size(u)))
#        return (x, y, z)

    def constructEquations(self):
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

    def getSolutionsForV6(self, usePrev=True):    
        real_sol = self.findEmbeddings(usePrev=usePrev)['real']
        coordinates = []
        
        for sol in real_sol:
            coordinates.append([sol['x6'].real, 
                               sol['y6'].real, 
                               sol['z6'].real])
        return coordinates

    
    def runSamplingPhiTheta(self, num_phi, num_theta):
        start = time.time()
        act_num = len(self.findEmbeddings()['real'])
        act_phi = self.getPhiRadian()
        act_theta = self.getThetaRadian()
       
        margin_degree = 5
        margin = margin_degree*math.pi/180.0
        l_phi,  r_phi = -math.pi/2.0 + margin, math.pi/2.0 - margin
        l_theta, r_theta = 0 + margin/2.0, math.pi - margin/2.0
        sols = self.runSamplingPhiThetaWithMargins(num_phi, num_theta,  l_phi,  r_phi, l_theta,  r_theta,  treshold=act_num)
        sols.append([act_phi, act_theta, act_num])
        end = time.time()
        self.printLog('time 1st round: '+str(end - start))
        self.printLog('Maximum number of embeddings in 1st round:')

        maximum = max([num for phi, theta, num in sols])
        self.printLog(str(maximum))
        max_positions = [ [phi, theta, num] for phi, theta, num in sols if num==maximum]
        step_phi = (r_phi-l_phi)/float(2*num_phi)
        step_theta = (r_theta-l_theta)/float(2*num_theta)
        argmax = []
        
        for phi, theta, num in max_positions:
            if num>=maximum:
                argmax.append([phi, theta])
            for phi2, theta2, num2 in self.runSamplingPhiThetaWithMargins(num_phi/4, num_theta/4,
                                                                          phi - step_phi,  phi + step_phi,
                                                                          theta - step_theta,  theta + step_theta,
                                                                          treshold=maximum):
                if num2>maximum:
                    self.setPhiRadian(phi2)
                    self.setThetaRadian(theta2)
                    self.printLog(str(num2)+',trying increase max', verbose=2)
                    num_check = len(self.findEmbeddings(usePrev=False)['real'])
                    if num_check>maximum:
                        maximum = num_check
                        argmax = []
                        self.printLog('Maximum increased to '+str(num_check), verbose=1)
                    argmax.append([phi2, theta2])
                elif num2==maximum:
                    argmax.append([phi2, theta2])
                    self.printLog(str(num2), verbose=2)
#        comparing minimums was here

        eps = 0.1
        for i in range(0, 100):
            labels = DBSCAN(eps=eps).fit_predict(argmax)
            if len([1 for el in labels if el==-1])>len(labels)/10:
                eps += 0.05
                self.printLog('Increasing eps for clustering')
            else:
                break
            if i==99:
                self.printLog('Clustering was not succesfull')
                labels = [0 for el in labels]
        
        clusters = [[] for i in range(0, max(list(labels)+[0])+1)]
        
        for label, point in zip(labels, argmax):
            if label>=0:
                clusters[label].append(point)
        
        print clusters
        
        res_graphs= []
        res_infos = []
        centers = []
        
        def dist_angle(a, b):
            return (a[0]-b[0])**2 + (a[1]-b[1])**2
        
        for cluster in clusters:
            s_x = 0
            s_y = 0
            for x, y in cluster:
                s_x += x
                s_y += y
            phi, theta = s_x/float(len(cluster)), s_y/float(len(cluster))
            self.setPhiRadian(phi)
            self.setThetaRadian(theta)
            n = len(self.findEmbeddings()['real'])
            if n < maximum:
                min_dist = dist_angle([phi, theta], cluster[0])
                self.setPhiRadian(cluster[0][0])
                self.setThetaRadian(cluster[0][1])
                for x, y in cluster:
                    d = dist_angle([phi, theta], [x, y])
                    if d < min_dist:
                        min_dist = d
                        self.setPhiRadian(x)
                        self.setThetaRadian(y)
                self.printLog('Center of cluster does not have maximum number of solutions \n -> nearest point chosen instead.')
            
            centers.append([self.getPhiRadian(), self.getThetaRadian()])
            res_graphs.append(GraphEmbeddingVangelis(lengths=self._lengths,  r26=self.getR26(), window=self._window))
            res_infos.append(str([self.getPhiRadian(), self.getThetaRadian()]))
        
        if maximum<act_num:
            self.setPhiRadian(act_phi)
            self.setThetaRadian(act_theta)
            res_graphs= [GraphEmbedding(lengths=self._lengths,  r26=self.getR26(), window=self._window)]
            res_infos = ['original kept:\n'+str([self.getPhiRadian(), self.getThetaRadian()])]
            centers = [[self.getPhiRadian(), self.getThetaRadian()]]
            sols.append([act_phi, act_theta, act_num])
            maximum = act_num
        
        if self._window:
            self._window.showClusters(clusters, centers)
            self._window.setGraphSequence(res_graphs, res_infos)

        end = time.time()
        self.printLog('Maximum number of embeddings:')
        self.printLog(str(maximum))

        self.printLog('time: '+str(end - start))
        
        return centers
    
    def getSamplingIterator(self, num_phi, num_theta,  l_phi,  r_phi, l_theta,  r_theta):
        step_phi = (r_phi-l_phi)/float(num_phi)
        step_theta = (r_theta-l_theta)/float(num_theta)

        def samples():
            phi = l_phi
            while (phi<r_phi+step_phi):
                theta = l_theta
                while theta<r_theta+step_theta:
                    yield [phi, theta]
                    theta += step_theta
                phi += step_phi
        
        return samples
    
    def runSamplingPhiThetaWithMargins(self, num_phi, num_theta,  l_phi,  r_phi, l_theta,  r_theta, treshold=0, animate=False):
        return self.runSampling(self.getSamplingIterator(num_phi, num_theta,  l_phi,  r_phi, l_theta,  r_theta), treshold=treshold, animate=animate)

    def runSampling(self, iterator, treshold=0, animate=False):
        solutions = []
        for phi, theta in iterator():
            self.setPhiRadian(phi)
            self.setThetaRadian(theta)
            if animate:
                self._window.update_graph2phi()
                self._window.update_graph2theta()
                self.computeIntersections()
                self._window.plotScene()
            num_real = len(self.findEmbeddings(errorMsg=False)['real'])
            self.printLog(str([phi, theta, num_real]), verbose=1)
            if num_real>=treshold:
                solutions.append([phi, theta, num_real])
        return solutions

    def runSamplingPhiThetaWithMargins_old(self, num_phi, num_theta,  l_phi,  r_phi, l_theta,  r_theta, treshold=0, animate=False):
        step_phi = (r_phi-l_phi)/float(num_phi)
        step_theta = (r_theta-l_theta)/float(num_theta)
        solutions = []
        for i in range(0, num_phi+1):
            phi = l_phi+step_phi*i
            self.setPhiRadian(phi)
            for j in range(0, num_theta+1):
                theta = l_theta + step_theta*j
                self.setThetaRadian(theta)
                if animate:
                    self._window.update_graph2phi()
                    self._window.update_graph2theta()
                    self.computeIntersections()
                    self._window.plotScene()
                num_real = len(self.findEmbeddings(errorMsg=False)['real'])
                self.printLog(str([phi, theta, num_real]), verbose=1)
                if num_real>=treshold:
                    solutions.append([phi, theta, num_real])
        return solutions

    def getIntersectionOfThreeSpheres(self, c1,c2,c3,r1,r2,r3):
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
 
#-------------------OLD PARTS----------------------------------------
#        max_minimums = []
#        res_graphs = []
#        res_infos = []
#        winner = ['dist', 'dist v6', 'dist 5', 'dist 10', 'dist 100', 'imag']+ ['avg '+ s for s in ['dist', 'dist v6', 'dist 5', 'dist 10', 'dist 100']] + ['length change']
#        for w in winner:
#            max_minimums.append(0)
#            res_graphs.append(None)
#            res_infos.append('')
#        for graph in graph_seq:
#            mins = graph.findMinDistanceOfSolutions(orig_l23, orig_r26)
#            info = ''
#            for k,  minimum in enumerate(mins):
#                info +=' min'+winner[k]+': '+str(minimum)+'\n'  

#            for k, minimum in enumerate(mins):
#                if max_minimums[k] < minimum:
#                    max_minimums[k] = minimum
#                    res_graphs[k] = graph
#                    res_infos[k] = 'Winner '+winner[k]+'\n'+ info


#    def findMinDistanceOfSolutions(self, orig_l23, orig_r26):
#        sols = self.findEmbeddings()
#        sols_real = sols['real']
#        min_imag = 1.0e20
#
#        combinations = []
#        for i in range(0,len(sols_real)):
#            for j in range(i+1,len(sols_real)):
#                combinations.append([forgetComplexPart(sols_real[i]), forgetComplexPart(sols_real[j])])
#        
#        
#        d = sorted([dist_solutions(sol1, sol2) for sol1, sol2 in combinations])
#        d_v6 = sorted([dist_v6(sol1, sol2) for sol1, sol2 in combinations])
#        d5 = sorted([dist_solutions_weigthed(sol1, sol2, 5) for sol1, sol2 in combinations])
#        d10 = sorted([dist_solutions_weigthed(sol1, sol2, 10) for sol1, sol2 in combinations])
#        d100 = sorted([dist_solutions_weigthed(sol1, sol2, 100) for sol1, sol2 in combinations])
#
#        
#        for sol in sols['complex']:
#            s_imag = sumImaginaryPart(sol)
#            if s_imag < min_imag:
#                min_imag = s_imag
#
#        def avgComparable(d, ratio=0.2):
#            res = 0
#            L = len(d)
#            for i in range(0, L-1):
#                res += d[i]
#                if d[0]/(d[i+1]+0.00001) < ratio:
#                    return res/float(i+1)
#    
#        return (d[0], d_v6[0], d5[0], d10[0],  d100[0],  min_imag,
#                avgComparable(d), avgComparable(d_v6), avgComparable(d5), avgComparable(d10),  avgComparable(d100),
#                math.exp(-(abs(self.getEdgeLength('23')-orig_l23)+ abs(self.getR26()-orig_r26)))
#                )
#
#def dist_solutions(sol1, sol2):
#    res = 0
#    for k in ['z4', 'y4', 'x4', 'y5', 'x5', 'z5', 'x6', 'z6', 'y6', 'x7', 'z7']:
#        res += (sol1[k] - sol2[k])**2
#    return np.sqrt(res)
#
#def dist_solutions_weigthed(sol1, sol2, weight):
#    res = 0
#    for k in ['z4', 'y4', 'x4', 'y5', 'x5', 'z5', 'x7', 'z7']:
#        res += (sol1[k] - sol2[k])**2
#    for k in ['x6', 'z6', 'y6']:
#        res += (weight*(sol1[k] - sol2[k]))**2
#    return np.sqrt(res)
#
#def sumImaginaryPart(sol):
#    res = 0
#    for k in ['z4', 'y4', 'x4', 'y5', 'x5', 'z5', 'x6', 'z6', 'y6', 'x7', 'z7']:
#        res += abs(sol[k].imag)
#    return res
#
#def dist_v6(sol1, sol2):
#    return self.dist([sol1[k] for k in ['x6', 'z6', 'y6']], [sol2[k] for k in ['x6', 'z6', 'y6']])
    
#def forgetComplexPart(sol):
#    res = {}
#    for k in ['z4', 'y4', 'x4', 'y5', 'x5', 'z5', 'x6', 'z6', 'y6', 'x7', 'z7']:
#        res[k] = sol[k].real
#    return res
