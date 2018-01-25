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




import numpy as np


import math
#import copy
#from sympy import symbols


from graphEmbedding import *


class GraphCouplerCurve(GraphEmbedding):
    def __init__(self, lengths={}, window=None):
        if not lengths:
            lengths = {'12': 1.99993774567597,
                         '13': 1.99476987780024,
                         '14': 2.003436460984393,
                         '15': 2.00289249524296,
                         '16': 2.000134247468136,
                         '23': 0.999614322089483,
                         '26': 1.001987710974071,
                         '27': 10.53609172287933,
                         '34': 1.00368644488060,
                         '37': 10.53631716364608,
                         '45': 1.001530148504854,
                         '47': 10.53572330314948,
                         '56': 0.995723616535744,
                         '57': 10.53627365999783,
                         '67': 10.53647884635266}
#            lengths = {'67': 6.7082039325, 
#                    '47': 8.0622577483, 
#                    '45': 8.33506448685, 
#                    '56': 15.0953602143, 
#                    '57': 13.0, 
#                    '14': 5.99552332995, 
#                    '15': 14.2478068488, 
#                    '16': 10.0498756211, 
#                    '34': 6.51062976985, 
#                    '12': 1.19993396795, 
#                    '13': 0.894427191, 
#                    '37': 9.2736184955, 
#                    '23': 0.385778797541, 
#                    '27': 9.5852672451,
#                        '26':11.05}
        
        super(GraphCouplerCurve, self).__init__(lengths, 'Max7vertices',  window=window)
        self.updateFixedTriangle()
        self.setRequiresRecomputing()

        self._branches = {}
        self._branches_mirror = {}
        self._samples = 0

    def setLengthsAndUpdateFixedTriangle(self, lengths, setRequiresRecomputing=True):
        self.setLengths(lengths)
        self.updateFixedTriangle()
        if setRequiresRecomputing:
            self.setRequiresRecomputing()

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
        self.setEdgeLength(float(r26), 2, 6)

    def getV1(self):
        return self._fixedTriangle[2]

    def getV2(self):
        return self._fixedTriangle[0]

    def getV3(self):
        return self._fixedTriangle[1]

    def getR26(self):
        return self.getEdgeLength('26')

    def getyV2(self):
        return self.getV2()[1]
    
    def setyV2(self, y2):
        newL23 = self.getEdgeLength(2, 3)-y2+self.getV2()[1]
        self.setEdgeLengthWithCorrespondingOnes(newL23,  [2, 3, 1, 7])
        self.updateFixedTriangle()

    def getPhiDegree(self):
        return self.getPhiRadian()/math.pi*180.0

    def getThetaDegree(self):
        return self.getThetaRadian()/math.pi*180.0

    def setPhiDegree(self, phi):
        self.setPhiRadian(phi*math.pi/180.0)

    def setPhiRadian_old(self, phi):
        v1 = self._fixedTriangle[2]
        y2 = v1[1] - math.tan(phi)*v1[0]
        
        v2 = self._fixedTriangle[0]
        
        newL23 = self.getEdgeLength(2, 3)-y2+v2[1]
        self.setEdgeLengthWithCorrespondingOnes(newL23,  [2, 3, 1, 7])
        self.updateFixedTriangle()
    
    def setPhiRadian(self, phi):
        self.setPhi([2, 3, 1, 7], phi)
        self.updateFixedTriangle()

    def setThetaRadian_old(self, theta):
        l12 = self.getEdgeLength('12')
        l16 = self.getEdgeLength('16')
        r26 = np.sqrt(l12**2+l16**2-2*l12*l16*math.cos(theta))
        self.setEdgeLength(float(r26), 2, 6)
    
    def setThetaRadian(self, theta):
        self.setTheta([2, 1, 6], theta)

    def setThetaDegree(self, phi):
        self.setThetaRadian(phi*math.pi/180.0) 

    
    def getPhiRadian(self):
        try:
            return math.asin((self.getV1()[1]-self.getV2()[1])/float(self.getEdgeLength('12')))
        except:
            self.printLog('Math error in Phi')
            return 0
    
    def getThetaRadian(self):
        try:
            r26 = self.getR26()
            l12 = self.getEdgeLength('12')
            l16 = self.getEdgeLength('16')
            return math.acos((-r26**2+l12**2+l16**2)/float(2*l12*l16))
        except:
            self.printLog('Math error in Theta ')
            return 0

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
        self.updateFixedTriangle()
        self._samples = N
        
        l14 = self.getEdgeLength(1, 4)
        l15 = self.getEdgeLength(1, 5)
        l16 = self.getEdgeLength(1, 6)
        l47 = self.getEdgeLength(4, 7)
        l57 = self.getEdgeLength(5, 7)
        l67 = self.getEdgeLength(6, 7)
        l34 = self.getEdgeLength(3, 4)
        l45 = self.getEdgeLength(4, 5)
        l56 = self.getEdgeLength(5, 6)
        
        v3 = self.getV3()
        v1 = self.getV1()
        
        self.minbound = 0
        self.maxbound = 0
        
        self._trace_v7 = []
        self._trace_v4 = [[], []]
        self._trace_v5 = [[[], []], [[], []]]
        self._trace_v6 = [[[[], []], [[], []]], [[[], []], [[], []]]]
        
        max_x7 = self.getAltitudeAndFoot(2, 3, 7)[0]
        
        for i in range(0,N):
            v7_x = float(max_x7*math.sin(2*math.pi*i/N))
            v7_z = float(max_x7*math.cos(2*math.pi*i/N))
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
                x_C, y_C, z_C = self.getV2()
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
            n = len(self.intersections)+len(self.intersections_mirror)
            self.printLog(str(n))
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
                return [self.getV1(), self.getV2(), self.getV3(), v4, v5, v6,  v7]
        else:
            self.printLog('Coupler curve must be computed first.')
    
#    def getSphere(self,  N):
#        u = np.linspace(0, 2 * np.pi, N)
#        v = np.linspace(0, np.pi, N)
#        x = self.getEdgeLength('16') * np.outer(np.cos(u), np.sin(v))
#        y = self.getEdgeLength('16') * np.outer(np.sin(u), np.sin(v))
#        z = self.getEdgeLength('16') * np.outer(np.ones(np.size(u)), np.cos(v))
#        x = np.add(x, self.getV1()[0]* np.ones(np.size(u)))
#        y = np.add(y, self.getV1()[1]* np.ones(np.size(u)))
#        z = np.add(z, self.getV1()[2]* np.ones(np.size(u)))
#        return (x, y, z)

    def getSolutionsForV6(self, usePrev=True):    
        real_sol = self.findEmbeddings(usePrev=usePrev)['real']
        coordinates = []
        
        for sol in real_sol:
            coordinates.append([sol['x6'].real, 
                               sol['y6'].real, 
                               sol['z6'].real])
        return coordinates

    
    

#    def getSamplingIterator_changingSelf(self, num_phi, num_theta,  l_phi,  r_phi, l_theta,  r_theta):
#        step_phi = (r_phi-l_phi)/float(num_phi)
#        step_theta = (r_theta-l_theta)/float(num_theta)
#
#        def samples():
#            phi = l_phi
#            while (phi<r_phi+step_phi):
#                theta = l_theta
#                while theta<r_theta+step_theta:
#                    self.setPhiRadian(phi)
#                    self.setThetaRadian(theta)
#                    yield self.getLengths()
#                    theta += step_theta
#                phi += step_phi
#        
#        return samples

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
            else:
                self.printLog(str(e))
 
 
#-------------------OLD PARTS----------------------------------------
#    def runSamplingPhiThetaWithMargins_old(self, num_phi, num_theta,  l_phi,  r_phi, l_theta,  r_theta, treshold=0, animate=False):
#        step_phi = (r_phi-l_phi)/float(num_phi)
#        step_theta = (r_theta-l_theta)/float(num_theta)
#        solutions = []
#        for i in range(0, num_phi+1):
#            phi = l_phi+step_phi*i
#            self.setPhiRadian(phi)
#            for j in range(0, num_theta+1):
#                theta = l_theta + step_theta*j
#                self.setThetaRadian(theta)
#                if animate:
#                    self._window.update_graph2phi()
#                    self._window.update_graph2theta()
#                    self.computeIntersections()
#                    self._window.plotScene()
#                num_real = len(self.findEmbeddings(errorMsg=False)['real'])
#                self.printLog(str([phi, theta, num_real]), verbose=1)
#                if num_real>=treshold:
#                    solutions.append([phi, theta, num_real])
#        return solutions




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
