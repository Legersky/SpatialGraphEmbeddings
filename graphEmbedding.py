import numpy as np
#import math
from phcpy.solver import solve
from phcpy.trackers import track
from phcpy.solutions import strsol2dict, is_real
from sympy import symbols


import math


class GraphEmbedding(object):
    def __init__(self, lengths, fixedTriangle, vertexWithFootAtOrigin=None,  window=None):
        self._window = window
#        self._equationsConstructor = equationsConstructor
        self._prevSystem = None
        self._prevSolutions = None
        
        self.setLengths(lengths)
        self._fixedTriangle_vertices = fixedTriangle
        self._vertexWithFootAtOrigin = vertexWithFootAtOrigin
        self.updateFixedTriangle()
    
    def setLengths(self, lengths):
        self._lengths = {}
        try:
            for e in lengths.keys():
                if e[0]=='L':
                    f=e[1:]
                else:
                    f=e
                self.setEdgeLength(float(lengths[e]), int(f[0]), int(f[1]))
                self.printLog(str(f)+': '+str(self.getEdgeLength(int(f[1]), int(f[0]))), verbose=2)
        except KeyError as er:
            self.printLog('Problem with setting lengths: '+str(er))
    
    def getEdgeLength(self, u, v=None):
        if v==None:
            return self.getEdgeLength(int(u[0]), int(u[1]))
        else:
            if u<v:
                return float(self._lengths[(u, v)])
            else:
                return float(self._lengths[(v, u)])
    
    def setEdgeLength(self, Luv, u, v):
        if u<v:
            self._lengths[(u, v)] = Luv
        else:
            self._lengths[(v, u)] = Luv
    
    def getLengths(self):
        return self._lengths

    def printLog(self, s, verbose=0):
        if self._window:
            self._window.printLog(str(s), verbose=verbose)
        else:
            print s

    def getEquations(self):
        return self.constructEquations()
    
    def getAltitudeAndFoot(self, u, v, w):
        ''' Returns altitude of triangle uvw from w and the distance of its foot from v'''
        Luv = self.getEdgeLength(u, v)
        Lvw = self.getEdgeLength(w, v)
        Luw = self.getEdgeLength(w, u)
        cos_alpha = (Lvw**2+Luv**2 - Luw**2)/float(Lvw*2*Luv)
        return [Lvw*math.sin(math.acos(cos_alpha)), cos_alpha*Lvw]
    
    def coordinatesOfTriangle(self, u, v, w, yshift=0):
        '''Returns coordinates of the tringle uvw so that it lies in x-y plane, u,v are on y-axis, y-coord. of u is yshift and v is in positive direction from u'''
        u_coor = [0, yshift, 0]
        v_coor = [0, yshift+self.getEdgeLength(u, v), 0]
        alt_w, foot_wv = self.getAltitudeAndFoot(u, v, w)
        w_coor = [alt_w, v_coor[1]-foot_wv, 0]
        return u_coor, v_coor, w_coor
    
    def updateFixedTriangle(self):
        '''Adjusts coordinates of the fixed triangle according to _lengths. If p!=None, the coordinate system is shifted so that foot of altitude from p in the triangle uvp is in the origin. (uvp must be in the graph)'''
        u, v, w = self._fixedTriangle_vertices
        if self._vertexWithFootAtOrigin!=None:
            yshift = -self.getAltitudeAndFoot(v, u, self._vertexWithFootAtOrigin)[1]
        else:
            yshift = 0
        self._fixedTriangle = self.coordinatesOfTriangle(u, v, w, yshift)
    
    def setEdgeLengthWithCorrespondingOnes(self, Luv,  u, v, w, p):
        '''Sets length of uv to Luv and also lengths of uw and up so that angles uvp and uvw preserves '''
        p_coord = self.coordinatesOfTriangle(u, v, p)[2]
        w_coord = self.coordinatesOfTriangle(u, v, w)[2]
        new_u = [0, self.getEdgeLength(u, v)-Luv, 0]
        self.setEdgeLength(Luv, u, v)
        self.setEdgeLength(self.dist(p_coord, new_u), p, u)
        self.setEdgeLength(self.dist(w_coord, new_u), w, u)

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
    
#    def getLengthsForPhiTheta(self, phi, theta):
##    def setThetaRadian(self, theta):
##        l12 = self.getEdgeLength('12')
##        l16 = self.getEdgeLength('16')
##        r26 = np.sqrt(l12**2+l16**2-2*l12*l16*math.cos(theta))
##        self.setR26(r26)     
#        
##        v1, _, v3, max_x7 = self.getFixedTriangleAndMaxX7(_lengths)
##        lengths = copy.copy(_lengths)     
##        x1, y1, z1 = v1
##        y2 = y1 - math.tan(phi)*x1
##
##        lengths['L12'] = self.dist(v1, [0, y2, 0])
##        lengths['L27'] = self.dist([max_x7, 0, 0], [0, y2, 0])
##        lengths['L23'] = self.dist(v3, [0, y2, 0])
##        
##        l12 = lengths['L12']
##        l16 = lengths['L16']
##        lengths['L26'] = np.sqrt(l12**2+l16**2-2*l12*l16*math.cos(theta))
##        return lengths
#
#        v1, v2, v3, max_x7 = self.getFixedTriangleAndMaxX7(self._lengths)
#        x1, y1, z1 = v1
#        y2 = y1 - math.tan(phi)*x1
#
#        v2 = [0, y2, 0]
#        
#        self._lengths['L12'] = self.dist(v1, v2)
#        self._lengths['L27'] = self.dist([max_x7, 0, 0], v2)
#        self._lengths['L23'] = self.dist(v3, v2)
#
#        self.setThetaRadian(theta)
#        return self.getLengths()

    def getFixedTriangleAndMaxX7(self):
        l12 = self._lengths['L12']
        l13 = self._lengths['L13']
        l27 = self._lengths['L27']
        l37 = self._lengths['L37']
        l23 = self._lengths['L23']
        
        theta1 = math.acos((-l13**2+l12**2+l23**2)/(2*l12*l23))
        theta7 = math.acos((-l37**2+l27**2+l23**2)/(2*l27*l23))

        max_x7 = math.sin(theta7)*l27
        x1 = math.sin(theta1)*l12

        y7 = math.cos(theta7)*l27
        y1 = math.cos(theta1)*l12
       
        v1 = [x1, y1-y7, 0]
        v2 = [0, -y7, 0]
        v3 = [0, -y7+l23, 0]
        
        return [v1, v2, v3, max_x7]

#    def getPhiTheta(self, lengths):
#        v1, v2, v3, max_x7 = self.getFixedTriangleAndMaxX7(lengths)
#        y1 = v1[1]
#        y2 = v2[1]
#        l26 = lengths['L26']
#        l12 = lengths['L12']
#        l16 = lengths['L16']
#        return [math.asin((y1-y2)/float(lengths['L12'])), math.acos((-l26**2+l12**2+l16**2)/float(2*l12*l16))]

    def constructEquations(self):
        '''system with correct mixed volume'''
        x4, y4, z4 = symbols('x4 y4 z4')
        x5, y5, z5 = symbols('x5 y5 z5')
        x6, y6, z6 = symbols('x6 y6 z6')
        x7, y7, z7 = symbols('x7 y7 z7')
#        L12 = lengths['12']
#        L13 = lengths['13']
        L14 = self.getEdgeLength('14')
        L15 = self.getEdgeLength('15')
        L16 = self.getEdgeLength('16')
    #        L27 = lengths['27']
        L37 = self.getEdgeLength('37')
        L47 = self.getEdgeLength('47')
        L57 = self.getEdgeLength('57')
        L67 = self.getEdgeLength('67')
        
#        L23 = lengths['23']
        L34 = self.getEdgeLength('34')
        L45 = self.getEdgeLength('45')
        L56 = self.getEdgeLength('56')
        L26 = self.getEdgeLength('26')
               
        X1, Y1, _ = self._fixedTriangle[2]
        _, Y2, _ = self._fixedTriangle[0]
        _, Y3, _ = self._fixedTriangle[1]

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








#        #        u = [0, y2, 0]
#        newL23 = self.getEdgeLength(2, 3)+y2-self.getV2()[1]
#        self.setEdgeLengthWithCorrespondingOnes(newL23,  2, 3, 1, 7)
#        self.updateFixedTriangle()
#        self._lengths[(1, 2)] = self.dist(self.getV1(), self.getV2())
#        self._lengths[(2, 7)] = self.dist([self._max_x7, 0, 0], self.getV2())
#        self._lengths[(2, 3)] = self.dist(self.getV3(), self.getV2())
        
#    def updateFixedTriangle(self):       
#        l12 = self.getEdgeLength('12')
#        l13 = self.getEdgeLength('13')
#        l27 = self.getEdgeLength('27')
#        l37 = self.getEdgeLength('37')
#        l23 = self.getEdgeLength('23')
#        
#        theta1 = math.acos((-l13**2+l12**2+l23**2)/(2*l12*l23))
#        theta7 = math.acos((-l37**2+l27**2+l23**2)/(2*l27*l23))
#
#        self._max_x7 = math.sin(theta7)*l27
#        x1 = math.sin(theta1)*l12
#
#        y7 = math.cos(theta7)*l27
#        y1 = math.cos(theta1)*l12
#
#        y2 = -y7
#        y3 = -y7+l23
#        y1 = y1-y7
#        
#        self._v1 = [x1, y1, 0]
#        self._v2 = [0, y2, 0]
#        self._v3 = [0, y3, 0]
