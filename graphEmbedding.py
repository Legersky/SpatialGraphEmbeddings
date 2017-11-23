import numpy as np
#import math
#from phcpy.solver import solve
#from phcpy.trackers import track
from phcpy.solutions import strsol2dict, is_real
    
from sympy import symbols

from random import random
import hashlib
import time
import os
import ast
import subprocess

import math


class GraphEmbedding(object):
    def __init__(self, lengths, graph_type, tmpFileName=None,  window=None):
        self._window = window
#        self._equationsConstructor = equationsConstructor
        self._prevSystem = None
        self._prevSolutions = None
        self._verbose = 1
        self._distSystem = False
        if graph_type == 'Max7vertices':
            self._fixedTriangle_vertices = [2, 3, 1]
            self._vertexWithFootAtOrigin = 7
            self.constructEquations = self.constructEquations_max7vertices
            self._numAllSolutions = 48
        elif graph_type == 'Max6vertices':
            self._fixedTriangle_vertices = [2, 3, 1]
            self._vertexWithFootAtOrigin = 4
            self.constructEquations = self.constructEquations_max6vertices
            self._numAllSolutions = 16
        elif graph_type == 'Max8vertices':
            self._fixedTriangle_vertices = [2, 3, 1]
            self._vertexWithFootAtOrigin = None
            self.constructEquations = self.constructEquations_max8vertices
            self._numAllSolutions = 160
        elif graph_type == 'Ring8vertices':
            self._fixedTriangle_vertices = [2, 3, 1]
            self._vertexWithFootAtOrigin = None
            self.constructEquations = self.constructEquations_ring8vertices
            self._numAllSolutions = 128
        elif graph_type == '7vert32a':
            self._fixedTriangle_vertices = [2, 3, 1]
            self._vertexWithFootAtOrigin = None
            self.constructEquations = self.constructEquations_7vert32a
            self._numAllSolutions = 32
        elif graph_type == '7vert32b':
            self._fixedTriangle_vertices = [2, 3, 1]
            self._vertexWithFootAtOrigin = None
            self.constructEquations = self.constructEquations_7vert32b
            self._numAllSolutions = 32
        elif graph_type == '7vert24':
            self._fixedTriangle_vertices = [2, 3, 1]
            self._vertexWithFootAtOrigin = None
            self.constructEquations = self.constructEquations_7vert24
            self._numAllSolutions = 24
        elif graph_type == '7vert16a':
            self._fixedTriangle_vertices = [2, 3, 1]
            self._vertexWithFootAtOrigin = None
            self.constructEquations = self.constructEquations_7vert16a
            self._numAllSolutions = 16
        elif graph_type == '7vert16b':
            self._fixedTriangle_vertices = [2, 3, 1]
            self._vertexWithFootAtOrigin = None
            self.constructEquations = self.constructEquations_7vert16b
            self._numAllSolutions = 16
        else:
            raise ValueError('Type %s not supported' % graph_type)
        self._graph_type = graph_type
        if tmpFileName==None:
            hash_object = hashlib.md5(str(time.time()).encode()+str(random()))
            self._fileNamePref = graph_type+'_'+str(hash_object.hexdigest())
        else:
            self._fileNamePref = graph_type+'_'+str(tmpFileName)
        self.setLengths(lengths)    

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
            self.updateFixedTriangle()
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
        if Luv<=0:
            raise ValueError('Length of '+str([u, v])+' cannot be set to '+str(Luv))
        if Luv>1e6:
            self._window.showError('Length'+str(Luv)+ ' of '+str([u, v])+' is too big')
            raise ValueError('Length'+str(Luv)+ ' of '+str([u, v])+' is too big')
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
            if verbose<= self._verbose:
                print s

    def getEquations(self):
        return self.constructEquations()
    
    def getAltitudeAndFoot(self, u, v, w):
        ''' Returns altitude of triangle uvw from w and the distance of its foot from v'''
        Luv = self.getEdgeLength(u, v)
        Lvw = self.getEdgeLength(w, v)
        Luw = self.getEdgeLength(w, u)
        cos_alpha = (Lvw**2+Luv**2 - Luw**2)/float(Lvw*2*Luv)
        if cos_alpha>=-1 and cos_alpha<=1:
            return [Lvw*math.sin(math.acos(cos_alpha)), cos_alpha*Lvw]
        else:
            raise ValueError('Altitude and foot for the triangle '+str([u, v, w])+' with lengths '+str([Luv, Lvw, Luw])+' is not defined.')
    
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
    
    def setEdgeLengthWithCorrespondingOnes(self, Luv,  uvwp):
        '''Sets length of uv to Luv and also lengths of uw and up so that angles uvp and uvw preserves '''
        u, v, w, p = uvwp
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

        if self._distSystem:
            y1_
        while True:
            if self._prevSystem and usePrev:
                filenameTmp = 'tmp/'+self._fileNamePref+'_prev.txt'
                with open(filenameTmp, 'w') as filePrev:
                    filePrev.write(str(self._prevSystem)+'\n')
                    filePrev.write(str(self._prevSolutions)+'\n')
                    
                process = subprocess.Popen(['python2','track.py', str(syst), filenameTmp, self._fileNamePref])
                process.wait()
                with open('tmp/'+self._fileNamePref+'track.txt','r') as file:
                    solutions_str = file.read()
                sols = ast.literal_eval(solutions_str)
                
                os.remove(filenameTmp)
                os.remove('tmp/'+self._fileNamePref+'track.txt')

            else:
                process = subprocess.Popen(['python2','solve.py', str(syst), self._fileNamePref])
                process.wait()
                with open('tmp/'+self._fileNamePref+'solve.txt','r') as file:
                    sols_str = file.read()
                sols = ast.literal_eval(sols_str)
                os.remove('tmp/'+self._fileNamePref+'solve.txt')

            result_real = []
            result_complex = []

            for sol in sols:
                soldic = strsol2dict(sol)
                if is_real(sol, tolerance):
                    result_real.append(soldic)
                else:
                    result_complex.append(soldic)

            if len(result_real)%4==0 and len(sols)==self._numAllSolutions:
                self._prevSystem = syst
                self._prevSolutions = sols
                return {'real':result_real, 'complex':result_complex}
            else:
                usePrev = False
                i += 1
                self.printLog('PHC failed, trying again: '+str(i), verbose=1)
            if i>=10:
                if self._window and errorMsg:
                    self._window.showError('PHC failed 10 times')
                print 'PHC failed 10 times'
                return {'real':[], 'complex':[]}

    def setPhi(self, uvwp, phi):
        margin_degree = 4
        margin = margin_degree*math.pi/180.0
        l_phi,  r_phi = -math.pi/2.0 + margin, math.pi/2.0 - margin
        if phi<=l_phi:
            self.printLog('Phi '+str(phi)+' too small -> set to '+str(l_phi))
            phi = l_phi
        if phi>=r_phi:
            self.printLog('Phi '+str(phi)+' too large -> set to '+str(r_phi))
            phi = r_phi
        u, v, w, p = uvwp
        altitude, foot_v = self.getAltitudeAndFoot(u, v, w)
        a = altitude*math.tan(phi)
        self.setEdgeLengthWithCorrespondingOnes(abs(foot_v+a), uvwp)

    def setTheta(self, uwc, theta):
        u, w, c = uwc
        Luw = self.getEdgeLength(u, w)
        Lwc = self.getEdgeLength(c, w)
        self.setEdgeLength(np.sqrt(Luw**2+Lwc**2-2*Luw*Lwc*math.cos(theta)), u, c)

    def setPhiTheta(self, uvwpc, phi, theta):
        self.setPhi(uvwpc[:-1], phi)
        self.setTheta([uvwpc[i] for i in [0, 2, 4]], theta)

    def getPhiTheta(self, uvwpc):      
        u, v, w, p, c = uvwpc
        
        foot_u = self.getAltitudeAndFoot(v, u, w)[1]
#        altitude_w = self.getAltitudeAndFoot(u, v, w)[0]
        Luw = self.getEdgeLength(u, w)
        Luc = self.getEdgeLength(u, c)
        Lcw = self.getEdgeLength(w, c)
#        return [math.acos(altitude_w/float(Luw)), math.acos((Luw**2+Lcw**2-Luc**2)/float(2*Luw*Lcw))]
        return [math.asin(foot_u/float(Luw)), math.acos((Luw**2+Lcw**2-Luc**2)/float(2*Luw*Lcw))]

    def getPhiTheta_old(self):
        v1, v2, v3 = self.coordinatesOfTriangle(1, 2, 3)
        y1 = v1[1]
        y2 = v2[1]
        l26 = self.getEdgeLength(2, 6)
        l12 = self.getEdgeLength(1, 2)
        l16 = self.getEdgeLength(1, 6)
        return [math.asin((y1-y2)/float(self.getEdgeLength(1, 2))), math.acos((-l26**2+l12**2+l16**2)/float(2*l12*l16))]

    def getEmbedding(self):
        sols = self.findEmbeddings()
        sol = sols['real'][0]
        if self._graph_type=='Max6vertices':
            res = [[0, 0, 0] for i in range(0, 6)]
            for k in sol:
                if k[0]=='x':
                    res[int(k[1])-1][0] = sol[k].real
                if k[0]=='y':
                    res[int(k[1])-1][1] = sol[k].real
                if k[0]=='z':
                    res[int(k[1])-1][2] = sol[k].real
            yshift = -self.getAltitudeAndFoot(3, 2, 4)[1]
            v2, v3, v1 = self.coordinatesOfTriangle(2, 3, 1, yshift)
            res[0] = v1
            res[1] = v2
            res[2] = v3
            return res
        elif self._graph_type=='Max7vertices':
            res = [[0, 0, 0] for i in range(0, 7)]
            for k in sol:
                if k[0]=='x':
                    res[int(k[1])-1][0] = sol[k].real
                if k[0]=='y':
                    res[int(k[1])-1][1] = sol[k].real
                if k[0]=='z':
                    res[int(k[1])-1][2] = sol[k].real
            yshift = -self.getAltitudeAndFoot(3, 2, 7)[1]
            v2, v3, v1 = self.coordinatesOfTriangle(2, 3, 1, yshift)
            res[0] = v1
            res[1] = v2
            res[2] = v3
            return res
        elif self._graph_type=='Max8vertices':
            res = [[0, 0, 0] for i in range(0, 8)]
            for k in sol:
                if k[0]=='x':
                    res[int(k[1])-1][0] = sol[k].real
                if k[0]=='y':
                    res[int(k[1])-1][1] = sol[k].real
                if k[0]=='z':
                    res[int(k[1])-1][2] = sol[k].real
            v2, v3, v1 = self.coordinatesOfTriangle(2, 3, 1)
            res[0] = v1
            res[1] = v2
            res[2] = v3
            return res

    def constructEquations_max7vertices(self):
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
        
#        self._vertexWithFootAtOrigin = 7
#        self.updateFixedTriangle()
#        X1, Y1, _ = self._fixedTriangle[2]
#        _, Y2, _ = self._fixedTriangle[0]
#        _, Y3, _ = self._fixedTriangle[1]
        try:
            yshift = -self.getAltitudeAndFoot(3, 2, 7)[1]
            v2, v3, v1 = self.coordinatesOfTriangle(2, 3, 1, yshift)
        except:
            raise TriangleInequalityError('Triangle inequality violated!')
        
        X1, Y1, _ = v1
        _, Y2, _ = v2
        _, Y3, _ = v3
        
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

    def constructEquations_max8vertices(self):
        '''system with correct mixed volume'''        
#        x1, y1, z1 = symbols('x1 y1 z1')
#        x3, y3, z3 = symbols('x3 y3 z3')
        x4, y4, z4 = symbols('x4 y4 z4')
        x5, y5, z5 = symbols('x5 y5 z5')
        x6, y6, z6 = symbols('x6 y6 z6')
        x7, y7, z7 = symbols('x7 y7 z7')
        x8, y8, z8 = symbols('x8 y8 z8')
        
        L12 = self.getEdgeLength('12')
        L13 = self.getEdgeLength('13')
        L14 = self.getEdgeLength('14')
        L15 = self.getEdgeLength('15')
        L16 = self.getEdgeLength('16')

        L27 = self.getEdgeLength('27')
        L37 = self.getEdgeLength('37')
        L47 = self.getEdgeLength('47')
        L57 = self.getEdgeLength('57')
        
        L28 = self.getEdgeLength('28')
        L58 = self.getEdgeLength('58')
        L68 = self.getEdgeLength('68')
        L78 = self.getEdgeLength('78')
        
        L23 = self.getEdgeLength('23')
        L34 = self.getEdgeLength('34')
        L45 = self.getEdgeLength('45')
        L56 = self.getEdgeLength('56')
        L26 = self.getEdgeLength('26')
        
        y1 = (-L12**2 + L13**2 - L23**2)/float(-2*L23)
        if L12**2  - y1**2 <0:
            raise TriangleInequalityError('Triangle inequality violated!')
        x1 = np.sqrt(L12**2  - y1**2 )
            
        eqs = [
#            L12**2 - x1**2 - y1**2 ,
            -(L23 - y4)**2 + L34**2 - x4**2 - z4**2 ,
            -L12**2 + L15**2 + 2*x1*x5 - x5**2 + 2*y1*y5 - y5**2 - z5**2 ,
            L26**2 - x6**2 - y6**2 - z6**2 ,
            L27**2 - x7**2 - y7**2 - z7**2 ,
            L28**2 - x8**2 - y8**2 - z8**2 ,
            L12**2 - L15**2 + L23**2 - L34**2 + L45**2 - 2*x1*x5 + 2*x4*x5 - 2*L23*y4 - 2*y1*y5 + 2*y4*y5 + 2*z4*z5 ,
            -L12**2 + L16**2 - L26**2 + 2*x1*x6 + 2*y1*y6 ,
            L12**2 - L15**2 - L26**2 + L56**2 - 2*x1*x5 + 2*x5*x6 - 2*y1*y5 + 2*y5*y6 + 2*z5*z6 ,
            L12**2 - L15**2 - L28**2 + L58**2 - 2*x1*x5 + 2*x5*x8 - 2*y1*y5 + 2*y5*y8 + 2*z5*z8 ,
            -L26**2 - L28**2 + L68**2 + 2*x6*x8 + 2*y6*y8 + 2*z6*z8 ,
            -L12**2 + L14**2 + L23**2 - L34**2 + 2*x1*x4 - 2*L23*y4 + 2*y1*y4 ,
#            -L12**2 + L13**2 - L23**2 + 2*L23*y1 ,
            -L23**2 - L27**2 + L37**2 + 2*L23*y7 ,
            L23**2 - L27**2 - L34**2 + L47**2 + 2*x4*x7 - 2*L23*y4 + 2*y4*y7 + 2*z4*z7 ,
            L12**2 - L15**2 - L27**2 + L57**2 - 2*x1*x5 + 2*x5*x7 - 2*y1*y5 + 2*y5*y7 + 2*z5*z7 ,
            -L27**2 - L28**2 + L78**2 + 2*x7*x8 + 2*y7*y8 + 2*z7*z8
        ]
        res = []
        for eq in eqs:
            res.append(str(eq)+';')
        return res

    def constructEquations_max8vertices_distSystem(self):
        '''system with correct mixed volume'''        
        y1, y2, y3, y4 = symbols('y1 y2 y3 y4')
        
        L12 = self.getEdgeLength('12')
        L13 = self.getEdgeLength('13')
        L14 = self.getEdgeLength('14')
        L15 = self.getEdgeLength('15')
        L16 = self.getEdgeLength('16')

        L27 = self.getEdgeLength('27')
        L37 = self.getEdgeLength('37')
        L47 = self.getEdgeLength('47')
        L57 = self.getEdgeLength('57')
        
        L28 = self.getEdgeLength('28')
        L58 = self.getEdgeLength('58')
        L68 = self.getEdgeLength('68')
        L78 = self.getEdgeLength('78')
        
        L23 = self.getEdgeLength('23')
        L34 = self.getEdgeLength('34')
        L45 = self.getEdgeLength('45')
        L56 = self.getEdgeLength('56')
        L26 = self.getEdgeLength('26')
        
       
        eqs = [
        (-L12**2*L34**2+2*L12**2*L34*L37+2*L12**2*L34*L47-L12**2*L37**2+2*L12**2*L37*L47-L12**2*L47**2+4*L12*L13*L23*L47-2*L12*L13*L27*L34
        +2*L12*L13*L27*L37-2*L12*L13*L27*L47-2*L12*L13*L34*L47+2*L12*L13*L34*y3-2*L12*L13*L37*L47-2*L12*L13*L37*y3+2*L12*L13*L47**2
        -2*L12*L13*L47*y3+2*L12*L14*L23*L34-2*L12*L14*L23*L37-2*L12*L14*L23*L47-2*L12*L14*L27*L34-2*L12*L14*L27*L37+2*L12*L14*L27*L47
        -2*L12*L14*L34*L37+2*L12*L14*L37**2-2*L12*L14*L37*L47+4*L12*L14*L37*y3-2*L12*L23*L34*L47-2*L12*L23*L34*y1-2*L12*L23*L37*L47
        +2*L12*L23*L37*y1+2*L12*L23*L47**2-2*L12*L23*L47*y1+2*L12*L27*L34**2-2*L12*L27*L34*L37-2*L12*L27*L34*L47+4*L12*L27*L34*y1
        +2*L12*L34**2*y1+4*L12*L34*L37*L47-2*L12*L34*L37*y1-2*L12*L34*L37*y3-2*L12*L34*L47*y1-2*L12*L34*y1*y3+2*L12*L37**2*y3
        -2*L12*L37*L47*y3-2*L12*L37*y1*y3+2*L12*L47*y1*y3-L13**2*L27**2+2*L13**2*L27*L47+2*L13**2*L27*y3-L13**2*L47**2+2*L13**2*L47*y3
        -L13**2*y3**2-2*L13*L14*L23*L27-2*L13*L14*L23*L47+2*L13*L14*L23*y3+2*L13*L14*L27**2+4*L13*L14*L27*L34-2*L13*L14*L27*L37-2*L13*L14*L27*L47
        -2*L13*L14*L27*y3+2*L13*L14*L37*L47-2*L13*L14*L37*y3-2*L13*L23*L27*L47+2*L13*L23*L27*y1+2*L13*L23*L47**2-2*L13*L23*L47*y1-2*L13*L23*L47*y3
        -2*L13*L23*y1*y3+2*L13*L27**2*L34-2*L13*L27*L34*L47-2*L13*L27*L34*y1-2*L13*L27*L34*y3-2*L13*L27*L37*y3+4*L13*L27*L47*y3-2*L13*L27*y1*y3
        +2*L13*L34*L47*y1-2*L13*L34*y1*y3-2*L13*L37*L47*y3+4*L13*L37*y1*y3+2*L13*L37*y3**2-2*L13*L47*y1*y3+2*L13*y1*y3**2-L14**2*L23**2
        +2*L14**2*L23*L27+2*L14**2*L23*L37-L14**2*L27**2+2*L14**2*L27*L37-L14**2*L37**2+2*L14*L23**2*L47+2*L14*L23**2*y1-2*L14*L23*L27*L34
        +4*L14*L23*L27*L37-2*L14*L23*L27*L47-2*L14*L23*L27*y1-2*L14*L23*L34*y1-2*L14*L23*L37*L47-2*L14*L23*L37*y1-2*L14*L23*L37*y3+4*L14*L23*L47*y1
        -2*L14*L23*y1*y3+2*L14*L27**2*L34-2*L14*L27*L34*L37-2*L14*L27*L34*y1-2*L14*L27*L37*y3+2*L14*L27*y1*y3+2*L14*L34*L37*y1+2*L14*L37**2*y3
        -2*L14*L37*y1*y3-L23**2*L47**2+2*L23**2*L47*y1-L23**2*y1**2+2*L23*L27*L34*L47-2*L23*L27*L34*y1-2*L23*L34*L47*y1+2*L23*L34*y1**2
        +4*L23*L34*y1*y3+2*L23*L37*L47*y3-2*L23*L37*y1*y3-2*L23*L47*y1*y3+2*L23*y1**2*y3-L27**2*L34**2+2*L27*L34**2*y1+2*L27*L34*L37*y3
        -2*L27*L34*y1*y3-L34**2*y1**2-2*L34*L37*y1*y3+2*L34*y1**2*y3-L37**2*y3**2+2*L37*y1*y3**2-y1**2*y3**2),
        (-L12**2*L45**2+2*L12**2*L45*L47+2*L12**2*L45*L57-L12**2*L47**2+2*L12**2*L47*L57-L12**2*L57**2-2*L12*L14*L27*L45+2*L12*L14*L27*L47
        -2*L12*L14*L27*L57-2*L12*L14*L45*L57+2*L12*L14*L45*y4-2*L12*L14*L47*L57-2*L12*L14*L47*y4+2*L12*L14*L57**2+4*L12*L14*L57*y3-2*L12*L14*L57*y4
        -2*L12*L15*L27*L45-2*L12*L15*L27*L47+2*L12*L15*L27*L57-2*L12*L15*L45*L47+2*L12*L15*L45*y3+2*L12*L15*L47**2-2*L12*L15*L47*L57-2*L12*L15*L47*y3
        +4*L12*L15*L47*y4-2*L12*L15*L57*y3+2*L12*L27*L45**2-2*L12*L27*L45*L47-2*L12*L27*L45*L57+4*L12*L27*L45*y1+2*L12*L45**2*y1+4*L12*L45*L47*L57
        -2*L12*L45*L47*y1-2*L12*L45*L47*y4-2*L12*L45*L57*y1-2*L12*L45*L57*y3-2*L12*L45*y1*y3-2*L12*L45*y1*y4+2*L12*L47**2*y4-2*L12*L47*L57*y3
        -2*L12*L47*L57*y4+2*L12*L47*y1*y3-2*L12*L47*y1*y4+2*L12*L57**2*y3-2*L12*L57*y1*y3+2*L12*L57*y1*y4-L14**2*L27**2+2*L14**2*L27*L57
        +2*L14**2*L27*y4-L14**2*L57**2+2*L14**2*L57*y4-L14**2*y4**2+2*L14*L15*L27**2+4*L14*L15*L27*L45-2*L14*L15*L27*L47-2*L14*L15*L27*L57
        -2*L14*L15*L27*y3-2*L14*L15*L27*y4+2*L14*L15*L47*L57-2*L14*L15*L47*y4-2*L14*L15*L57*y3+2*L14*L15*y3*y4+2*L14*L27**2*L45-2*L14*L27*L45*L57
        -2*L14*L27*L45*y1-2*L14*L27*L45*y4-2*L14*L27*L47*y4-2*L14*L27*L57*y3+4*L14*L27*L57*y4+2*L14*L27*y1*y3-2*L14*L27*y1*y4+2*L14*L45*L57*y1
        -2*L14*L45*y1*y4-2*L14*L47*L57*y4+4*L14*L47*y1*y4+2*L14*L47*y4**2+2*L14*L57**2*y3-2*L14*L57*y1*y3-2*L14*L57*y1*y4-2*L14*L57*y3*y4
        -2*L14*y1*y3*y4+2*L14*y1*y4**2-L15**2*L27**2+2*L15**2*L27*L47+2*L15**2*L27*y3-L15**2*L47**2+2*L15**2*L47*y3-L15**2*y3**2+2*L15*L27**2*L45
        -2*L15*L27*L45*L47-2*L15*L27*L45*y1-2*L15*L27*L45*y3+4*L15*L27*L47*y3-2*L15*L27*L47*y4-2*L15*L27*L57*y3-2*L15*L27*y1*y3+2*L15*L27*y1*y4+2*L15*L45*L47*y1
        -2*L15*L45*y1*y3+2*L15*L47**2*y4-2*L15*L47*L57*y3-2*L15*L47*y1*y3-2*L15*L47*y1*y4-2*L15*L47*y3*y4+4*L15*L57*y1*y3+2*L15*L57*y3**2+2*L15*y1*y3**2
        -2*L15*y1*y3*y4-L27**2*L45**2+2*L27*L45**2*y1+2*L27*L45*L47*y4+2*L27*L45*L57*y3-2*L27*L45*y1*y3-2*L27*L45*y1*y4-L45**2*y1**2-2*L45*L47*y1*y4
        -2*L45*L57*y1*y3+2*L45*y1**2*y3+2*L45*y1**2*y4+4*L45*y1*y3*y4-L47**2*y4**2+2*L47*L57*y3*y4-2*L47*y1*y3*y4+2*L47*y1*y4**2-L57**2*y3**2+2*L57*y1*y3**2
        -2*L57*y1*y3*y4-y1**2*y3**2+2*y1**2*y3*y4-y1**2*y4**2),
        (-L12**2*L56**2+2*L12**2*L56*L58+2*L12**2*L56*L68-L12**2*L58**2+2*L12**2*L58*L68-L12**2*L68**2+2*L12*L15*L26*L56-2*L12*L15*L26*L58
        -2*L12*L15*L26*L68-2*L12*L15*L28*L56+2*L12*L15*L28*L58-2*L12*L15*L28*L68-2*L12*L15*L56*L68-2*L12*L15*L58*L68+2*L12*L15*L68**2
        +4*L12*L15*L68*y4+4*L12*L16*L26*L58-2*L12*L16*L28*L56-2*L12*L16*L28*L58+2*L12*L16*L28*L68-2*L12*L16*L56*L58+2*L12*L16*L56*y4
        +2*L12*L16*L58**2-2*L12*L16*L58*L68-2*L12*L16*L58*y4-2*L12*L16*L68*y4-2*L12*L26*L56*L58-2*L12*L26*L56*y2+2*L12*L26*L58**2
        -2*L12*L26*L58*L68-2*L12*L26*L58*y2+2*L12*L26*L68*y2+2*L12*L28*L56**2-2*L12*L28*L56*L58-2*L12*L28*L56*L68+4*L12*L28*L56*y2
        +2*L12*L56**2*y2+4*L12*L56*L58*L68-2*L12*L56*L58*y2-2*L12*L56*L68*y2-2*L12*L56*L68*y4-2*L12*L56*y2*y4-2*L12*L58*L68*y4
        +2*L12*L58*y2*y4+2*L12*L68**2*y4-2*L12*L68*y2*y4-L15**2*L26**2+2*L15**2*L26*L28+2*L15**2*L26*L68-L15**2*L28**2+2*L15**2*L28*L68
        -L15**2*L68**2-2*L15*L16*L26*L28-2*L15*L16*L26*L58+2*L15*L16*L26*y4+2*L15*L16*L28**2+4*L15*L16*L28*L56-2*L15*L16*L28*L58-2*L15*L16*L28*L68
        -2*L15*L16*L28*y4+2*L15*L16*L58*L68-2*L15*L16*L68*y4+2*L15*L26**2*L58+2*L15*L26**2*y2-2*L15*L26*L28*L56-2*L15*L26*L28*L58+4*L15*L26*L28*L68
        -2*L15*L26*L28*y2-2*L15*L26*L56*y2-2*L15*L26*L58*L68+4*L15*L26*L58*y2-2*L15*L26*L68*y2-2*L15*L26*L68*y4-2*L15*L26*y2*y4+2*L15*L28**2*L56
        -2*L15*L28*L56*L68-2*L15*L28*L56*y2-2*L15*L28*L68*y4+2*L15*L28*y2*y4+2*L15*L56*L68*y2+2*L15*L68**2*y4-2*L15*L68*y2*y4-L16**2*L28**2
        +2*L16**2*L28*L58+2*L16**2*L28*y4-L16**2*L58**2+2*L16**2*L58*y4-L16**2*y4**2-2*L16*L26*L28*L58+2*L16*L26*L28*y2+2*L16*L26*L58**2
        -2*L16*L26*L58*y2-2*L16*L26*L58*y4-2*L16*L26*y2*y4+2*L16*L28**2*L56-2*L16*L28*L56*L58-2*L16*L28*L56*y2-2*L16*L28*L56*y4+4*L16*L28*L58*y4
        -2*L16*L28*L68*y4-2*L16*L28*y2*y4+2*L16*L56*L58*y2-2*L16*L56*y2*y4-2*L16*L58*L68*y4-2*L16*L58*y2*y4+4*L16*L68*y2*y4+2*L16*L68*y4**2+2*L16*y2*y4**2
        -L26**2*L58**2+2*L26**2*L58*y2-L26**2*y2**2+2*L26*L28*L56*L58-2*L26*L28*L56*y2-2*L26*L56*L58*y2+2*L26*L56*y2**2+4*L26*L56*y2*y4+2*L26*L58*L68*y4
        -2*L26*L58*y2*y4-2*L26*L68*y2*y4+2*L26*y2**2*y4-L28**2*L56**2+2*L28*L56**2*y2+2*L28*L56*L68*y4-2*L28*L56*y2*y4-L56**2*y2**2-2*L56*L68*y2*y4+2*L56*y2**2*y4
        -L68**2*y4**2+2*L68*y2*y4**2-y2**2*y4**2),
        (-L12**2*L57**2+2*L12**2*L57*L58+2*L12**2*L57*L78-L12**2*L58**2+2*L12**2*L58*L78-L12**2*L78**2+2*L12*L15*L27*L57-2*L12*L15*L27*L58-2*L12*L15*L27*L78
        -2*L12*L15*L28*L57+2*L12*L15*L28*L58-2*L12*L15*L28*L78-2*L12*L15*L57*L78-2*L12*L15*L58*L78+2*L12*L15*L78**2+4*L12*L15*L78*y4-2*L12*L27*L57*L58-2*L12*L27*L57*y2
        +2*L12*L27*L58**2-2*L12*L27*L58*L78+4*L12*L27*L58*y1-2*L12*L27*L58*y2+2*L12*L27*L78*y2+2*L12*L28*L57**2-2*L12*L28*L57*L58-2*L12*L28*L57*L78-2*L12*L28*L57*y1
        +4*L12*L28*L57*y2-2*L12*L28*L58*y1+2*L12*L28*L78*y1+2*L12*L57**2*y2+4*L12*L57*L58*L78-2*L12*L57*L58*y1-2*L12*L57*L58*y2-2*L12*L57*L78*y2-2*L12*L57*L78*y4
        +2*L12*L57*y1*y4-2*L12*L57*y2*y4+2*L12*L58**2*y1-2*L12*L58*L78*y1-2*L12*L58*L78*y4-2*L12*L58*y1*y4+2*L12*L58*y2*y4+2*L12*L78**2*y4-2*L12*L78*y1*y4
        -2*L12*L78*y2*y4-L15**2*L27**2+2*L15**2*L27*L28+2*L15**2*L27*L78-L15**2*L28**2+2*L15**2*L28*L78-L15**2*L78**2+2*L15*L27**2*L58+2*L15*L27**2*y2
        -2*L15*L27*L28*L57-2*L15*L27*L28*L58+4*L15*L27*L28*L78-2*L15*L27*L28*y1-2*L15*L27*L28*y2-2*L15*L27*L57*y2-2*L15*L27*L58*L78-2*L15*L27*L58*y1
        +4*L15*L27*L58*y2-2*L15*L27*L78*y2-2*L15*L27*L78*y4+2*L15*L27*y1*y4-2*L15*L27*y2*y4+2*L15*L28**2*L57+2*L15*L28**2*y1-2*L15*L28*L57*L78+4*L15*L28*L57*y1
        -2*L15*L28*L57*y2-2*L15*L28*L58*y1-2*L15*L28*L78*y1-2*L15*L28*L78*y4-2*L15*L28*y1*y4+2*L15*L28*y2*y4+2*L15*L57*L78*y2+2*L15*L58*L78*y1+2*L15*L78**2*y4
        -2*L15*L78*y1*y4-2*L15*L78*y2*y4-L27**2*L58**2+2*L27**2*L58*y2-L27**2*y2**2+2*L27*L28*L57*L58-2*L27*L28*L57*y2-2*L27*L28*L58*y1+2*L27*L28*y1*y2
        -2*L27*L57*L58*y2+2*L27*L57*y2**2+4*L27*L57*y2*y4+2*L27*L58**2*y1+2*L27*L58*L78*y4-2*L27*L58*y1*y2-2*L27*L58*y1*y4-2*L27*L58*y2*y4-2*L27*L78*y2*y4
        -2*L27*y1*y2*y4+2*L27*y2**2*y4-L28**2*L57**2+2*L28**2*L57*y1-L28**2*y1**2+2*L28*L57**2*y2-2*L28*L57*L58*y1+2*L28*L57*L78*y4-2*L28*L57*y1*y2
        -2*L28*L57*y1*y4-2*L28*L57*y2*y4+2*L28*L58*y1**2+4*L28*L58*y1*y4-2*L28*L78*y1*y4+2*L28*y1**2*y4-2*L28*y1*y2*y4-L57**2*y2**2+2*L57*L58*y1*y2
        -2*L57*L78*y2*y4-2*L57*y1*y2*y4+2*L57*y2**2*y4-L58**2*y1**2-2*L58*L78*y1*y4+2*L58*y1**2*y4-2*L58*y1*y2*y4-L78**2*y4**2+4*L78*y1*y2*y4+2*L78*y1*y4**2
        +2*L78*y2*y4**2-y1**2*y4**2+2*y1*y2*y4**2-y2**2*y4**2)
        ]
        res = []
        for eq in eqs:
            res.append(str(eq)+';')
        return res

    def inequalities_max8vertices(self):
        L12 = self.getEdgeLength('12')
        L13 = self.getEdgeLength('13')
        L14 = self.getEdgeLength('14')
        L15 = self.getEdgeLength('15')
        L16 = self.getEdgeLength('16')

        L27 = self.getEdgeLength('27')
        L37 = self.getEdgeLength('37')
        L47 = self.getEdgeLength('47')
        L57 = self.getEdgeLength('57')

        L28 = self.getEdgeLength('28')
        L58 = self.getEdgeLength('58')
        L68 = self.getEdgeLength('68')
        L78 = self.getEdgeLength('78')

        L23 = self.getEdgeLength('23')
        L34 = self.getEdgeLength('34')
        L45 = self.getEdgeLength('45')
        L56 = self.getEdgeLength('56')
        L26 = self.getEdgeLength('26')
        cond_y1_left=[L12+L27-2*math.sqrt(L12*L27), L13+L37-2*math.sqrt(L13*L37), L14+L47-2*math.sqrt(L14*L47), L15+L57-2*math.sqrt(L15*L57),
                       -(1/4)*(-2*L12*L23+2*L12*L27-2*L12*L37-2*L13*L23-2*L13*L27+2*L13*L37+2*L23**2-2*L23*L27-2*L23*L37
                               +2*math.sqrt(L12**2*L23**2-2*L12**2*L23*L27-2*L12**2*L23*L37+L12**2*L27**2-2*L12**2*L27*L37
                                            +L12**2*L37**2-2*L12*L13*L23**2+4*L12*L13*L23*L27+4*L12*L13*L23*L37-2*L12*L13*L27**2+4*L12*L13*L27*L37
                                            -2*L12*L13*L37**2-2*L12*L23**3+4*L12*L23**2*L27+4*L12*L23**2*L37-2*L12*L23*L27**2+4*L12*L23*L27*L37
                                            -2*L12*L23*L37**2+L13**2*L23**2-2*L13**2*L23*L27-2*L13**2*L23*L37+L13**2*L27**2-2*L13**2*L27*L37+L13**2*L37**2
                                            -2*L13*L23**3+4*L13*L23**2*L27+4*L13*L23**2*L37-2*L13*L23*L27**2+4*L13*L23*L27*L37-2*L13*L23*L37**2+L23**4
                                            -2*L23**3*L27-2*L23**3*L37+L23**2*L27**2-2*L23**2*L27*L37+L23**2*L37**2))/L23,
                        -(1/4)*(-2*L13*L34+2*L13*L37-2*L13*L47-2*L14*L34-2*L14*L37+2*L14*L47+2*L34**2-2*L34*L37-2*L34*L47
                                +2*math.sqrt(L13**2*L34**2-2*L13**2*L34*L37-2*L13**2*L34*L47+L13**2*L37**2-2*L13**2*L37*L47+L13**2*L47**2
                                             -2*L13*L14*L34**2+4*L13*L14*L34*L37+4*L13*L14*L34*L47-2*L13*L14*L37**2+4*L13*L14*L37*L47-2*L13*L14*L47**2
                                             -2*L13*L34**3+4*L13*L34**2*L37+4*L13*L34**2*L47-2*L13*L34*L37**2+4*L13*L34*L37*L47-2*L13*L34*L47**2
                                             +L14**2*L34**2-2*L14**2*L34*L37-2*L14**2*L34*L47+L14**2*L37**2-2*L14**2*L37*L47+L14**2*L47**2
                                             -2*L14*L34**3+4*L14*L34**2*L37+4*L14*L34**2*L47-2*L14*L34*L37**2+4*L14*L34*L37*L47-2*L14*L34*L47**2
                                             +L34**4-2*L34**3*L37-2*L34**3*L47+L34**2*L37**2-2*L34**2*L37*L47+L34**2*L47**2))/L34,
                        -(1/4)*(-2*L14*L45+2*L14*L47-2*L14*L57-2*L15*L45-2*L15*L47+2*L15*L57+2*L45**2-2*L45*L47-2*L45*L57
                                +2*math.sqrt(L14**2*L45**2-2*L14**2*L45*L47-2*L14**2*L45*L57+L14**2*L47**2-2*L14**2*L47*L57+L14**2*L57**2
                                             -2*L14*L15*L45**2+4*L14*L15*L45*L47+4*L14*L15*L45*L57-2*L14*L15*L47**2+4*L14*L15*L47*L57
                                             -2*L14*L15*L57**2-2*L14*L45**3+4*L14*L45**2*L47+4*L14*L45**2*L57-2*L14*L45*L47**2+4*L14*L45*L47*L57
                                             -2*L14*L45*L57**2+L15**2*L45**2-2*L15**2*L45*L47-2*L15**2*L45*L57+L15**2*L47**2-2*L15**2*L47*L57
                                             +L15**2*L57**2-2*L15*L45**3+4*L15*L45**2*L47+4*L15*L45**2*L57-2*L15*L45*L47**2+4*L15*L45*L47*L57
                                             -2*L15*L45*L57**2+L45**4-2*L45**3*L47-2*L45**3*L57+L45**2*L47**2-2*L45**2*L47*L57+L45**2*L57**2))/L45]

        cond_y1_right=[L12+L27+2*math.sqrt(L12*L27), L13+L37+2*math.sqrt(L13*L37), L14+L47+2*math.sqrt(L14*L47), L15+L57+2*math.sqrt(L15*L57), 
                       -(1/4)*(-2*L12*L23+2*L12*L27-2*L12*L37-2*L13*L23-2*L13*L27+2*L13*L37+2*L23**2-2*L23*L27-2*L23*L37
                               -2*math.sqrt(L12**2*L23**2-2*L12**2*L23*L27-2*L12**2*L23*L37+L12**2*L27**2-2*L12**2*L27*L37+L12**2*L37**2-2*L12*L13*L23**2
                                            +4*L12*L13*L23*L27+4*L12*L13*L23*L37-2*L12*L13*L27**2+4*L12*L13*L27*L37-2*L12*L13*L37**2-2*L12*L23**3
                                            +4*L12*L23**2*L27+4*L12*L23**2*L37-2*L12*L23*L27**2+4*L12*L23*L27*L37-2*L12*L23*L37**2+L13**2*L23**2
                                            -2*L13**2*L23*L27-2*L13**2*L23*L37+L13**2*L27**2-2*L13**2*L27*L37+L13**2*L37**2-2*L13*L23**3+4*L13*L23**2*L27
                                            +4*L13*L23**2*L37-2*L13*L23*L27**2+4*L13*L23*L27*L37-2*L13*L23*L37**2+L23**4-2*L23**3*L27-2*L23**3*L37
                                            +L23**2*L27**2-2*L23**2*L27*L37+L23**2*L37**2))/L23,
                        -(1/4)*(-2*L13*L34+2*L13*L37-2*L13*L47-2*L14*L34-2*L14*L37+2*L14*L47+2*L34**2-2*L34*L37-2*L34*L47
                                -2*math.sqrt(L13**2*L34**2-2*L13**2*L34*L37-2*L13**2*L34*L47+L13**2*L37**2-2*L13**2*L37*L47+L13**2*L47**2
                                             -2*L13*L14*L34**2+4*L13*L14*L34*L37+4*L13*L14*L34*L47-2*L13*L14*L37**2+4*L13*L14*L37*L47-2*L13*L14*L47**2
                                             -2*L13*L34**3+4*L13*L34**2*L37+4*L13*L34**2*L47-2*L13*L34*L37**2+4*L13*L34*L37*L47-2*L13*L34*L47**2
                                             +L14**2*L34**2-2*L14**2*L34*L37-2*L14**2*L34*L47+L14**2*L37**2-2*L14**2*L37*L47+L14**2*L47**2-2*L14*L34**3
                                             +4*L14*L34**2*L37+4*L14*L34**2*L47-2*L14*L34*L37**2+4*L14*L34*L37*L47-2*L14*L34*L47**2+L34**4-2*L34**3*L37
                                             -2*L34**3*L47+L34**2*L37**2-2*L34**2*L37*L47+L34**2*L47**2))/L34,
                        -(1/4)*(-2*L14*L45+2*L14*L47-2*L14*L57-2*L15*L45-2*L15*L47+2*L15*L57+2*L45**2-2*L45*L47-2*L45*L57
                                -2*math.sqrt(L14**2*L45**2-2*L14**2*L45*L47-2*L14**2*L45*L57+L14**2*L47**2-2*L14**2*L47*L57+L14**2*L57**2
                                             -2*L14*L15*L45**2+4*L14*L15*L45*L47+4*L14*L15*L45*L57-2*L14*L15*L47**2+4*L14*L15*L47*L57-2*L14*L15*L57**2
                                             -2*L14*L45**3+4*L14*L45**2*L47+4*L14*L45**2*L57-2*L14*L45*L47**2+4*L14*L45*L47*L57-2*L14*L45*L57**2
                                             +L15**2*L45**2-2*L15**2*L45*L47-2*L15**2*L45*L57+L15**2*L47**2-2*L15**2*L47*L57+L15**2*L57**2
                                             -2*L15*L45**3+4*L15*L45**2*L47+4*L15*L45**2*L57-2*L15*L45*L47**2+4*L15*L45*L47*L57-2*L15*L45*L57**2+L45**4
                                             -2*L45**3*L47-2*L45**3*L57+L45**2*L47**2-2*L45**2*L47*L57+L45**2*L57**2))/L45]

        Int_y1=[max(cond_y1_left),min(cond_y1_right)]

        ##y1 is embeddable iff Int_y1[0]<=y1<=Int_y1[1]

        cond_y4_left=[L12+L15-2*math.sqrt(L12*L15), L26+L56-2*math.sqrt(L26*L56), L27+L57-2*math.sqrt(L27*L57),
                       -(1/4)*(2*L12*L15-2*L12*L16-2*L12*L56-2*L15*L16-2*L15*L26+2*L16**2-2*L16*L26-2*L16*L56+2*L26*L56
                               +2*math.sqrt(L12**2*L15**2-2*L12**2*L15*L16-2*L12**2*L15*L56+L12**2*L16**2-2*L12**2*L16*L56+L12**2*L56**2-2*L12*L15**2*L16
                                            -2*L12*L15**2*L26+4*L12*L15*L16**2+4*L12*L15*L16*L26+4*L12*L15*L16*L56+4*L12*L15*L26*L56-2*L12*L16**3
                                            -2*L12*L16**2*L26+4*L12*L16**2*L56+4*L12*L16*L26*L56-2*L12*L16*L56**2-2*L12*L26*L56**2+L15**2*L16**2
                                            -2*L15**2*L16*L26+L15**2*L26**2-2*L15*L16**3+4*L15*L16**2*L26-2*L15*L16**2*L56-2*L15*L16*L26**2
                                            +4*L15*L16*L26*L56-2*L15*L26**2*L56+L16**4-2*L16**3*L26-2*L16**3*L56+L16**2*L26**2+4*L16**2*L26*L56+L16**2*L56**2-2*L16*L26**2*L56-2*L16*L26*L56**2+L26**2*L56**2))/L16, -(1/4)*(2*L26*L56-2*L26*L58-2*L26*L68-2*L28*L56+2*L28*L58-2*L28*L68-2*L56*L68-2*L58*L68+2*L68**2+2*math.sqrt(L26**2*L56**2-2*L26**2*L56*L58-2*L26**2*L56*L68+L26**2*L58**2-2*L26**2*L58*L68+L26**2*L68**2-2*L26*L28*L56**2+4*L26*L28*L56*L58+4*L26*L28*L56*L68-2*L26*L28*L58**2+4*L26*L28*L58*L68-2*L26*L28*L68**2-2*L26*L56**2*L68+4*L26*L56*L58*L68+4*L26*L56*L68**2-2*L26*L58**2*L68+4*L26*L58*L68**2-2*L26*L68**3+L28**2*L56**2-2*L28**2*L56*L58-2*L28**2*L56*L68+L28**2*L58**2-2*L28**2*L58*L68+L28**2*L68**2-2*L28*L56**2*L68+4*L28*L56*L58*L68+4*L28*L56*L68**2-2*L28*L58**2*L68+4*L28*L58*L68**2-2*L28*L68**3+L56**2*L68**2-2*L56*L58*L68**2-2*L56*L68**3+L58**2*L68**2-2*L58*L68**3+L68**4))/L68, -(1/4)*(2*L27*L57-2*L27*L58-2*L27*L78-2*L28*L57+2*L28*L58-2*L28*L78-2*L57*L78-2*L58*L78+2*L78**2+2*math.sqrt(L27**2*L57**2-2*L27**2*L57*L58-2*L27**2*L57*L78+L27**2*L58**2-2*L27**2*L58*L78+L27**2*L78**2-2*L27*L28*L57**2+4*L27*L28*L57*L58+4*L27*L28*L57*L78-2*L27*L28*L58**2+4*L27*L28*L58*L78-2*L27*L28*L78**2-2*L27*L57**2*L78+4*L27*L57*L58*L78+4*L27*L57*L78**2-2*L27*L58**2*L78+4*L27*L58*L78**2-2*L27*L78**3+L28**2*L57**2-2*L28**2*L57*L58-2*L28**2*L57*L78+L28**2*L58**2-2*L28**2*L58*L78+L28**2*L78**2-2*L28*L57**2*L78+4*L28*L57*L58*L78+4*L28*L57*L78**2-2*L28*L58**2*L78+4*L28*L58*L78**2-2*L28*L78**3+L57**2*L78**2-2*L57*L58*L78**2-2*L57*L78**3+L58**2*L78**2-2*L58*L78**3+L78**4))/L78]

        cond_y4_right=[L12+L15+2*math.sqrt(L12*L15), L26+L56+2*math.sqrt(L26*L56), L27+L57+2*math.sqrt(L27*L57), -(1/4)*(2*L12*L15-2*L12*L16-2*L12*L56-2*L15*L16-2*L15*L26+2*L16**2-2*L16*L26-2*L16*L56+2*L26*L56-2*math.sqrt(L12**2*L15**2-2*L12**2*L15*L16-2*L12**2*L15*L56+L12**2*L16**2-2*L12**2*L16*L56+L12**2*L56**2-2*L12*L15**2*L16-2*L12*L15**2*L26+4*L12*L15*L16**2+4*L12*L15*L16*L26+4*L12*L15*L16*L56+4*L12*L15*L26*L56-2*L12*L16**3-2*L12*L16**2*L26+4*L12*L16**2*L56+4*L12*L16*L26*L56-2*L12*L16*L56**2-2*L12*L26*L56**2+L15**2*L16**2-2*L15**2*L16*L26+L15**2*L26**2-2*L15*L16**3+4*L15*L16**2*L26-2*L15*L16**2*L56-2*L15*L16*L26**2+4*L15*L16*L26*L56-2*L15*L26**2*L56+L16**4-2*L16**3*L26-2*L16**3*L56+L16**2*L26**2+4*L16**2*L26*L56+L16**2*L56**2-2*L16*L26**2*L56-2*L16*L26*L56**2+L26**2*L56**2))/L16, -(1/4)*(2*L26*L56-2*L26*L58-2*L26*L68-2*L28*L56+2*L28*L58-2*L28*L68-2*L56*L68-2*L58*L68+2*L68**2-2*math.sqrt(L26**2*L56**2-2*L26**2*L56*L58-2*L26**2*L56*L68+L26**2*L58**2-2*L26**2*L58*L68+L26**2*L68**2-2*L26*L28*L56**2+4*L26*L28*L56*L58+4*L26*L28*L56*L68-2*L26*L28*L58**2+4*L26*L28*L58*L68-2*L26*L28*L68**2-2*L26*L56**2*L68+4*L26*L56*L58*L68+4*L26*L56*L68**2-2*L26*L58**2*L68+4*L26*L58*L68**2-2*L26*L68**3+L28**2*L56**2-2*L28**2*L56*L58-2*L28**2*L56*L68+L28**2*L58**2-2*L28**2*L58*L68+L28**2*L68**2-2*L28*L56**2*L68+4*L28*L56*L58*L68+4*L28*L56*L68**2-2*L28*L58**2*L68+4*L28*L58*L68**2-2*L28*L68**3+L56**2*L68**2-2*L56*L58*L68**2-2*L56*L68**3+L58**2*L68**2-2*L58*L68**3+L68**4))/L68, -(1/4)*(2*L27*L57-2*L27*L58-2*L27*L78-2*L28*L57+2*L28*L58-2*L28*L78-2*L57*L78-2*L58*L78+2*L78**2-2*math.sqrt(L27**2*L57**2-2*L27**2*L57*L58-2*L27**2*L57*L78+L27**2*L58**2-2*L27**2*L58*L78+L27**2*L78**2-2*L27*L28*L57**2+4*L27*L28*L57*L58+4*L27*L28*L57*L78-2*L27*L28*L58**2+4*L27*L28*L58*L78-2*L27*L28*L78**2-2*L27*L57**2*L78+4*L27*L57*L58*L78+4*L27*L57*L78**2-2*L27*L58**2*L78+4*L27*L58*L78**2-2*L27*L78**3+L28**2*L57**2-2*L28**2*L57*L58-2*L28**2*L57*L78+L28**2*L58**2-2*L28**2*L58*L78+L28**2*L78**2-2*L28*L57**2*L78+4*L28*L57*L58*L78+4*L28*L57*L78**2-2*L28*L58**2*L78+4*L28*L58*L78**2-2*L28*L78**3+L57**2*L78**2-2*L57*L58*L78**2-2*L57*L78**3+L58**2*L78**2-2*L58*L78**3+L78**4))/L78]

        Int_y4=[max(cond_y4_left),min(cond_y4_right)]

        ##y1 is embeddable iff Int_y4[0]<=y4<=Int_y4[1]
        ## the mechanism is embeddable if both y1/y4 are so
        ## should we use small tolerance (eg <10^(-5))?

    def constructEquations_ring8vertices(self):
        '''system with correct mixed volume'''    
        x4, y4, z4 = symbols('x4 y4 z4')
        x5, y5, z5 = symbols('x5 y5 z5')
        x6, y6, z6 = symbols('x6 y6 z6')
        x7, y7, z7 = symbols('x7 y7 z7')
        x8, y8, z8 = symbols('x8 y8 z8')
        
        L12 = self.getEdgeLength(1, 2)
        L13 = self.getEdgeLength(1, 3)
        L14 = self.getEdgeLength(1, 4)
        L15 = self.getEdgeLength(1, 5)
        L16 = self.getEdgeLength(1, 6)
        L17 = self.getEdgeLength(1, 7)
        L23 = self.getEdgeLength(2, 3)
        L27 = self.getEdgeLength(2, 7)
        L28 = self.getEdgeLength(2, 8)
        L34 = self.getEdgeLength(3, 4)
        L38 = self.getEdgeLength(3, 8)
        L45 = self.getEdgeLength(4, 5)
        L48 = self.getEdgeLength(4, 8)
        L56 = self.getEdgeLength(5, 6)
        L58 = self.getEdgeLength(5, 8)
        L68 = self.getEdgeLength(6, 8)
        L78 = self.getEdgeLength(7, 8)
        L67 = self.getEdgeLength(6, 7)
        
        y1 = (-L12**2 + L13**2 - L23**2)/float(-2*L23)
        if L12**2  - y1**2 <0:
            raise TriangleInequalityError('Triangle inequality violated!')
        x1 = np.sqrt(L12**2  - y1**2 )
        
        eqs = [
                -L12**2 + L14**2 + 2*x1*x4 - x4**2 + 2*y1*y4 - y4**2 - z4**2 ,
                -L12**2 + L15**2 + 2*x1*x5 - x5**2 + 2*y1*y5 - y5**2 - z5**2 ,
                -L12**2 + L16**2 + 2*x1*x6 - x6**2 + 2*y1*y6 - y6**2 - z6**2 ,
                L27**2 - x7**2 - y7**2 - z7**2 ,
                L28**2 - x8**2 - y8**2 - z8**2 ,
#                -L12**2 + L13**2 - L23**2 + 2*L23*y1 ,
                -L12**2 + L17**2 - L27**2 + 2*x1*x7 + 2*y1*y7 ,
                L12**2 - L14**2 - L23**2 + L34**2 - 2*x1*x4 + 2*L23*y4 - 2*y1*y4 ,
                -L23**2 - L28**2 + L38**2 + 2*L23*y8 ,
                2*L12**2 - L14**2 - L15**2 + L45**2 - 2*x1*x4 - 2*x1*x5 + 2*x4*x5 - 2*y1*y4 - 2*y1*y5 + 2*y4*y5 + 2*z4*z5 ,
                L12**2 - L14**2 - L28**2 + L48**2 - 2*x1*x4 + 2*x4*x8 - 2*y1*y4 + 2*y4*y8 + 2*z4*z8 ,
                2*L12**2 - L15**2 - L16**2 + L56**2 - 2*x1*x5 - 2*x1*x6 + 2*x5*x6 - 2*y1*y5 - 2*y1*y6 + 2*y5*y6 + 2*z5*z6 ,
                L12**2 - L15**2 - L28**2 + L58**2 - 2*x1*x5 + 2*x5*x8 - 2*y1*y5 + 2*y5*y8 + 2*z5*z8 ,
                L12**2 - L16**2 - L28**2 + L68**2 - 2*x1*x6 + 2*x6*x8 - 2*y1*y6 + 2*y6*y8 + 2*z6*z8 ,
                -L27**2 - L28**2 + L78**2 + 2*x7*x8 + 2*y7*y8 + 2*z7*z8 ,
                L12**2 - L16**2 - L27**2 + L67**2 - 2*x1*x6 + 2*x6*x7 - 2*y1*y6 + 2*y6*y7 + 2*z6*z7 
        ]
        res = []
        for eq in eqs:
            res.append(str(eq)+';')
        return res



    def constructEquations_max6vertices(self):
        '''system with correct mixed volume'''        
        x4, y4, z4 = symbols('x4 y4 z4')
        x5, y5, z5 = symbols('x5 y5 z5')
        x6, y6, z6 = symbols('x6 y6 z6')
        
#        L12 = self.getEdgeLength(1, 2)
#        L23 = self.getEdgeLength(2, 3)
        L34 = self.getEdgeLength(3, 4)
        L45 = self.getEdgeLength(4, 5)
        L56 = self.getEdgeLength(5, 6)
        L16 = self.getEdgeLength(1, 6)
#        L13 = self.getEdgeLength(1, 3)
#        L24 = self.getEdgeLength(2, 4)
        L35 = self.getEdgeLength(3, 5)
        L46 = self.getEdgeLength(4, 6)
        L15 = self.getEdgeLength(1, 5)
        L26 = self.getEdgeLength(2, 6)
#        
        try:
            yshift = -self.getAltitudeAndFoot(3, 2, 4)[1]
            v2, v3, v1 = self.coordinatesOfTriangle(2, 3, 1, yshift)
        except:
            raise TriangleInequalityError('Triangle inequality violated!')
        X1, Y1, _ = v1
        _, Y2, _ = v2
        _, Y3, _ = v3
        
        eqs = [
            L34**2 - Y3**2 - x4**2 - z4**2 ,
            L35**2 - (Y3 - y5)**2 - x5**2 - z5**2 ,
            L26**2 - (Y2 - y6)**2 - x6**2 - z6**2 ,
            L16**2 - L26**2 - X1**2 - Y1**2 + Y2**2 + 2*X1*x6 + 2*Y1*y6 - 2*Y2*y6 ,
            L15**2 - L35**2 - X1**2 - Y1**2 + Y3**2 + 2*X1*x5 + 2*Y1*y5 - 2*Y3*y5 ,
            -L34**2 - L35**2 + L45**2 + 2*Y3**2 + 2*x4*x5 - 2*Y3*y5 + 2*z4*z5 ,
            -L26**2 - L35**2 + L56**2 + Y2**2 + Y3**2 + 2*x5*x6 - 2*Y3*y5 - 2*Y2*y6 + 2*y5*y6 + 2*z5*z6 ,
            -L26**2 - L34**2 + L46**2 + Y2**2 + Y3**2 + 2*x4*x6 - 2*Y2*y6 + 2*z4*z6 ,
        ]
        res = []
        for eq in eqs:
            res.append(str(eq)+';')
        return res

    def constructEquations_7vert32a(self):
        x1 , y1 , z1 , s1 = symbols('x1 y1 z1 s1')
        x2 , y2 , z2 , s2 = symbols('x2 y2 z2 s2')
        x3 , y3 , z3 , s3 = symbols('x3 y3 z3 s3')
        x4 , y4 , z4 , s4 = symbols('x4 y4 z4 s4')
        x5 , y5 , z5 , s5 = symbols('x5 y5 z5 s5')
        x6 , y6 , z6 , s6 = symbols('x6 y6 z6 s6')
        x7 , y7 , z7 , s7 = symbols('x7 y7 z7 s7')

        L12 = self.getEdgeLength(1, 2)
        L13 = self.getEdgeLength(1, 3)
        L14 = self.getEdgeLength(1, 4)
        L16 = self.getEdgeLength(1, 6)
        L23 = self.getEdgeLength(2, 3)
        L24 = self.getEdgeLength(2, 4)
        L25 = self.getEdgeLength(2, 5)
        L34 = self.getEdgeLength(3, 4)
        L35 = self.getEdgeLength(3, 5)
        L36 = self.getEdgeLength(3, 6)
        L37 = self.getEdgeLength(3, 7)
        L47 = self.getEdgeLength(4, 7)
        L56 = self.getEdgeLength(5, 6)
        L57 = self.getEdgeLength(5, 7)
        L67 = self.getEdgeLength(6, 7)

        y1 = (-L12**2 + L13**2 - L23**2)/float(-2*L23)
        if L12**2  - y1**2 <0:
            raise TriangleInequalityError('Triangle inequality violated!')
        x1 = np.sqrt(L12**2  - y1**2 )

        eqs = [
#        -L12**2 + L13**2 - L23**2 + 2*L23*y1 ,
        -L12**2 + L14**2 + 2*x1*x4 + 2*y1*y4 - s4 ,
        -L12**2 + L16**2 + 2*x1*x6 + 2*y1*y6 - s6 ,
        L24**2 - s4 ,
        L25**2 - s5 ,
        -L23**2 + L34**2 + 2*L23*y4 - s4 ,
        -L23**2 + L35**2 + 2*L23*y5 - s5 ,
        -L23**2 + L36**2 + 2*L23*y6 - s6 ,
        -L23**2 + L37**2 + 2*L23*y7 - s7 ,
        L47**2 + 2*x4*x7 + 2*y4*y7 + 2*z4*z7 - s4 - s7 ,
        L56**2 + 2*x5*x6 + 2*y5*y6 + 2*z5*z6 - s5 - s6 ,
        L57**2 + 2*x5*x7 + 2*y5*y7 + 2*z5*z7 - s5 - s7 ,
        L67**2 + 2*x6*x7 + 2*y6*y7 + 2*z6*z7 - s6 - s7 ,
#        L12**2 - x1**2 - y1**2 ,
        -x4**2 - y4**2 - z4**2 + s4 ,
        -x5**2 - y5**2 - z5**2 + s5 ,
        -x6**2 - y6**2 - z6**2 + s6 ,
        -x7**2 - y7**2 - z7**2 + s7 ,
        ]
        
        res = []
        for eq in eqs:
            res.append(str(eq)+';')
        return res
    
    def constructEquations_7vert32b(self):
        x1 , y1 , z1 , s1 = symbols('x1 y1 z1 s1')
        x2 , y2 , z2 , s2 = symbols('x2 y2 z2 s2')
        x3 , y3 , z3 , s3 = symbols('x3 y3 z3 s3')
        x4 , y4 , z4 , s4 = symbols('x4 y4 z4 s4')
        x5 , y5 , z5 , s5 = symbols('x5 y5 z5 s5')
        x6 , y6 , z6 , s6 = symbols('x6 y6 z6 s6')
        x7 , y7 , z7 , s7 = symbols('x7 y7 z7 s7')

        L12 = self.getEdgeLength(1, 2)
        L13 = self.getEdgeLength(1, 3)
        L14 = self.getEdgeLength(1, 4)
        L15 = self.getEdgeLength(1, 5)
        L23 = self.getEdgeLength(2, 3)
        L25 = self.getEdgeLength(2, 5)
        L26 = self.getEdgeLength(2, 6)
        L27 = self.getEdgeLength(2, 7)
        L34 = self.getEdgeLength(3, 4)
        L36 = self.getEdgeLength(3, 6)
        L37 = self.getEdgeLength(3, 7)
        L45 = self.getEdgeLength(4, 5)
        L47 = self.getEdgeLength(4, 7)
        L56 = self.getEdgeLength(5, 6)
        L67 = self.getEdgeLength(6, 7)

        y1 = (-L12**2 + L13**2 - L23**2)/float(-2*L23)
        if L12**2  - y1**2 <0:
            raise TriangleInequalityError('Triangle inequality violated!')
        x1 = np.sqrt(L12**2  - y1**2 )

        eqs = [
#        -L12**2 + L13**2 - L23**2 + 2*L23*y1 ,
        -L12**2 + L14**2 + 2*x1*x4 + 2*y1*y4 - s4 ,
        -L12**2 + L15**2 + 2*x1*x5 + 2*y1*y5 - s5 ,
        L25**2 - s5 ,
        L26**2 - s6 ,
        L27**2 - s7 ,
        -L23**2 + L34**2 + 2*L23*y4 - s4 ,
        -L23**2 + L36**2 + 2*L23*y6 - s6 ,
        -L23**2 + L37**2 + 2*L23*y7 - s7 ,
        L45**2 + 2*x4*x5 + 2*y4*y5 + 2*z4*z5 - s4 - s5 ,
        L47**2 + 2*x4*x7 + 2*y4*y7 + 2*z4*z7 - s4 - s7 ,
        L56**2 + 2*x5*x6 + 2*y5*y6 + 2*z5*z6 - s5 - s6 ,
        L67**2 + 2*x6*x7 + 2*y6*y7 + 2*z6*z7 - s6 - s7 ,
#        L12**2 - x1**2 - y1**2 ,
        -x4**2 - y4**2 - z4**2 + s4 ,
        -x5**2 - y5**2 - z5**2 + s5 ,
        -x6**2 - y6**2 - z6**2 + s6 ,
        -x7**2 - y7**2 - z7**2 + s7 ,
        ]

        res = []
        for eq in eqs:
            res.append(str(eq)+';')
        return res

    def constructEquations_7vert24(self):
        x1 , y1 , z1 , s1 = symbols('x1 y1 z1 s1')
        x2 , y2 , z2 , s2 = symbols('x2 y2 z2 s2')
        x3 , y3 , z3 , s3 = symbols('x3 y3 z3 s3')
        x4 , y4 , z4 , s4 = symbols('x4 y4 z4 s4')
        x5 , y5 , z5 , s5 = symbols('x5 y5 z5 s5')
        x6 , y6 , z6 , s6 = symbols('x6 y6 z6 s6')
        x7 , y7 , z7 , s7 = symbols('x7 y7 z7 s7')

        L12 = self.getEdgeLength(1, 2)
        L13 = self.getEdgeLength(1, 3)
        L14 = self.getEdgeLength(1, 4)
        L15 = self.getEdgeLength(1, 5)
        L23 = self.getEdgeLength(2, 3)
        L25 = self.getEdgeLength(2, 5)
        L26 = self.getEdgeLength(2, 6)
        L27 = self.getEdgeLength(2, 7)
        L34 = self.getEdgeLength(3, 4)
        L36 = self.getEdgeLength(3, 6)
        L37 = self.getEdgeLength(3, 7)
        L46 = self.getEdgeLength(4, 6)
        L47 = self.getEdgeLength(4, 7)
        L56 = self.getEdgeLength(5, 6)
        L57 = self.getEdgeLength(5, 7)

        y1 = (-L12**2 + L13**2 - L23**2)/float(-2*L23)
        if L12**2  - y1**2 <0:
            raise TriangleInequalityError('Triangle inequality violated!')
        x1 = np.sqrt(L12**2  - y1**2 )

        eqs = [
#        -L12**2 + L13**2 - L23**2 + 2*L23*y1 ,
        -L12**2 + L14**2 + 2*x1*x4 + 2*y1*y4 - s4 ,
        -L12**2 + L15**2 + 2*x1*x5 + 2*y1*y5 - s5 ,
        L25**2 - s5 ,
        L26**2 - s6 ,
        L27**2 - s7 ,
        -L23**2 + L34**2 + 2*L23*y4 - s4 ,
        -L23**2 + L36**2 + 2*L23*y6 - s6 ,
        -L23**2 + L37**2 + 2*L23*y7 - s7 ,
        L46**2 + 2*x4*x6 + 2*y4*y6 + 2*z4*z6 - s4 - s6 ,
        L47**2 + 2*x4*x7 + 2*y4*y7 + 2*z4*z7 - s4 - s7 ,
        L56**2 + 2*x5*x6 + 2*y5*y6 + 2*z5*z6 - s5 - s6 ,
        L57**2 + 2*x5*x7 + 2*y5*y7 + 2*z5*z7 - s5 - s7 ,
#        L12**2 - x1**2 - y1**2 ,
        -x4**2 - y4**2 - z4**2 + s4 ,
        -x5**2 - y5**2 - z5**2 + s5 ,
        -x6**2 - y6**2 - z6**2 + s6 ,
        -x7**2 - y7**2 - z7**2 + s7 ,
        ]

        res = []
        for eq in eqs:
            res.append(str(eq)+';')
        return res

    def constructEquations_7vert16a(self):
        x1 , y1 , z1 , s1 = symbols('x1 y1 z1 s1')
        x2 , y2 , z2 , s2 = symbols('x2 y2 z2 s2')
        x3 , y3 , z3 , s3 = symbols('x3 y3 z3 s3')
        x4 , y4 , z4 , s4 = symbols('x4 y4 z4 s4')
        x5 , y5 , z5 , s5 = symbols('x5 y5 z5 s5')
        x6 , y6 , z6 , s6 = symbols('x6 y6 z6 s6')
        x7 , y7 , z7 , s7 = symbols('x7 y7 z7 s7')

        L12 = self.getEdgeLength(1, 2)
        L13 = self.getEdgeLength(1, 3)
        L16 = self.getEdgeLength(1, 6)
        L17 = self.getEdgeLength(1, 7)
        L23 = self.getEdgeLength(2, 3)
        L24 = self.getEdgeLength(2, 4)
        L25 = self.getEdgeLength(2, 5)
        L34 = self.getEdgeLength(3, 4)
        L35 = self.getEdgeLength(3, 5)
        L36 = self.getEdgeLength(3, 6)
        L37 = self.getEdgeLength(3, 7)
        L46 = self.getEdgeLength(4, 6)
        L47 = self.getEdgeLength(4, 7)
        L56 = self.getEdgeLength(5, 6)
        L57 = self.getEdgeLength(5, 7)

        y1 = (-L12**2 + L13**2 - L23**2)/float(-2*L23)
        if L12**2  - y1**2 <0:
            raise TriangleInequalityError('Triangle inequality violated!')
        x1 = np.sqrt(L12**2  - y1**2 )

        eqs = [
#        -L12**2 + L13**2 - L23**2 + 2*L23*y1 ,
        -L12**2 + L16**2 + 2*x1*x6 + 2*y1*y6 - s6 ,
        -L12**2 + L17**2 + 2*x1*x7 + 2*y1*y7 - s7 ,
        L24**2 - s4 ,
        L25**2 - s5 ,
        -L23**2 + L34**2 + 2*L23*y4 - s4 ,
        -L23**2 + L35**2 + 2*L23*y5 - s5 ,
        -L23**2 + L36**2 + 2*L23*y6 - s6 ,
        -L23**2 + L37**2 + 2*L23*y7 - s7 ,
        L46**2 + 2*x4*x6 + 2*y4*y6 + 2*z4*z6 - s4 - s6 ,
        L47**2 + 2*x4*x7 + 2*y4*y7 + 2*z4*z7 - s4 - s7 ,
        L56**2 + 2*x5*x6 + 2*y5*y6 + 2*z5*z6 - s5 - s6 ,
        L57**2 + 2*x5*x7 + 2*y5*y7 + 2*z5*z7 - s5 - s7 ,
#        L12**2 - x1**2 - y1**2 ,
        -x4**2 - y4**2 - z4**2 + s4 ,
        -x5**2 - y5**2 - z5**2 + s5 ,
        -x6**2 - y6**2 - z6**2 + s6 ,
        -x7**2 - y7**2 - z7**2 + s7 ,
        ]

        res = []
        for eq in eqs:
            res.append(str(eq)+';')
        return res

    def constructEquations_7vert16b(self):
        x1 , y1 , z1 , s1 = symbols('x1 y1 z1 s1')
        x2 , y2 , z2 , s2 = symbols('x2 y2 z2 s2')
        x3 , y3 , z3 , s3 = symbols('x3 y3 z3 s3')
        x4 , y4 , z4 , s4 = symbols('x4 y4 z4 s4')
        x5 , y5 , z5 , s5 = symbols('x5 y5 z5 s5')
        x6 , y6 , z6 , s6 = symbols('x6 y6 z6 s6')
        x7 , y7 , z7 , s7 = symbols('x7 y7 z7 s7')

        L12 = self.getEdgeLength(1, 2)
        L13 = self.getEdgeLength(1, 3)
        L14 = self.getEdgeLength(1, 4)
        L15 = self.getEdgeLength(1, 5)
        L23 = self.getEdgeLength(2, 3)
        L25 = self.getEdgeLength(2, 5)
        L26 = self.getEdgeLength(2, 6)
        L27 = self.getEdgeLength(2, 7)
        L35 = self.getEdgeLength(3, 5)
        L36 = self.getEdgeLength(3, 6)
        L37 = self.getEdgeLength(3, 7)
        L45 = self.getEdgeLength(4, 5)
        L46 = self.getEdgeLength(4, 6)
        L47 = self.getEdgeLength(4, 7)
        L67 = self.getEdgeLength(6, 7)

        y1 = (-L12**2 + L13**2 - L23**2)/float(-2*L23)
        if L12**2  - y1**2 <0:
            raise TriangleInequalityError('Triangle inequality violated!')
        x1 = np.sqrt(L12**2  - y1**2 )

        eqs = [
#        -L12**2 + L13**2 - L23**2 + 2*L23*y1 ,
        -L12**2 + L14**2 + 2*x1*x4 + 2*y1*y4 - s4 ,
        -L12**2 + L15**2 + 2*x1*x5 + 2*y1*y5 - s5 ,
        L25**2 - s5 ,
        L26**2 - s6 ,
        L27**2 - s7 ,
        -L23**2 + L35**2 + 2*L23*y5 - s5 ,
        -L23**2 + L36**2 + 2*L23*y6 - s6 ,
        -L23**2 + L37**2 + 2*L23*y7 - s7 ,
        L45**2 + 2*x4*x5 + 2*y4*y5 + 2*z4*z5 - s4 - s5 ,
        L46**2 + 2*x4*x6 + 2*y4*y6 + 2*z4*z6 - s4 - s6 ,
        L47**2 + 2*x4*x7 + 2*y4*y7 + 2*z4*z7 - s4 - s7 ,
        L67**2 + 2*x6*x7 + 2*y6*y7 + 2*z6*z7 - s6 - s7 ,
#        L12**2 - x1**2 - y1**2 ,
        -x4**2 - y4**2 - z4**2 + s4 ,
        -x5**2 - y5**2 - z5**2 + s5 ,
        -x6**2 - y6**2 - z6**2 + s6 ,
        -x7**2 - y7**2 - z7**2 + s7 ,
        ]

        res = []
        for eq in eqs:
            res.append(str(eq)+';')
        return res

class TriangleInequalityError(ValueError):
    def __init__(self, errorMsg):
        super(TriangleInequalityError, self).__init__(errorMsg)
