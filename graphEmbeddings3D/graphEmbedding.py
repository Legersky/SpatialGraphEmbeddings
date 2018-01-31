#Copyright (C) 2018 Jan Legersky
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

from phcpy.solutions import strsol2dict, is_real
    
from sympy import symbols


import hashlib
import time
import os
import ast
import subprocess

import math

class GraphEmbedding(object):
    '''
    This class implements the computation of spatial embeddings for a graph :math:`G` with edge lengths :math:`\\mathbb{d}`

    Supported graphs:
    
    - :math:`G_{16}`
        - `graph_type='Max6vertices'`
        - edges: `{(1, 3), (5, 6), (2, 6), (2, 3), (3, 5), (1, 2), (4, 6), (1, 5), (4, 5), (1, 6), (3, 4), (2, 4)}`
    - :math:`G_{48}`
        - `graph_type='Max7vertices'`
        - edges: `{(2, 7), (4, 7), (1, 3), (4, 5), (1, 4), (5, 6), (2, 6), (1, 6), (3, 7), (1, 2), (6, 7), (5, 7), (1, 5), (2, 3), (3, 4)}`
    - :math:`G_{32a}`
        - `graph_type='7vert32a'`
        - edges: `{(4, 7), (1, 3), (5, 6), (1, 4), (1, 6), (3, 7), (2, 5), (3, 5), (1, 2), (6, 7), (5, 7), (3, 6), (2, 3), (3, 4), (2, 4)}`
    - :math:`G_{32b}`
        - `graph_type='7vert32b'`
        - edges: `{(2, 7), (4, 7), (2, 6), (4, 5), (1, 4), (5, 6), (1, 3), (2, 3), (3, 7), (2, 5), (1, 2), (6, 7), (1, 5), (3, 6), (3, 4)}`
    - :math:`G_{24}`
        - `graph_type='7vert24'`
        - edges: `{(2, 7), (4, 7), (2, 6), (5, 6), (1, 4), (1, 3), (2, 3), (3, 7), (2, 5), (1, 2), (4, 6), (5, 7), (1, 5), (3, 6), (3, 4)}`
    - :math:`G_{16a}`
        - `graph_type='7vert16a'`
        - edges: `{(4, 7), (1, 3), (5, 6), (1, 6), (3, 7), (2, 5), (3, 5), (1, 2), (4, 6), (5, 7), (3, 6), (1, 7), (2, 3), (3, 4), (2, 4)}`
    - :math:`G_{16b}`
        - `graph_type='7vert16b'`
        - edges: `{(2, 7), (4, 7), (2, 6), (4, 5), (1, 4), (1, 3), (2, 3), (3, 7), (2, 5), (3, 5), (1, 2), (6, 7), (4, 6), (1, 5), (3, 6)}`
    - :math:`G_{160}`
        - `graph_type='Max8vertices'`, or `'Max8vertices_distSyst'` for using distance system instead of sphere equations (faster but often inaccurate)
        - edges: `{(2, 7), (3, 2), (2, 6), (6, 8), (7, 8), (6, 1), (3, 1), (2, 8), (4, 7), (2, 1), (5, 8), (4, 3), (5, 1), (5, 4), (3, 7), (4, 1), (6, 5), (5, 7)}`
    - :math:`G_{128}`
        - `graph_type='Ring8vertices'`
        - edges: `{(1, 2), (2, 7), (5, 6), (1, 3), (6, 7), (6, 8), (4, 8), (4, 5), (2, 8), (7, 8), (1, 4), (3, 8), (1, 5), (1, 6), (1, 7), (2, 3), (3, 4), (5, 8)}`
    '''
#    .. image:: http://jan.legersky.cz/public_files/spatialGraphEmbeddings/graphs_7and8vert.png
#       :width: 70 %
#       :alt: Supported 7 and 8-vertex graphs

    def __init__(self, lengths, graph_type, tmpFileName=None,  window=None):
        '''
        Inputs:
        
        - `lengths` is a dictionary with edge lengths of graph given by `graph_type`
        - `tmpFileName` is used for temporary files used during computations. If `None`, random hash is used.
        
        For implementation of a new graph, method **constructEquations_newGraph** must be implemented and *__init__* modified accordingly.
        Moreover, constructor of :py:mod:`AlgRealEmbeddings` must be modified.
        
        Optionally, method **getEmbedding** and function **getEdgeLengthsByEmbedding** should be implemented.
        '''
        self._window = window
        self._prevSystem = None
        self._prevSolutions = None
        self._verbose = 1
        self._fixedTriangle_vertices = [2, 3, 1]
        self._vertexWithFootAtOrigin = None
        if graph_type == 'Max7vertices':
            self._vertexWithFootAtOrigin = 7
            self.constructEquations = self.constructEquations_max7vertices
            self._numAllSolutions = 48
        elif graph_type == 'Max6vertices':
            self._vertexWithFootAtOrigin = 4
            self.constructEquations = self.constructEquations_max6vertices
            self._numAllSolutions = 16
        elif graph_type == 'Max8vertices':
            self.constructEquations = self.constructEquations_max8vertices
            self._numAllSolutions = 160
        elif graph_type == 'Max8vertices_distSyst':
            self.constructEquations = self.constructEquations_max8vertices_distSystem
            self._numAllSolutions = 80
        elif graph_type == 'Ring8vertices':
            self.constructEquations = self.constructEquations_ring8vertices
            self._numAllSolutions = 128
        elif graph_type == '7vert32a':
            self.constructEquations = self.constructEquations_7vert32a
            self._numAllSolutions = 32
        elif graph_type == '7vert32b':
            self.constructEquations = self.constructEquations_7vert32b
            self._numAllSolutions = 32
        elif graph_type == '7vert24':
            self.constructEquations = self.constructEquations_7vert24
            self._numAllSolutions = 24
        elif graph_type == '7vert16a':
            self.constructEquations = self.constructEquations_7vert16a
            self._numAllSolutions = 16
        elif graph_type == '7vert16b':
            self.constructEquations = self.constructEquations_7vert16b
            self._numAllSolutions = 16
        else:
            raise ValueError('Type %s not supported' % graph_type)
        self._graph_type = graph_type
        if tmpFileName==None:
            from random import random
            hash_object = hashlib.md5(str(time.time()).encode()+str(random()))
            self._fileNamePref = graph_type+'_'+str(hash_object.hexdigest())
        else:
            self._fileNamePref = graph_type+'_'+str(tmpFileName)
        self.setLengths(lengths)    

    def setLengths(self, lengths):
        '''
        Set edge lengths to `lengths`.
        '''
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
        '''
        Return length of edge `uv`.
        '''
        if v==None:
            return self.getEdgeLength(int(u[0]), int(u[1]))
        else:
            if u<v:
                return float(self._lengths[(u, v)])
            else:
                return float(self._lengths[(v, u)])
    
    def setEdgeLength(self, Luv, u, v):
        '''
        Set length of edge `uv` to `Luv`.
        '''
        if Luv<=10e-6:
            raise ValueError('Length of '+str([u, v])+' cannot be set to '+str(Luv))
        if Luv>1e6:
            self._window.showError('Length'+str(Luv)+ ' of '+str([u, v])+' is too big')
            raise ValueError('Length'+str(Luv)+ ' of '+str([u, v])+' is too big')
        if u<v:
            self._lengths[(u, v)] = Luv
        else:
            self._lengths[(v, u)] = Luv
    
    def getLengths(self):
        '''
        Return dictionary of edge lengths.
        '''
        return self._lengths

    def printLog(self, s, verbose=0):
        if self._window:
            self._window.printLog(str(s), verbose=verbose)
        else:
            if verbose<= self._verbose:
                print s

    def getEquations(self):
        '''
        Return sphere equations of the graph corresponding to current edge lengths.
        '''
        return self.constructEquations()
    
    def getAltitudeAndFoot(self, u, v, w):
#        ''' 
#        Return altitude of triangle `uvw` from `w` and the distance of its foot from `v`.
#        '''
        Luv = self.getEdgeLength(u, v)
        Lvw = self.getEdgeLength(w, v)
        Luw = self.getEdgeLength(w, u)
        cos_alpha = (Lvw**2+Luv**2 - Luw**2)/float(Lvw*2*Luv)
        if cos_alpha>=-1 and cos_alpha<=1:
            return [Lvw*math.sin(math.acos(cos_alpha)), cos_alpha*Lvw]
        else:
            raise TriangleInequalityError('Altitude and foot for the triangle '+str([u, v, w])+' with lengths '+str([Luv, Lvw, Luw])+' is not defined.')
    
    def coordinatesOfTriangle(self, u, v, w, yshift=0):
#        '''
#        Return coordinates of the tringle uvw so that it lies in x-y plane, u,v are on y-axis, y-coord. of u is yshift and v is in positive direction from u
#        '''
        u_coor = [0, yshift, 0]
        v_coor = [0, yshift+self.getEdgeLength(u, v), 0]
        alt_w, foot_wv = self.getAltitudeAndFoot(u, v, w)
        w_coor = [alt_w, v_coor[1]-foot_wv, 0]
        return u_coor, v_coor, w_coor
    
    def updateFixedTriangle(self):
#        '''
#        Adjusts coordinates of the fixed triangle according to _lengths. 
#        If `self._vertexWithFootAtOrigin` is not `None`, then the coordinate system is shifted so that foot of altitude from p in the triangle uvp is in the origin. (uvp must be in the graph)
#        '''
        u, v, w = self._fixedTriangle_vertices
        if self._vertexWithFootAtOrigin!=None:
            yshift = -self.getAltitudeAndFoot(v, u, self._vertexWithFootAtOrigin)[1]
        else:
            yshift = 0
        self._fixedTriangle = self.coordinatesOfTriangle(u, v, w, yshift)
    
    def setEdgeLengthWithCorrespondingOnes(self, Luv,  uvwp):
#        '''Sets length of uv to Luv and also lengths of uw and up so that angles uvp and uvw preserves '''
        u, v, w, p = uvwp
        p_coord = self.coordinatesOfTriangle(u, v, p)[2]
        w_coord = self.coordinatesOfTriangle(u, v, w)[2]
        new_u = [0, self.getEdgeLength(u, v)-Luv, 0]
        self.setEdgeLength(Luv, u, v)
        self.setEdgeLength(self.dist(p_coord, new_u), p, u)
        self.setEdgeLength(self.dist(w_coord, new_u), w, u)

    def dist(self, u, v):
        return float(np.sqrt( (u[0]-v[0])**2 + (u[1]-v[1])**2 + (u[2]-v[2])**2))

    
    def findEmbeddings(self, tolerance=1.0e-15,  errorMsg=True, usePrev=True):
        '''
        Compute embeddings of the graph compatible with the current edge lengths and fixed triangle and return them as dictionary 
        `{['real']: listRealEmbeddings, ['complex']: listComplexEmbeddings}`. Embeddings are considered real if the imaginary part of all coordinates is smaller than `tolerance`.
        
        Package ``phcpy`` is used for the computation. If `usePrev=True`, then the solutions are tracked from ones from the previous call, if there was any.
        '''
        syst = self.getEquations()
        i = 0

        if self._graph_type == 'Max8vertices_distSyst':
            num_conjugated = 2
            try:
                y1_left,  y1_right, y4_left,  y4_right = self.inequalities_max8vertices()
            except TriangleInequalityError:
                self.printLog('intervals could not be computed')
                return {'real':[], 'complex':[None for _ in range(0, 80)]}
        else:
            num_conjugated = 4

        while True:
            if self._prevSystem and usePrev:
                filenameTmp = os.path.dirname(os.path.realpath(__file__))+'/../tmp/'+self._fileNamePref+'_prev.txt'
                with open(filenameTmp, 'w') as filePrev:
                    filePrev.write(str(self._prevSystem)+'\n')
                    filePrev.write(str(self._prevSolutions)+'\n')
                    
                process = subprocess.Popen(['python2',os.path.dirname(os.path.realpath(__file__))+'/track.py', str(syst), filenameTmp, self._fileNamePref])
                process.wait()
                with open(os.path.dirname(os.path.realpath(__file__))+'/../tmp/'+self._fileNamePref+'track.txt','r') as file:
                    solutions_str = file.read()
                sols = ast.literal_eval(solutions_str)
                
                os.remove(filenameTmp)
                os.remove(os.path.dirname(os.path.realpath(__file__))+'/../tmp/'+self._fileNamePref+'track.txt')

            else:
                process = subprocess.Popen(['python2',os.path.dirname(os.path.realpath(__file__))+'/solve.py', str(syst), self._fileNamePref])
                process.wait()
                with open(os.path.dirname(os.path.realpath(__file__))+'/../tmp/'+self._fileNamePref+'solve.txt','r') as file:
                    sols_str = file.read()
                sols = ast.literal_eval(sols_str)
                os.remove(os.path.dirname(os.path.realpath(__file__))+'/../tmp/'+self._fileNamePref+'solve.txt')

            result_real = []
            result_complex = []

            for sol in sols:
                soldic = strsol2dict(sol)
                if self._graph_type != 'Max8vertices_distSyst':
                    if is_real(sol, tolerance):
                        result_real.append(soldic)
                    else:
                        result_complex.append(soldic)
                else:
                    if is_real(sol, tolerance) and soldic['y1'].real>=y1_left and soldic['y1'].real<=y1_right and soldic['y4'].real>=y4_left and soldic['y4'].real<=y4_right:
                        result_real.append(soldic)
                    else:
                        result_complex.append(soldic)          

            if len(result_real)%num_conjugated==0 and len(sols)==self._numAllSolutions:
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
        '''
        Set edge lengths so that the angles :math:`\\phi` and :math:`\\theta` in the subgraph :math:`(u,v,w,p,c)` given by 5-tuple `uvwpc` are `phi` and `theta`.
        '''
        self.setPhi(uvwpc[:-1], phi)
        self.setTheta([uvwpc[i] for i in [0, 2, 4]], theta)

    def getPhiTheta(self, uvwpc):      
        '''
        Return angles :math:`\\phi` and :math:`\\theta` in the subgraph :math:`(u,v,w,p,c)` given by 5-tuple `uvwpc`.
        '''
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
        '''
        Return one of the real embeddings comaptible with the current edge lengths.
        '''
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
        elif self._graph_type=='Ring8vertices':
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
        else:
            raise NotImplementedError('Method getEmbedding is not implemented for graph '+ self._graph_type)

    def constructEquations_max7vertices(self):
#        '''system with correct mixed volume'''        
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
#        '''system with correct mixed volume'''        
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
#        '''system with correct mixed volume'''        
        y1, y2, y3, y4 = symbols('y1 y2 y3 y4')
        
        c12 = self.getEdgeLength('12')**2
        c13 = self.getEdgeLength('13')**2
        c14 = self.getEdgeLength('14')**2
        c15 = self.getEdgeLength('15')**2
        c16 = self.getEdgeLength('16')**2

        c27 = self.getEdgeLength('27')**2
        c37 = self.getEdgeLength('37')**2
        c47 = self.getEdgeLength('47')**2
        c57 = self.getEdgeLength('57')**2
        
        c28 = self.getEdgeLength('28')**2
        c58 = self.getEdgeLength('58')**2
        c68 = self.getEdgeLength('68')**2
        c78 = self.getEdgeLength('78')**2
        
        c23 = self.getEdgeLength('23')**2
        c34 = self.getEdgeLength('34')**2
        c45 = self.getEdgeLength('45')**2
        c56 = self.getEdgeLength('56')**2
        c26 = self.getEdgeLength('26')**2
        
       
        eqs = [
        (-c12**2*c34**2+2*c12**2*c34*c37+2*c12**2*c34*c47-c12**2*c37**2+2*c12**2*c37*c47-c12**2*c47**2+4*c12*c13*c23*c47-2*c12*c13*c27*c34
        +2*c12*c13*c27*c37-2*c12*c13*c27*c47-2*c12*c13*c34*c47+2*c12*c13*c34*y3-2*c12*c13*c37*c47-2*c12*c13*c37*y3+2*c12*c13*c47**2
        -2*c12*c13*c47*y3+2*c12*c14*c23*c34-2*c12*c14*c23*c37-2*c12*c14*c23*c47-2*c12*c14*c27*c34-2*c12*c14*c27*c37+2*c12*c14*c27*c47
        -2*c12*c14*c34*c37+2*c12*c14*c37**2-2*c12*c14*c37*c47+4*c12*c14*c37*y3-2*c12*c23*c34*c47-2*c12*c23*c34*y1-2*c12*c23*c37*c47
        +2*c12*c23*c37*y1+2*c12*c23*c47**2-2*c12*c23*c47*y1+2*c12*c27*c34**2-2*c12*c27*c34*c37-2*c12*c27*c34*c47+4*c12*c27*c34*y1
        +2*c12*c34**2*y1+4*c12*c34*c37*c47-2*c12*c34*c37*y1-2*c12*c34*c37*y3-2*c12*c34*c47*y1-2*c12*c34*y1*y3+2*c12*c37**2*y3
        -2*c12*c37*c47*y3-2*c12*c37*y1*y3+2*c12*c47*y1*y3-c13**2*c27**2+2*c13**2*c27*c47+2*c13**2*c27*y3-c13**2*c47**2+2*c13**2*c47*y3
        -c13**2*y3**2-2*c13*c14*c23*c27-2*c13*c14*c23*c47+2*c13*c14*c23*y3+2*c13*c14*c27**2+4*c13*c14*c27*c34-2*c13*c14*c27*c37-2*c13*c14*c27*c47
        -2*c13*c14*c27*y3+2*c13*c14*c37*c47-2*c13*c14*c37*y3-2*c13*c23*c27*c47+2*c13*c23*c27*y1+2*c13*c23*c47**2-2*c13*c23*c47*y1-2*c13*c23*c47*y3
        -2*c13*c23*y1*y3+2*c13*c27**2*c34-2*c13*c27*c34*c47-2*c13*c27*c34*y1-2*c13*c27*c34*y3-2*c13*c27*c37*y3+4*c13*c27*c47*y3-2*c13*c27*y1*y3
        +2*c13*c34*c47*y1-2*c13*c34*y1*y3-2*c13*c37*c47*y3+4*c13*c37*y1*y3+2*c13*c37*y3**2-2*c13*c47*y1*y3+2*c13*y1*y3**2-c14**2*c23**2
        +2*c14**2*c23*c27+2*c14**2*c23*c37-c14**2*c27**2+2*c14**2*c27*c37-c14**2*c37**2+2*c14*c23**2*c47+2*c14*c23**2*y1-2*c14*c23*c27*c34
        +4*c14*c23*c27*c37-2*c14*c23*c27*c47-2*c14*c23*c27*y1-2*c14*c23*c34*y1-2*c14*c23*c37*c47-2*c14*c23*c37*y1-2*c14*c23*c37*y3+4*c14*c23*c47*y1
        -2*c14*c23*y1*y3+2*c14*c27**2*c34-2*c14*c27*c34*c37-2*c14*c27*c34*y1-2*c14*c27*c37*y3+2*c14*c27*y1*y3+2*c14*c34*c37*y1+2*c14*c37**2*y3
        -2*c14*c37*y1*y3-c23**2*c47**2+2*c23**2*c47*y1-c23**2*y1**2+2*c23*c27*c34*c47-2*c23*c27*c34*y1-2*c23*c34*c47*y1+2*c23*c34*y1**2
        +4*c23*c34*y1*y3+2*c23*c37*c47*y3-2*c23*c37*y1*y3-2*c23*c47*y1*y3+2*c23*y1**2*y3-c27**2*c34**2+2*c27*c34**2*y1+2*c27*c34*c37*y3
        -2*c27*c34*y1*y3-c34**2*y1**2-2*c34*c37*y1*y3+2*c34*y1**2*y3-c37**2*y3**2+2*c37*y1*y3**2-y1**2*y3**2),
        (-c12**2*c45**2+2*c12**2*c45*c47+2*c12**2*c45*c57-c12**2*c47**2+2*c12**2*c47*c57-c12**2*c57**2-2*c12*c14*c27*c45+2*c12*c14*c27*c47
        -2*c12*c14*c27*c57-2*c12*c14*c45*c57+2*c12*c14*c45*y4-2*c12*c14*c47*c57-2*c12*c14*c47*y4+2*c12*c14*c57**2+4*c12*c14*c57*y3-2*c12*c14*c57*y4
        -2*c12*c15*c27*c45-2*c12*c15*c27*c47+2*c12*c15*c27*c57-2*c12*c15*c45*c47+2*c12*c15*c45*y3+2*c12*c15*c47**2-2*c12*c15*c47*c57-2*c12*c15*c47*y3
        +4*c12*c15*c47*y4-2*c12*c15*c57*y3+2*c12*c27*c45**2-2*c12*c27*c45*c47-2*c12*c27*c45*c57+4*c12*c27*c45*y1+2*c12*c45**2*y1+4*c12*c45*c47*c57
        -2*c12*c45*c47*y1-2*c12*c45*c47*y4-2*c12*c45*c57*y1-2*c12*c45*c57*y3-2*c12*c45*y1*y3-2*c12*c45*y1*y4+2*c12*c47**2*y4-2*c12*c47*c57*y3
        -2*c12*c47*c57*y4+2*c12*c47*y1*y3-2*c12*c47*y1*y4+2*c12*c57**2*y3-2*c12*c57*y1*y3+2*c12*c57*y1*y4-c14**2*c27**2+2*c14**2*c27*c57
        +2*c14**2*c27*y4-c14**2*c57**2+2*c14**2*c57*y4-c14**2*y4**2+2*c14*c15*c27**2+4*c14*c15*c27*c45-2*c14*c15*c27*c47-2*c14*c15*c27*c57
        -2*c14*c15*c27*y3-2*c14*c15*c27*y4+2*c14*c15*c47*c57-2*c14*c15*c47*y4-2*c14*c15*c57*y3+2*c14*c15*y3*y4+2*c14*c27**2*c45-2*c14*c27*c45*c57
        -2*c14*c27*c45*y1-2*c14*c27*c45*y4-2*c14*c27*c47*y4-2*c14*c27*c57*y3+4*c14*c27*c57*y4+2*c14*c27*y1*y3-2*c14*c27*y1*y4+2*c14*c45*c57*y1
        -2*c14*c45*y1*y4-2*c14*c47*c57*y4+4*c14*c47*y1*y4+2*c14*c47*y4**2+2*c14*c57**2*y3-2*c14*c57*y1*y3-2*c14*c57*y1*y4-2*c14*c57*y3*y4
        -2*c14*y1*y3*y4+2*c14*y1*y4**2-c15**2*c27**2+2*c15**2*c27*c47+2*c15**2*c27*y3-c15**2*c47**2+2*c15**2*c47*y3-c15**2*y3**2+2*c15*c27**2*c45
        -2*c15*c27*c45*c47-2*c15*c27*c45*y1-2*c15*c27*c45*y3+4*c15*c27*c47*y3-2*c15*c27*c47*y4-2*c15*c27*c57*y3-2*c15*c27*y1*y3+2*c15*c27*y1*y4+2*c15*c45*c47*y1
        -2*c15*c45*y1*y3+2*c15*c47**2*y4-2*c15*c47*c57*y3-2*c15*c47*y1*y3-2*c15*c47*y1*y4-2*c15*c47*y3*y4+4*c15*c57*y1*y3+2*c15*c57*y3**2+2*c15*y1*y3**2
        -2*c15*y1*y3*y4-c27**2*c45**2+2*c27*c45**2*y1+2*c27*c45*c47*y4+2*c27*c45*c57*y3-2*c27*c45*y1*y3-2*c27*c45*y1*y4-c45**2*y1**2-2*c45*c47*y1*y4
        -2*c45*c57*y1*y3+2*c45*y1**2*y3+2*c45*y1**2*y4+4*c45*y1*y3*y4-c47**2*y4**2+2*c47*c57*y3*y4-2*c47*y1*y3*y4+2*c47*y1*y4**2-c57**2*y3**2+2*c57*y1*y3**2
        -2*c57*y1*y3*y4-y1**2*y3**2+2*y1**2*y3*y4-y1**2*y4**2),
        (-c12**2*c56**2+2*c12**2*c56*c58+2*c12**2*c56*c68-c12**2*c58**2+2*c12**2*c58*c68-c12**2*c68**2+2*c12*c15*c26*c56-2*c12*c15*c26*c58
        -2*c12*c15*c26*c68-2*c12*c15*c28*c56+2*c12*c15*c28*c58-2*c12*c15*c28*c68-2*c12*c15*c56*c68-2*c12*c15*c58*c68+2*c12*c15*c68**2
        +4*c12*c15*c68*y4+4*c12*c16*c26*c58-2*c12*c16*c28*c56-2*c12*c16*c28*c58+2*c12*c16*c28*c68-2*c12*c16*c56*c58+2*c12*c16*c56*y4
        +2*c12*c16*c58**2-2*c12*c16*c58*c68-2*c12*c16*c58*y4-2*c12*c16*c68*y4-2*c12*c26*c56*c58-2*c12*c26*c56*y2+2*c12*c26*c58**2
        -2*c12*c26*c58*c68-2*c12*c26*c58*y2+2*c12*c26*c68*y2+2*c12*c28*c56**2-2*c12*c28*c56*c58-2*c12*c28*c56*c68+4*c12*c28*c56*y2
        +2*c12*c56**2*y2+4*c12*c56*c58*c68-2*c12*c56*c58*y2-2*c12*c56*c68*y2-2*c12*c56*c68*y4-2*c12*c56*y2*y4-2*c12*c58*c68*y4
        +2*c12*c58*y2*y4+2*c12*c68**2*y4-2*c12*c68*y2*y4-c15**2*c26**2+2*c15**2*c26*c28+2*c15**2*c26*c68-c15**2*c28**2+2*c15**2*c28*c68
        -c15**2*c68**2-2*c15*c16*c26*c28-2*c15*c16*c26*c58+2*c15*c16*c26*y4+2*c15*c16*c28**2+4*c15*c16*c28*c56-2*c15*c16*c28*c58-2*c15*c16*c28*c68
        -2*c15*c16*c28*y4+2*c15*c16*c58*c68-2*c15*c16*c68*y4+2*c15*c26**2*c58+2*c15*c26**2*y2-2*c15*c26*c28*c56-2*c15*c26*c28*c58+4*c15*c26*c28*c68
        -2*c15*c26*c28*y2-2*c15*c26*c56*y2-2*c15*c26*c58*c68+4*c15*c26*c58*y2-2*c15*c26*c68*y2-2*c15*c26*c68*y4-2*c15*c26*y2*y4+2*c15*c28**2*c56
        -2*c15*c28*c56*c68-2*c15*c28*c56*y2-2*c15*c28*c68*y4+2*c15*c28*y2*y4+2*c15*c56*c68*y2+2*c15*c68**2*y4-2*c15*c68*y2*y4-c16**2*c28**2
        +2*c16**2*c28*c58+2*c16**2*c28*y4-c16**2*c58**2+2*c16**2*c58*y4-c16**2*y4**2-2*c16*c26*c28*c58+2*c16*c26*c28*y2+2*c16*c26*c58**2
        -2*c16*c26*c58*y2-2*c16*c26*c58*y4-2*c16*c26*y2*y4+2*c16*c28**2*c56-2*c16*c28*c56*c58-2*c16*c28*c56*y2-2*c16*c28*c56*y4+4*c16*c28*c58*y4
        -2*c16*c28*c68*y4-2*c16*c28*y2*y4+2*c16*c56*c58*y2-2*c16*c56*y2*y4-2*c16*c58*c68*y4-2*c16*c58*y2*y4+4*c16*c68*y2*y4+2*c16*c68*y4**2+2*c16*y2*y4**2
        -c26**2*c58**2+2*c26**2*c58*y2-c26**2*y2**2+2*c26*c28*c56*c58-2*c26*c28*c56*y2-2*c26*c56*c58*y2+2*c26*c56*y2**2+4*c26*c56*y2*y4+2*c26*c58*c68*y4
        -2*c26*c58*y2*y4-2*c26*c68*y2*y4+2*c26*y2**2*y4-c28**2*c56**2+2*c28*c56**2*y2+2*c28*c56*c68*y4-2*c28*c56*y2*y4-c56**2*y2**2-2*c56*c68*y2*y4+2*c56*y2**2*y4
        -c68**2*y4**2+2*c68*y2*y4**2-y2**2*y4**2),
        (-c12**2*c57**2+2*c12**2*c57*c58+2*c12**2*c57*c78-c12**2*c58**2+2*c12**2*c58*c78-c12**2*c78**2+2*c12*c15*c27*c57-2*c12*c15*c27*c58-2*c12*c15*c27*c78
        -2*c12*c15*c28*c57+2*c12*c15*c28*c58-2*c12*c15*c28*c78-2*c12*c15*c57*c78-2*c12*c15*c58*c78+2*c12*c15*c78**2+4*c12*c15*c78*y4-2*c12*c27*c57*c58-2*c12*c27*c57*y2
        +2*c12*c27*c58**2-2*c12*c27*c58*c78+4*c12*c27*c58*y1-2*c12*c27*c58*y2+2*c12*c27*c78*y2+2*c12*c28*c57**2-2*c12*c28*c57*c58-2*c12*c28*c57*c78-2*c12*c28*c57*y1
        +4*c12*c28*c57*y2-2*c12*c28*c58*y1+2*c12*c28*c78*y1+2*c12*c57**2*y2+4*c12*c57*c58*c78-2*c12*c57*c58*y1-2*c12*c57*c58*y2-2*c12*c57*c78*y2-2*c12*c57*c78*y4
        +2*c12*c57*y1*y4-2*c12*c57*y2*y4+2*c12*c58**2*y1-2*c12*c58*c78*y1-2*c12*c58*c78*y4-2*c12*c58*y1*y4+2*c12*c58*y2*y4+2*c12*c78**2*y4-2*c12*c78*y1*y4
        -2*c12*c78*y2*y4-c15**2*c27**2+2*c15**2*c27*c28+2*c15**2*c27*c78-c15**2*c28**2+2*c15**2*c28*c78-c15**2*c78**2+2*c15*c27**2*c58+2*c15*c27**2*y2
        -2*c15*c27*c28*c57-2*c15*c27*c28*c58+4*c15*c27*c28*c78-2*c15*c27*c28*y1-2*c15*c27*c28*y2-2*c15*c27*c57*y2-2*c15*c27*c58*c78-2*c15*c27*c58*y1
        +4*c15*c27*c58*y2-2*c15*c27*c78*y2-2*c15*c27*c78*y4+2*c15*c27*y1*y4-2*c15*c27*y2*y4+2*c15*c28**2*c57+2*c15*c28**2*y1-2*c15*c28*c57*c78+4*c15*c28*c57*y1
        -2*c15*c28*c57*y2-2*c15*c28*c58*y1-2*c15*c28*c78*y1-2*c15*c28*c78*y4-2*c15*c28*y1*y4+2*c15*c28*y2*y4+2*c15*c57*c78*y2+2*c15*c58*c78*y1+2*c15*c78**2*y4
        -2*c15*c78*y1*y4-2*c15*c78*y2*y4-c27**2*c58**2+2*c27**2*c58*y2-c27**2*y2**2+2*c27*c28*c57*c58-2*c27*c28*c57*y2-2*c27*c28*c58*y1+2*c27*c28*y1*y2
        -2*c27*c57*c58*y2+2*c27*c57*y2**2+4*c27*c57*y2*y4+2*c27*c58**2*y1+2*c27*c58*c78*y4-2*c27*c58*y1*y2-2*c27*c58*y1*y4-2*c27*c58*y2*y4-2*c27*c78*y2*y4
        -2*c27*y1*y2*y4+2*c27*y2**2*y4-c28**2*c57**2+2*c28**2*c57*y1-c28**2*y1**2+2*c28*c57**2*y2-2*c28*c57*c58*y1+2*c28*c57*c78*y4-2*c28*c57*y1*y2
        -2*c28*c57*y1*y4-2*c28*c57*y2*y4+2*c28*c58*y1**2+4*c28*c58*y1*y4-2*c28*c78*y1*y4+2*c28*y1**2*y4-2*c28*y1*y2*y4-c57**2*y2**2+2*c57*c58*y1*y2
        -2*c57*c78*y2*y4-2*c57*y1*y2*y4+2*c57*y2**2*y4-c58**2*y1**2-2*c58*c78*y1*y4+2*c58*y1**2*y4-2*c58*y1*y2*y4-c78**2*y4**2+4*c78*y1*y2*y4+2*c78*y1*y4**2
        +2*c78*y2*y4**2-y1**2*y4**2+2*y1*y2*y4**2-y2**2*y4**2)        ]

        res = []
        for eq in eqs:
            res.append(str(eq)+';')
        return res

    def inequalities_max8vertices(self):
        c12 = self.getEdgeLength('12')**2
        c13 = self.getEdgeLength('13')**2
        c14 = self.getEdgeLength('14')**2
        c15 = self.getEdgeLength('15')**2
        c16 = self.getEdgeLength('16')**2

        c27 = self.getEdgeLength('27')**2
        c37 = self.getEdgeLength('37')**2
        c47 = self.getEdgeLength('47')**2
        c57 = self.getEdgeLength('57')**2
        
        c28 = self.getEdgeLength('28')**2
        c58 = self.getEdgeLength('58')**2
        c68 = self.getEdgeLength('68')**2
        c78 = self.getEdgeLength('78')**2
        
        c23 = self.getEdgeLength('23')**2
        c34 = self.getEdgeLength('34')**2
        c45 = self.getEdgeLength('45')**2
        c56 = self.getEdgeLength('56')**2
        c26 = self.getEdgeLength('26')**2
        
        triang_ineqs = [c12**2-2*c12*c13-2*c12*c23+c13**2-2*c13*c23+c23**2,
                        c12**2-2*c12*c16-2*c12*c26+c16**2-2*c16*c26+c26**2, 
                        c13**2-2*c13*c14-2*c13*c34+c14**2-2*c14*c34+c34**2, 
                        c14**2-2*c14*c15-2*c14*c45+c15**2-2*c15*c45+c45**2,
                        c15**2-2*c15*c16-2*c15*c56+c16**2-2*c16*c56+c56**2, 
                        c23**2-2*c23*c27-2*c23*c37+c27**2-2*c27*c37+c37**2, 
                        c34**2-2*c34*c37-2*c34*c47+c37**2-2*c37*c47+c47**2, 
                        c45**2-2*c45*c47-2*c45*c57+c47**2-2*c47*c57+c57**2]
        #all these equations should be <=0
        if max(triang_ineqs)>0:
            raise TriangleInequalityError('Triangle inequality violated!')
        
        try:
            cond_y1_left=[c12+c27-2*math.sqrt(c12*c27), c13+c37-2*math.sqrt(c13*c37), c14+c47-2*math.sqrt(c14*c47), c15+c57-2*math.sqrt(c15*c57),
                           -(1/4.0)*(-2*c12*c23+2*c12*c27-2*c12*c37-2*c13*c23-2*c13*c27+2*c13*c37+2*c23**2-2*c23*c27-2*c23*c37
                                   +2*math.sqrt(c12**2*c23**2-2*c12**2*c23*c27-2*c12**2*c23*c37+c12**2*c27**2-2*c12**2*c27*c37
                                                +c12**2*c37**2-2*c12*c13*c23**2+4*c12*c13*c23*c27+4*c12*c13*c23*c37-2*c12*c13*c27**2+4*c12*c13*c27*c37
                                                -2*c12*c13*c37**2-2*c12*c23**3+4*c12*c23**2*c27+4*c12*c23**2*c37-2*c12*c23*c27**2+4*c12*c23*c27*c37
                                                -2*c12*c23*c37**2+c13**2*c23**2-2*c13**2*c23*c27-2*c13**2*c23*c37+c13**2*c27**2-2*c13**2*c27*c37+c13**2*c37**2
                                                -2*c13*c23**3+4*c13*c23**2*c27+4*c13*c23**2*c37-2*c13*c23*c27**2+4*c13*c23*c27*c37-2*c13*c23*c37**2+c23**4
                                                -2*c23**3*c27-2*c23**3*c37+c23**2*c27**2-2*c23**2*c27*c37+c23**2*c37**2))/c23,
                            -(1/4.0)*(-2*c13*c34+2*c13*c37-2*c13*c47-2*c14*c34-2*c14*c37+2*c14*c47+2*c34**2-2*c34*c37-2*c34*c47
                                    +2*math.sqrt(c13**2*c34**2-2*c13**2*c34*c37-2*c13**2*c34*c47+c13**2*c37**2-2*c13**2*c37*c47+c13**2*c47**2
                                                 -2*c13*c14*c34**2+4*c13*c14*c34*c37+4*c13*c14*c34*c47-2*c13*c14*c37**2+4*c13*c14*c37*c47-2*c13*c14*c47**2
                                                 -2*c13*c34**3+4*c13*c34**2*c37+4*c13*c34**2*c47-2*c13*c34*c37**2+4*c13*c34*c37*c47-2*c13*c34*c47**2
                                                 +c14**2*c34**2-2*c14**2*c34*c37-2*c14**2*c34*c47+c14**2*c37**2-2*c14**2*c37*c47+c14**2*c47**2
                                                 -2*c14*c34**3+4*c14*c34**2*c37+4*c14*c34**2*c47-2*c14*c34*c37**2+4*c14*c34*c37*c47-2*c14*c34*c47**2
                                                 +c34**4-2*c34**3*c37-2*c34**3*c47+c34**2*c37**2-2*c34**2*c37*c47+c34**2*c47**2))/c34,
                            -(1/4.0)*(-2*c14*c45+2*c14*c47-2*c14*c57-2*c15*c45-2*c15*c47+2*c15*c57+2*c45**2-2*c45*c47-2*c45*c57
                                    +2*math.sqrt(c14**2*c45**2-2*c14**2*c45*c47-2*c14**2*c45*c57+c14**2*c47**2-2*c14**2*c47*c57+c14**2*c57**2
                                                 -2*c14*c15*c45**2+4*c14*c15*c45*c47+4*c14*c15*c45*c57-2*c14*c15*c47**2+4*c14*c15*c47*c57
                                                 -2*c14*c15*c57**2-2*c14*c45**3+4*c14*c45**2*c47+4*c14*c45**2*c57-2*c14*c45*c47**2+4*c14*c45*c47*c57
                                                 -2*c14*c45*c57**2+c15**2*c45**2-2*c15**2*c45*c47-2*c15**2*c45*c57+c15**2*c47**2-2*c15**2*c47*c57
                                                 +c15**2*c57**2-2*c15*c45**3+4*c15*c45**2*c47+4*c15*c45**2*c57-2*c15*c45*c47**2+4*c15*c45*c47*c57
                                                 -2*c15*c45*c57**2+c45**4-2*c45**3*c47-2*c45**3*c57+c45**2*c47**2-2*c45**2*c47*c57+c45**2*c57**2))/c45]

            cond_y1_right=[c12+c27+2*math.sqrt(c12*c27), c13+c37+2*math.sqrt(c13*c37), c14+c47+2*math.sqrt(c14*c47), c15+c57+2*math.sqrt(c15*c57), 
                           -(1/4.0)*(-2*c12*c23+2*c12*c27-2*c12*c37-2*c13*c23-2*c13*c27+2*c13*c37+2*c23**2-2*c23*c27-2*c23*c37
                                   -2*math.sqrt(c12**2*c23**2-2*c12**2*c23*c27-2*c12**2*c23*c37+c12**2*c27**2-2*c12**2*c27*c37+c12**2*c37**2-2*c12*c13*c23**2
                                                +4*c12*c13*c23*c27+4*c12*c13*c23*c37-2*c12*c13*c27**2+4*c12*c13*c27*c37-2*c12*c13*c37**2-2*c12*c23**3
                                                +4*c12*c23**2*c27+4*c12*c23**2*c37-2*c12*c23*c27**2+4*c12*c23*c27*c37-2*c12*c23*c37**2+c13**2*c23**2
                                                -2*c13**2*c23*c27-2*c13**2*c23*c37+c13**2*c27**2-2*c13**2*c27*c37+c13**2*c37**2-2*c13*c23**3+4*c13*c23**2*c27
                                                +4*c13*c23**2*c37-2*c13*c23*c27**2+4*c13*c23*c27*c37-2*c13*c23*c37**2+c23**4-2*c23**3*c27-2*c23**3*c37
                                                +c23**2*c27**2-2*c23**2*c27*c37+c23**2*c37**2))/c23,
                            -(1/4.0)*(-2*c13*c34+2*c13*c37-2*c13*c47-2*c14*c34-2*c14*c37+2*c14*c47+2*c34**2-2*c34*c37-2*c34*c47
                                    -2*math.sqrt(c13**2*c34**2-2*c13**2*c34*c37-2*c13**2*c34*c47+c13**2*c37**2-2*c13**2*c37*c47+c13**2*c47**2
                                                 -2*c13*c14*c34**2+4*c13*c14*c34*c37+4*c13*c14*c34*c47-2*c13*c14*c37**2+4*c13*c14*c37*c47-2*c13*c14*c47**2
                                                 -2*c13*c34**3+4*c13*c34**2*c37+4*c13*c34**2*c47-2*c13*c34*c37**2+4*c13*c34*c37*c47-2*c13*c34*c47**2
                                                 +c14**2*c34**2-2*c14**2*c34*c37-2*c14**2*c34*c47+c14**2*c37**2-2*c14**2*c37*c47+c14**2*c47**2-2*c14*c34**3
                                                 +4*c14*c34**2*c37+4*c14*c34**2*c47-2*c14*c34*c37**2+4*c14*c34*c37*c47-2*c14*c34*c47**2+c34**4-2*c34**3*c37
                                                 -2*c34**3*c47+c34**2*c37**2-2*c34**2*c37*c47+c34**2*c47**2))/c34,
                            -(1/4.0)*(-2*c14*c45+2*c14*c47-2*c14*c57-2*c15*c45-2*c15*c47+2*c15*c57+2*c45**2-2*c45*c47-2*c45*c57
                                    -2*math.sqrt(c14**2*c45**2-2*c14**2*c45*c47-2*c14**2*c45*c57+c14**2*c47**2-2*c14**2*c47*c57+c14**2*c57**2
                                                 -2*c14*c15*c45**2+4*c14*c15*c45*c47+4*c14*c15*c45*c57-2*c14*c15*c47**2+4*c14*c15*c47*c57-2*c14*c15*c57**2
                                                 -2*c14*c45**3+4*c14*c45**2*c47+4*c14*c45**2*c57-2*c14*c45*c47**2+4*c14*c45*c47*c57-2*c14*c45*c57**2
                                                 +c15**2*c45**2-2*c15**2*c45*c47-2*c15**2*c45*c57+c15**2*c47**2-2*c15**2*c47*c57+c15**2*c57**2
                                                 -2*c15*c45**3+4*c15*c45**2*c47+4*c15*c45**2*c57-2*c15*c45*c47**2+4*c15*c45*c47*c57-2*c15*c45*c57**2+c45**4
                                                 -2*c45**3*c47-2*c45**3*c57+c45**2*c47**2-2*c45**2*c47*c57+c45**2*c57**2))/c45]

            cond_y4_left=[c12+c15-2*math.sqrt(c12*c15), c26+c56-2*math.sqrt(c26*c56), c27+c57-2*math.sqrt(c27*c57),
                           -(1/4.0)*(2*c12*c15-2*c12*c16-2*c12*c56-2*c15*c16-2*c15*c26+2*c16**2-2*c16*c26-2*c16*c56+2*c26*c56
                                   +2*math.sqrt(c12**2*c15**2-2*c12**2*c15*c16-2*c12**2*c15*c56+c12**2*c16**2-2*c12**2*c16*c56+c12**2*c56**2-2*c12*c15**2*c16
                                                -2*c12*c15**2*c26+4*c12*c15*c16**2+4*c12*c15*c16*c26+4*c12*c15*c16*c56+4*c12*c15*c26*c56-2*c12*c16**3
                                                -2*c12*c16**2*c26+4*c12*c16**2*c56+4*c12*c16*c26*c56-2*c12*c16*c56**2-2*c12*c26*c56**2+c15**2*c16**2
                                                -2*c15**2*c16*c26+c15**2*c26**2-2*c15*c16**3+4*c15*c16**2*c26-2*c15*c16**2*c56-2*c15*c16*c26**2
                                                +4*c15*c16*c26*c56-2*c15*c26**2*c56+c16**4-2*c16**3*c26-2*c16**3*c56+c16**2*c26**2+4*c16**2*c26*c56+c16**2*c56**2-2*c16*c26**2*c56-2*c16*c26*c56**2+c26**2*c56**2))/c16, -(1/4)*(2*c26*c56-2*c26*c58-2*c26*c68-2*c28*c56+2*c28*c58-2*c28*c68-2*c56*c68-2*c58*c68+2*c68**2+2*math.sqrt(c26**2*c56**2-2*c26**2*c56*c58-2*c26**2*c56*c68+c26**2*c58**2-2*c26**2*c58*c68+c26**2*c68**2-2*c26*c28*c56**2+4*c26*c28*c56*c58+4*c26*c28*c56*c68-2*c26*c28*c58**2+4*c26*c28*c58*c68-2*c26*c28*c68**2-2*c26*c56**2*c68+4*c26*c56*c58*c68+4*c26*c56*c68**2-2*c26*c58**2*c68+4*c26*c58*c68**2-2*c26*c68**3+c28**2*c56**2-2*c28**2*c56*c58-2*c28**2*c56*c68+c28**2*c58**2-2*c28**2*c58*c68+c28**2*c68**2-2*c28*c56**2*c68+4*c28*c56*c58*c68+4*c28*c56*c68**2-2*c28*c58**2*c68+4*c28*c58*c68**2-2*c28*c68**3+c56**2*c68**2-2*c56*c58*c68**2-2*c56*c68**3+c58**2*c68**2-2*c58*c68**3+c68**4))/c68, -(1/4)*(2*c27*c57-2*c27*c58-2*c27*c78-2*c28*c57+2*c28*c58-2*c28*c78-2*c57*c78-2*c58*c78+2*c78**2+2*math.sqrt(c27**2*c57**2-2*c27**2*c57*c58-2*c27**2*c57*c78+c27**2*c58**2-2*c27**2*c58*c78+c27**2*c78**2-2*c27*c28*c57**2+4*c27*c28*c57*c58+4*c27*c28*c57*c78-2*c27*c28*c58**2+4*c27*c28*c58*c78-2*c27*c28*c78**2-2*c27*c57**2*c78+4*c27*c57*c58*c78+4*c27*c57*c78**2-2*c27*c58**2*c78+4*c27*c58*c78**2-2*c27*c78**3+c28**2*c57**2-2*c28**2*c57*c58-2*c28**2*c57*c78+c28**2*c58**2-2*c28**2*c58*c78+c28**2*c78**2-2*c28*c57**2*c78+4*c28*c57*c58*c78+4*c28*c57*c78**2-2*c28*c58**2*c78+4*c28*c58*c78**2-2*c28*c78**3+c57**2*c78**2-2*c57*c58*c78**2-2*c57*c78**3+c58**2*c78**2-2*c58*c78**3+c78**4))/c78]

            cond_y4_right=[c12+c15+2*math.sqrt(c12*c15), c26+c56+2*math.sqrt(c26*c56), c27+c57+2*math.sqrt(c27*c57), 
                           -(1/4.0)*(2*c12*c15-2*c12*c16-2*c12*c56-2*c15*c16-2*c15*c26+2*c16**2-2*c16*c26-2*c16*c56+2*c26*c56
                                   -2*math.sqrt(c12**2*c15**2-2*c12**2*c15*c16-2*c12**2*c15*c56+c12**2*c16**2-2*c12**2*c16*c56+c12**2*c56**2
                                                -2*c12*c15**2*c16-2*c12*c15**2*c26+4*c12*c15*c16**2+4*c12*c15*c16*c26+4*c12*c15*c16*c56+4*c12*c15*c26*c56
                                                -2*c12*c16**3-2*c12*c16**2*c26+4*c12*c16**2*c56+4*c12*c16*c26*c56-2*c12*c16*c56**2-2*c12*c26*c56**2+c15**2*c16**2
                                                -2*c15**2*c16*c26+c15**2*c26**2-2*c15*c16**3+4*c15*c16**2*c26-2*c15*c16**2*c56-2*c15*c16*c26**2+4*c15*c16*c26*c56
                                                -2*c15*c26**2*c56+c16**4-2*c16**3*c26-2*c16**3*c56+c16**2*c26**2+4*c16**2*c26*c56+c16**2*c56**2-2*c16*c26**2*c56
                                                -2*c16*c26*c56**2+c26**2*c56**2))/c16,
                            -(1/4.0)*(2*c26*c56-2*c26*c58-2*c26*c68-2*c28*c56+2*c28*c58-2*c28*c68-2*c56*c68-2*c58*c68+2*c68**2
                                    -2*math.sqrt(c26**2*c56**2-2*c26**2*c56*c58-2*c26**2*c56*c68+c26**2*c58**2-2*c26**2*c58*c68+c26**2*c68**2-2*c26*c28*c56**2
                                                 +4*c26*c28*c56*c58+4*c26*c28*c56*c68-2*c26*c28*c58**2+4*c26*c28*c58*c68-2*c26*c28*c68**2-2*c26*c56**2*c68
                                                 +4*c26*c56*c58*c68+4*c26*c56*c68**2-2*c26*c58**2*c68+4*c26*c58*c68**2-2*c26*c68**3+c28**2*c56**2-2*c28**2*c56*c58
                                                 -2*c28**2*c56*c68+c28**2*c58**2-2*c28**2*c58*c68+c28**2*c68**2-2*c28*c56**2*c68+4*c28*c56*c58*c68+4*c28*c56*c68**2
                                                 -2*c28*c58**2*c68+4*c28*c58*c68**2-2*c28*c68**3+c56**2*c68**2-2*c56*c58*c68**2-2*c56*c68**3+c58**2*c68**2-2*c58*c68**3+c68**4))/c68,
                            -(1/4.0)*(2*c27*c57-2*c27*c58-2*c27*c78-2*c28*c57+2*c28*c58-2*c28*c78-2*c57*c78-2*c58*c78+2*c78**2
                                    -2*math.sqrt(c27**2*c57**2-2*c27**2*c57*c58-2*c27**2*c57*c78+c27**2*c58**2-2*c27**2*c58*c78+c27**2*c78**2-2*c27*c28*c57**2+4*c27*c28*c57*c58
                                                 +4*c27*c28*c57*c78-2*c27*c28*c58**2+4*c27*c28*c58*c78-2*c27*c28*c78**2-2*c27*c57**2*c78+4*c27*c57*c58*c78+4*c27*c57*c78**2
                                                 -2*c27*c58**2*c78+4*c27*c58*c78**2-2*c27*c78**3+c28**2*c57**2-2*c28**2*c57*c58-2*c28**2*c57*c78+c28**2*c58**2-2*c28**2*c58*c78
                                                 +c28**2*c78**2-2*c28*c57**2*c78+4*c28*c57*c58*c78+4*c28*c57*c78**2-2*c28*c58**2*c78+4*c28*c58*c78**2-2*c28*c78**3+c57**2*c78**2
                                                 -2*c57*c58*c78**2-2*c57*c78**3+c58**2*c78**2-2*c58*c78**3+c78**4))/c78]

        except:
            raise TriangleInequalityError('Triangular or tetrangular inequality violated!')
        
        #y1 is embeddable iff Int_y1[0]<=y1<=Int_y1[1]
        #y4 is embeddable iff Int_y4[0]<=y4<=Int_y4[1]
        # the mechanism is embeddable if both y1/y4 are so
        # should we use small tolerance (eg <10^(-5))?
        return [max(cond_y1_left),min(cond_y1_right), max(cond_y4_left),min(cond_y4_right)]

    def constructEquations_ring8vertices(self):
#        '''system with correct mixed volume'''    
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
#        '''system with correct mixed volume'''        
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
    '''
    Exception raised if a tringle inequality is violated.
    '''
    def __init__(self, errorMsg):
        super(TriangleInequalityError, self).__init__(errorMsg)

def getEdgeLengthsByEmbedding(graph_type, vert_coordinates):
    '''
    Return edge lengths for `graph_type` obtained by taking corresponding distances of vertices given by `vert_coordinates`.
    '''
    def dist( u, v):
        return float(np.sqrt( (u[0]-v[0])**2 + (u[1]-v[1])**2 + (u[2]-v[2])**2))
    if graph_type=='Max6vertices':
        v1, v2, v3, v4, v5, v6 = vert_coordinates
        lengths = {(1, 2) : dist(v1,v2), 
            (2, 3) : dist(v2,v3), 
            (3, 4) : dist(v3,v4), 
            (4, 5) : dist(v4,v5), 
            (5, 6) : dist(v5,v6), 
            (1, 6) : dist(v1,v6), 
            (1, 3) : dist(v1,v3), 
            (2, 4) : dist(v2,v4), 
            (3, 5) : dist(v3,v5), 
            (4, 6) : dist(v4,v6), 
            (1, 5) : dist(v1,v5), 
            (2, 6) : dist(v2,v6)}
    elif graph_type=='7vert16a':
        v1, v2, v3, v4, v5, v6, v7 = vert_coordinates
        lengths = {(1, 2) : dist(v1,v2),
            (1, 3) : dist(v1,v3),
            (1, 6) : dist(v1,v6),
            (1, 7) : dist(v1,v7),
            (2, 3) : dist(v2,v3),
            (2, 4) : dist(v2,v4),
            (2, 5) : dist(v2,v5),
            (3, 4) : dist(v3,v4),
            (3, 5) : dist(v3,v5),
            (3, 6) : dist(v3,v6),
            (3, 7) : dist(v3,v7),
            (4, 6) : dist(v4,v6),
            (4, 7) : dist(v4,v7),
            (5, 6) : dist(v5,v6),
            (5, 7) : dist(v5,v7)}
    elif graph_type=='7vert16b':
        v1, v2, v3, v4, v5, v6, v7 = vert_coordinates     
        lengths = {(1, 2) : dist(v1,v2),
            (1, 3) : dist(v1,v3),
            (1, 4) : dist(v1,v4),
            (1, 5) : dist(v1,v5),
            (2, 3) : dist(v2,v3),
            (2, 5) : dist(v2,v5),
            (2, 6) : dist(v2,v6),
            (2, 7) : dist(v2,v7),
            (3, 5) : dist(v3,v5),
            (3, 6) : dist(v3,v6),
            (3, 7) : dist(v3,v7),
            (4, 5) : dist(v4,v5),
            (4, 6) : dist(v4,v6),
            (4, 7) : dist(v4,v7),
            (6, 7) : dist(v6,v7)}
    elif graph_type=='7vert24':
        v1, v2, v3, v4, v5, v6, v7 = vert_coordinates          
        lengths = {(1, 2) : dist(v1,v2),
            (1, 3) : dist(v1,v3),
            (1, 4) : dist(v1,v4),
            (1, 5) : dist(v1,v5),
            (2, 3) : dist(v2,v3),
            (2, 5) : dist(v2,v5),
            (2, 6) : dist(v2,v6),
            (2, 7) : dist(v2,v7),
            (3, 4) : dist(v3,v4),
            (3, 6) : dist(v3,v6),
            (3, 7) : dist(v3,v7),
            (4, 6) : dist(v4,v6),
            (4, 7) : dist(v4,v7),
            (5, 6) : dist(v5,v6),
            (5, 7) : dist(v5,v7)}
    elif graph_type=='7vert32a':
        v1, v2, v3, v4, v5, v6, v7 = vert_coordinates
        lengths = {(1, 2) : dist(v1,v2),
            (1, 3) : dist(v1,v3),
            (1, 4) : dist(v1,v4),
            (1, 6) : dist(v1,v6),
            (2, 3) : dist(v2,v3),
            (2, 4) : dist(v2,v4),
            (2, 5) : dist(v2,v5),
            (3, 4) : dist(v3,v4),
            (3, 5) : dist(v3,v5),
            (3, 6) : dist(v3,v6),
            (3, 7) : dist(v3,v7),
            (4, 7) : dist(v4,v7),
            (5, 6) : dist(v5,v6),
            (5, 7) : dist(v5,v7),
            (6, 7) : dist(v6,v7)}
    elif graph_type=='7vert32b':
        v1, v2, v3, v4, v5, v6, v7 = vert_coordinates        
        lengths = {(1, 2) : dist(v1,v2),
            (1, 3) : dist(v1,v3),
            (1, 4) : dist(v1,v4),
            (1, 5) : dist(v1,v5),
            (2, 3) : dist(v2,v3),
            (2, 5) : dist(v2,v5),
            (2, 6) : dist(v2,v6),
            (2, 7) : dist(v2,v7),
            (3, 4) : dist(v3,v4),
            (3, 6) : dist(v3,v6),
            (3, 7) : dist(v3,v7),
            (4, 5) : dist(v4,v5),
            (4, 7) : dist(v4,v7),
            (5, 6) : dist(v5,v6),
            (6, 7) : dist(v6,v7)}
    elif graph_type=='Max7vertices':
        v1, v2, v3, v4, v5, v6, v7 = vert_coordinates        
        lengths = {(1, 2) : dist(v1,v2),
            (1, 3) : dist(v1,v3),
            (1, 4) : dist(v1,v4),
            (1, 5) : dist(v1,v5),
            (1, 6) : dist(v1,v6),
            (2, 7) : dist(v7,v2),
            (3, 7) : dist(v7,v3),
            (4, 7) : dist(v7,v4),
            (5, 7) : dist(v7,v5),
            (6, 7) : dist(v7,v6),
            (2, 3) : dist(v2,v3), 
            (3, 4) : dist(v3,v4), 
            (4, 5) : dist(v4,v5), 
            (5, 6) : dist(v5,v6), 
            (2, 6) : dist(v2,v6)}
    elif graph_type=='Ring8vertices':
        v1, v2, v3, v4, v5, v6, v7, v8 = vert_coordinates        
        lengths = {
            (1, 2) : dist(v1,v2),
            (1, 3) : dist(v1,v3),
            (1, 4) : dist(v1,v4),
            (1, 5) : dist(v1,v5),
            (1, 6) : dist(v1,v6),
            (1, 7) : dist(v1,v7),
            (2, 3) : dist(v2,v3),
            (2, 7) : dist(v2,v7),
            (2, 8) : dist(v2,v8),
            (3, 4) : dist(v3,v4),
            (3, 8) : dist(v3,v8),
            (4, 5) : dist(v4,v5),
            (4, 8) : dist(v4,v8),
            (5, 6) : dist(v5,v6),
            (5, 8) : dist(v5,v8),
            (6, 8) : dist(v6,v8),
            (7, 8) : dist(v7,v8),
            (6, 7) : dist(v6,v7)}
    elif graph_type=='Max8vertices' or graph_type=='Max8vertices_distSyst':
        v1, v2, v3, v4, v5, v6, v7, v8 = vert_coordinates         
        lengths = {
            (2, 1) : dist(v2,v1),
            (2, 7) : dist(v2,v7),
            (2, 6) : dist(v2,v6),
            (3, 1) : dist(v3,v1),
            (3, 7) : dist(v3,v7),
            (3, 2) : dist(v3,v2),
            (4, 1) : dist(v4,v1),
            (4, 7) : dist(v4,v7),
            (4, 3) : dist(v4,v3),
            (5, 1) : dist(v5,v1),
            (5, 7) : dist(v5,v7),
            (5, 4) : dist(v5,v4),
            (6, 1) : dist(v6,v1),
            (6, 5) : dist(v6,v5),
            (5, 8) : dist(v5,v8),
            (6, 8) : dist(v6,v8),
            (7, 8) : dist(v7,v8),
            (2, 8) : dist(v2,v8)}
    else:
        raise NotImplementedError('Method getEdgeLengthsByEmbedding is not implemented for graph '+ graph_type)
    return lengths
