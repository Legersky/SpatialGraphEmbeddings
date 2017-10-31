import time
import math
import copy
import pickle
from sklearn.cluster import DBSCAN

from graphEmbedding import *
from graphCouplerCurve import *

class AlgRealEmbeddings(object):
    def __init__(self, lengths, fixedTriangle, window=None):
        self._window = window
        self.starting_lengths = lengths
        self._fixedTriangle = fixedTriangle

    def getSamplingIterator(self, num_phi, num_theta,  l_phi,  r_phi, l_theta,  r_theta):
        step_phi = (r_phi-l_phi)/float(num_phi)
        step_theta = (r_theta-l_theta)/float(num_theta)
        phi = l_phi
        
        margin_degree = 4
        margin = margin_degree*math.pi/180.0
        l_ext,  r_ext = -math.pi/2.0 + margin, math.pi/2.0 - margin
        while (phi<r_phi+step_phi):
            theta = l_theta
            while theta<r_theta+step_theta:
                if phi>=l_ext and phi<=r_ext:
                    yield phi, theta
                theta += step_theta
            phi += step_phi

    def runSamplingPhiThetaWithMargins(self, starting_graph, num_phi, num_theta,  l_phi,  r_phi, l_theta,  r_theta, uvwpc, treshold=0):
        return self.runSampling(starting_graph, self.getSamplingIterator(num_phi, num_theta,  l_phi,  r_phi, l_theta,  r_theta), uvwpc, treshold=treshold)

    def runSampling(self, starting_graph, iterator, uvwpc, treshold=0, animate=False):
        solutions = []
        if self._window and animate:
            N = self._window.spinBoxSamples.value()
            graphCouplerCurve = GraphCouplerCurve(self.starting_lengths, window=self._window)
            graphCouplerCurve.computeCouplerCurve(N)
            self._window.setActiveGraph(graphCouplerCurve)

        for phi, theta in iterator:
            try:
                starting_graph.setPhiTheta(uvwpc, phi, theta)
                num_real = len(starting_graph.findEmbeddings(errorMsg=False)['real'])
            except Exception as e:
                self.printLog(str(e))
                num_real = 0
            
            if self._window and animate:
                graphCouplerCurve.setPhiTheta([2, 3, 1, 7, 6], phi, theta)
                graphCouplerCurve.updateFixedTriangle()
                self._window.update_graph2phi()
                self._window.update_graph2theta()
                graphCouplerCurve.computeIntersections()
                self._window.plotScene()
                
            self.printLog(str([phi, theta, num_real]), verbose=1)

            if num_real>=treshold:
                solutions.append([phi, theta, num_real])              
        return solutions

    def computeSamplingPhiTheta(self, starting_lengths, num_phi, num_theta, uvwpc):
        starting_graph = GraphEmbedding(starting_lengths, self._fixedTriangle, window=self._window)
        
        start = time.time()
        act_num = len(starting_graph.findEmbeddings()['real'])
        act_phi, act_theta = starting_graph.getPhiTheta(uvwpc)
       
        margin_degree = 5
        margin = margin_degree*math.pi/180.0
        l_phi,  r_phi = -math.pi/2.0 + margin, math.pi/2.0 - margin
        l_theta, r_theta = 0 + margin/2.0, math.pi - margin/2.0
        
        sols = self.runSamplingPhiThetaWithMargins(starting_graph, num_phi, num_theta,  l_phi,  r_phi, l_theta,  r_theta, uvwpc, treshold=act_num)
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
            for phi2, theta2, num2 in self.runSamplingPhiThetaWithMargins(starting_graph, num_phi/4, num_theta/4,
                                                                          phi - step_phi,  phi + step_phi,
                                                                          theta - step_theta,  theta + step_theta, uvwpc, 
                                                                          treshold=maximum):
                if num2>maximum:
                    maximum = num2
                    argmax = []
                    self.printLog('Maximum increased to '+str(num2), verbose=1)
                    argmax.append([phi2, theta2])
                elif num2==maximum:
                    argmax.append([phi2, theta2])
                    self.printLog(str(num2), verbose=2)
        
        if len(argmax)==1:
            clusters = [argmax]
        else:
            self.printLog('Clustering phase starts:')
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
        
        res_lengths= []
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
            phi_c, theta_c = s_x/float(len(cluster)), s_y/float(len(cluster))
            starting_graph.setPhiTheta(uvwpc, phi_c, theta_c)
            n = len(starting_graph.findEmbeddings()['real'])
            if n < maximum:
                min_dist = dist_angle([phi_c, theta_c], cluster[0])
                phi_tmp = cluster[0][0]
                theta_tmp = cluster[0][1]
                for x, y in cluster:
                    d = dist_angle([phi_c, theta_c], [x, y])
                    if d < min_dist:
                        min_dist = d
                        phi_tmp = x
                        theta_tmp = y
                self.printLog('Center of cluster does not have maximum number of solutions \n -> nearest point chosen instead.')
                phi_c = phi_tmp
                theta_c = theta_tmp
                starting_graph.setPhiTheta(uvwpc, phi_c, theta_c)

            centers.append([phi_c, theta_c])
            
            res_lengths.append(copy.copy(starting_graph.getLengths()))
            res_infos.append(str([phi_c, theta_c]) + '\n embeddings: '+str(maximum))
        
#        if maximum<act_num:
#            res_lengths= [lengths]
#            self.printLog('original kept:\n'+str([act_phi, act_theta]))
#            res_infos = ['original kept:\n'+str([act_phi, act_theta])]
#            centers = [[act_phi, act_theta]]
#            maximum = act_num

        end = time.time()
        self.printLog('Maximum number of embeddings:')
        self.printLog(str(maximum))

        self.printLog('time: '+str(end - start))
        
        return [clusters, centers, res_lengths, res_infos, maximum]

    def runSamplingPhiTheta(self, starting_lengths, num_phi, num_theta, uvwpc):
        [clusters, centers, res_lengths, res_infos, maximum] = self.computeSamplingPhiTheta(starting_lengths, num_phi, num_theta, uvwpc)        
        if self._window:
            self._window.showClusters(clusters, centers)
            self._window.setGraphSequence([GraphCouplerCurve(lengths, window=self._window) for lengths in res_lengths], res_infos)

    def findMoreEmbeddings(self, starting_lengths, num_phi, num_theta, combinations,  prev_max, required_num, name=''):
        self._max_found = False
        self._num_phi = num_phi
        self._num_theta = num_theta
        self._combinations = combinations
        self._required_num = required_num
        self._reachedMaxs = []
        self._actMaximum = 0
        
        previous_steps = [['', 0]]
        previous_lengths = [copy.copy(starting_lengths)]
        self.findMoreEmbeddings_recursion(starting_lengths, previous_steps, previous_lengths, prev_max)
        
        
        report  = ['Reached maximum: ', self._actMaximum, 'Steps and lengths:', self._reachedMaxs]
        for r in report:
            self.printLog(r)                                                
        if self._window:
            self._window.showDialog(report)
        
        import hashlib
        hash_object = hashlib.md5(str(self._reachedMaxs).encode())
        
        fileName = './res/'+name+'_'+str(self._actMaximum)+'_embd_'+str(hash_object.hexdigest())+'.p'
        self.printLog('Results saved to: '+fileName)
        pickle.dump([self._actMaximum, self._reachedMaxs], open(fileName, 'wb'))

    def findMoreEmbeddings_recursion(self, starting_lengths, previous_steps, previous_lengths, prev_max):
        for uvwpc in self._combinations:
            if previous_steps[-1][0] != uvwpc:
                self.printLog('Actual step: ' + str(uvwpc))
                
                lengths_tmp = copy.copy(starting_lengths)
                [clusters, centers, res_lengths, res_infos, maximum] = self.computeSamplingPhiTheta(lengths_tmp, self._num_phi, self._num_theta, uvwpc)

                cl = 0                
                for lengths in res_lengths:
                    if not self._max_found:
                        cl +=1
                        steps = previous_steps+[[uvwpc, cl]]
                        all_lengths = previous_lengths + [copy.copy(lengths)]
                        
                        if maximum>self._actMaximum:
                            self._reachedMaxs = []
                            self._actMaximum = maximum
                        
                        if maximum==self._actMaximum:
                            self._reachedMaxs.append([steps, all_lengths])
                        
                        if maximum == self._required_num:
                            report  = [
                                        str(self._required_num)+' EMBEDDINGS FOUND:', 
                                        'Applied steps:',
                                        steps[1:],
                                        'All lengths:',
                                        all_lengths
                                         ]
                            for r in report:
                                self.printLog(r)                                                
                            if self._window:
                                self._window.showDialog(report)
                            self._max_found = True
                        
                        elif maximum>prev_max:
                            report  = [
                                        'MAXIMUM INCREASED to '+str(maximum), 
                                        'Applied steps:',
                                        steps[1:],
                                        'New lengths:',
                                        lengths
                                         ]
                            for r in report:
                                self.printLog(r)                                                
                            if self._window:
                                self._window.showDialog(report)
                            self.findMoreEmbeddings_recursion(lengths, steps, all_lengths, maximum)

    def printLog(self, s, verbose=0):
        if self._window:
            self._window.printLog(str(s), verbose=verbose)
        else:
            print s
