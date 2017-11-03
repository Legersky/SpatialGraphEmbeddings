import time
import math
import copy
#import pickle
from sklearn.cluster import DBSCAN
import hashlib
import sys
from random import random
import os
#import memory_profiler

from graphEmbedding import *
from graphCouplerCurve import *

class AlgRealEmbeddings(object):
    def __init__(self, graph_type, num_phi=20, num_theta=20, window=None, name=None):
        self._window = window
        self._graph_type = graph_type
        if graph_type == 'Max7vertices':
            self._numAllSol = 48
            self._combinations = [[2, 3, 1, 7, 6], [3, 7, 2, 4, 1], [5, 6, 1, 7, 4], [6, 1, 5, 2, 7], [3, 4, 1, 7, 2], [6, 7, 5, 2, 1], [4, 3, 1, 7, 5],
                                  [5, 1, 4, 6, 7], [2, 6, 1, 7, 3], [5, 7, 4, 6, 1], [3, 2, 1, 7, 4], [4, 1, 3, 5, 7], [6, 5, 1, 7, 2], [2, 7, 6, 3, 1],
                                  [3, 1, 2, 4, 7], [5, 4, 1, 7, 6], [4, 7, 3, 5, 1], [6, 2, 1, 7, 5], [2, 1, 6, 3, 7], [4, 5, 1, 7, 3]]
#            self._combinations = [[5, 6, 1, 7, 4], [4, 5, 1, 7, 3], [5, 7, 4, 6, 1], [2, 6, 1, 7, 3], [2, 7, 6, 3, 1], [3, 1, 2, 4, 7], [4, 1, 3, 5, 7],
#                                  [5, 4, 1, 7, 6], [4, 3, 1, 7, 5], [6, 7, 5, 2, 1], [4, 7, 3, 5, 1], [5, 1, 4, 6, 7], [3, 4, 1, 7, 2], [2, 1, 6, 3, 7],
#                                  [6, 5, 1, 7, 2], [3, 2, 1, 7, 4], [2, 3, 1, 7, 6], [6, 2, 1, 7, 5], [6, 1, 5, 2, 7], [3, 7, 2, 4, 1]]
        elif graph_type == 'Max6vertices':
            self._numAllSol = 16
            self._combinations =  [[2, 6, 1, 4, 3], [5, 4, 3, 6, 1], [5, 3, 1, 4, 6], [2, 1, 6, 3, 4], [5, 6, 1, 4, 3], [2, 4, 6, 3, 1], [6, 2, 1, 4, 5], [3, 1, 2, 5, 4],
                                   [6, 5, 1, 4, 2], [3, 4, 2, 5, 1], [2, 3, 1, 4, 6], [5, 1, 3, 6, 4], [3, 5, 1, 4, 2], [6, 4, 5, 2, 1], [3, 2, 1, 4, 5], [6, 1, 5, 2, 4]]
#            self._combinations = [[3, 2, 1, 4, 5], [6, 1, 5, 2, 4], [2, 4, 6, 3, 1], [2, 6, 1, 4, 3], [5, 1, 3, 6, 4], [5, 4, 3, 6, 1], [5, 3, 1, 4, 6], [5, 6, 1, 4, 3],
#                                  [6, 5, 1, 4, 2], [2, 3, 1, 4, 6], [3, 5, 1, 4, 2], [6, 4, 5, 2, 1], [3, 4, 2, 5, 1], [6, 2, 1, 4, 5], [3, 1, 2, 5, 4], [2, 1, 6, 3, 4]]
        elif graph_type == 'Max8vertices':
            self._numAllSol = 160
            self._combinations = [[3, 7, 2, 4, 1], [6, 8, 2, 5, 1], [4, 5, 1, 7, 3], [8, 7, 2, 5, 6], [4, 3, 1, 7, 5], [6, 2, 1, 8, 5], [3, 2, 1, 7, 4], [8, 5, 7, 6, 2], 
                                  [4, 1, 3, 5, 7], [6, 5, 1, 8, 2], [3, 1, 2, 4, 7], [8, 2, 7, 6, 5], [4, 7, 3, 5, 1], [6, 1, 2, 5, 8], [3, 4, 1, 7, 2], [8, 6, 2, 5, 7]]
        else:
            raise ValueError('Type %s not supported' % graph_type)
        
        hash_object = hashlib.md5(str(time.time()).encode()+str(random()))
        if name:
            self._fileNamePref = name+'_'+str(hash_object.hexdigest())[0:8]
        else:
            self._fileNamePref = graph_type+'_'+str(hash_object.hexdigest())[0:8]
        
        self._num_phi = num_phi
        self._num_theta = num_theta
        self._num_phi_second = num_phi/8
        self._num_theta_second = num_theta/8
        
        self.verbose = 1

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

    def runSampling_slower(self, starting_graph, iterator, uvwpc, treshold=0, animate=False):
        solutions = []
        if self._window and animate:
            N = self._window.spinBoxSamples.value()
            graphCouplerCurve = GraphCouplerCurve(starting_graph.getLengths(), window=self._window)
            graphCouplerCurve.computeCouplerCurve(N)
            self._window.setActiveGraph(graphCouplerCurve)

        self.printLog('...', newLine=False)
        
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

            if num_real>=treshold:
                self.printLog(str([phi, theta, num_real]), verbose=2)
                solutions.append([phi, theta, num_real])              
        return solutions

    def runSampling(self, starting_graph, iterator, uvwpc, treshold=0, animate=False):
        if self._window and animate:
            N = self._window.spinBoxSamples.value()
            graphCouplerCurve = GraphCouplerCurve(starting_graph.getLengths(), window=self._window)
            graphCouplerCurve.computeCouplerCurve(N)
            self._window.setActiveGraph(graphCouplerCurve)
        self.printLog('...', newLine=False)
        eqs = []
        phiThetas = []
        for phi, theta in iterator:
            try:
                starting_graph.setPhiTheta(uvwpc, phi, theta)
                eqs.append(str(starting_graph.getEquations()))
                phiThetas.append([phi, theta])
            except ValueError:
                pass
#            num_real = len(starting_graph.findEmbeddings(errorMsg=False)['real'])
            
            if self._window and animate:
                graphCouplerCurve.setPhiTheta([2, 3, 1, 7, 6], phi, theta)
                graphCouplerCurve.updateFixedTriangle()
                self._window.update_graph2phi()
                self._window.update_graph2theta()
                graphCouplerCurve.computeIntersections()
                self._window.plotScene()
        
        
        
        filenameTmp = 'tmp/'+self._fileNamePref
        with open(filenameTmp+'_prev.txt', 'w') as filePrev:
            filePrev.write(str(starting_graph._prevSystem)+'\n')
            filePrev.write(str(starting_graph._prevSolutions)+'\n')
            filePrev.write(str(self._numAllSol)+'\n')
        
        N = len(eqs)/2
        with open(filenameTmp+'_1_eqs.txt','w') as fileEqs:
            fileEqs.write('\n'.join(eqs[0:N]))
        
        process_1 = subprocess.Popen(['python', 'numReal.py', filenameTmp+'_1_', filenameTmp+'_prev.txt'])
        
        with open(filenameTmp+'_2_eqs.txt','w') as fileEqs:
            fileEqs.write('\n'.join(eqs[N:]))
        
        process_2 = subprocess.Popen(['python', 'numReal.py', filenameTmp+'_2_', filenameTmp+'_prev.txt'])
         
        
        process_1.wait()
        process_2.wait()
        
        with open(filenameTmp+'_1_'+'numReal.txt','r') as file:
            nums_real_str_1 = file.read()
        
        nums_real_1 = ast.literal_eval(nums_real_str_1)
        
        with open(filenameTmp+'_2_'+'numReal.txt','r') as file:
            nums_real_str_2 = file.read()
        
        nums_real_2 = ast.literal_eval(nums_real_str_2)
        
        
        solutions = []
        for i, num_real in enumerate(nums_real_1+nums_real_2):
            if num_real>=treshold:
                self.printLog(str(phiThetas[i]+ [num_real]), verbose=2)
                solutions.append(phiThetas[i]+ [num_real])              
        return solutions

    def runSampling_nonparallel(self, starting_graph, iterator, uvwpc, treshold=0, animate=False):
        if self._window and animate:
            N = self._window.spinBoxSamples.value()
            graphCouplerCurve = GraphCouplerCurve(starting_graph.getLengths(), window=self._window)
            graphCouplerCurve.computeCouplerCurve(N)
            self._window.setActiveGraph(graphCouplerCurve)
        self.printLog('...', newLine=False)
        eqs = []
        phiThetas = []
        for phi, theta in iterator:
            try:
                starting_graph.setPhiTheta(uvwpc, phi, theta)
                eqs.append(str(starting_graph.getEquations()))
                phiThetas.append([phi, theta])
            except TriangleInequalityError:
                pass
#            num_real = len(starting_graph.findEmbeddings(errorMsg=False)['real'])
            
            if self._window and animate:
                graphCouplerCurve.setPhiTheta([2, 3, 1, 7, 6], phi, theta)
                graphCouplerCurve.updateFixedTriangle()
                self._window.update_graph2phi()
                self._window.update_graph2theta()
                graphCouplerCurve.computeIntersections()
                self._window.plotScene()
        
        with open('tmp/'+self._fileNamePref+'eqs.txt','w') as fileEqs:
            fileEqs.write('\n'.join(eqs))
        
        filenameTmp = 'tmp/'+self._fileNamePref
        with open(filenameTmp+'_prev.txt', 'w') as filePrev:
            filePrev.write(str(starting_graph._prevSystem)+'\n')
            filePrev.write(str(starting_graph._prevSolutions)+'\n')
            filePrev.write(str(self._numAllSol)+'\n')
            
        process = subprocess.Popen(['python', 'numReal.py', 'tmp/'+self._fileNamePref, filenameTmp+'_prev.txt'])
        process.wait()
        
        with open('tmp/'+self._fileNamePref+'numReal.txt','r') as file:
            nums_real_str = file.read()
        
        nums_real = ast.literal_eval(nums_real_str)
        
        solutions = []
        for i, num_real in enumerate(nums_real):
            if num_real>=treshold:
                self.printLog(str(phiThetas[i]+ [num_real]), verbose=2)
                solutions.append(phiThetas[i]+ [num_real])              
        return solutions
    
    def sampleToGetMoreEmbd(self, starting_lengths, uvwpc):
        start = time.time()
        starting_graph = GraphEmbedding(copy.copy(starting_lengths), self._graph_type, window=self._window, tmpFileName=self._fileNamePref)
        argmax, maximum = self.computeSamplingPhiTheta(starting_graph, uvwpc)
        res = self.clusterPhiTheta(argmax, maximum, starting_graph, uvwpc)
        
        end = time.time()
        self.printLog('time: '+str(end - start))
        return res
        
    def computeSamplingPhiTheta(self, starting_graph, uvwpc):        
        act_num = len(starting_graph.findEmbeddings()['real'])
        act_phi, act_theta = starting_graph.getPhiTheta(uvwpc)
       
        margin_degree = 5
        margin = margin_degree*math.pi/180.0
        l_phi,  r_phi = -math.pi/2.0 + margin, math.pi/2.0 - margin
        l_theta, r_theta = 0 + margin/2.0, math.pi - margin/2.0
        
        sols = self.runSamplingPhiThetaWithMargins(starting_graph, self._num_phi, self._num_theta,  l_phi,  r_phi, l_theta,  r_theta, uvwpc, treshold=act_num)
        sols.append([act_phi, act_theta, act_num])

        self.printLog('Maximum number of embeddings in 1st round:')
        maximum = max([num for phi, theta, num in sols])
        self.printLog(str(maximum))
        max_positions = [ [phi, theta, num] for phi, theta, num in sols if num==maximum]
        
        step_phi = (r_phi-l_phi)/float(2*self._num_phi)
        step_theta = (r_theta-l_theta)/float(2*self._num_theta)
        argmax = []
        
        for phi, theta, num in max_positions:
            if num>=maximum:
                argmax.append([phi, theta])
            for phi2, theta2, num2 in self.runSamplingPhiThetaWithMargins(starting_graph, self._num_phi_second, self._num_theta_second,
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
        
        return argmax, maximum
        
    def clusterPhiTheta(self, argmax, maximum, starting_graph, uvwpc):
        if len(argmax)==1:
            clusters = [argmax]
        else:
            self.printLog('\nClustering phase')
            eps = 0.1
            for i in range(0, 100):
                labels = DBSCAN(eps=eps).fit_predict(argmax)
                if len([1 for el in labels if el==-1])>len(labels)/10:
                    eps += 0.05
                    self.printLog('Increasing eps for clustering', 2)
                else:
                    break
                if i==99:
                    self.printLog('Clustering was not succesfull')
                    labels = [0 for el in labels]
            
            clusters = [[] for i in range(0, max(list(labels)+[0])+1)]
            
            for label, point in zip(labels, argmax):
                if label>=0:
                    clusters[label].append(point)

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
            try:
                n = len(starting_graph.findEmbeddings()['real'])
            except ValueError:
                n = 0
            if n < maximum or n==0:
                phi_tmp = cluster[0][0]
                theta_tmp = cluster[0][1]
                min_dist = dist_angle([phi_c, theta_c], [phi_tmp,  theta_tmp])
                for x, y in cluster:
                    d = dist_angle([phi_c, theta_c], [x, y])
                    if d < min_dist:
                        min_dist = d
                        phi_tmp = x
                        theta_tmp = y
                self.printLog('Center of cluster does not have maximum number of solutions \n -> nearest point chosen instead.')
                phi_c = phi_tmp
                theta_c = theta_tmp
                

            centers.append([phi_c, theta_c])
            starting_graph.setPhiTheta(uvwpc, phi_c, theta_c)
            res_lengths.append(copy.copy(starting_graph.getLengths()))
            res_infos.append(str([phi_c, theta_c]) + '\n embeddings: '+str(maximum))
        self.printLog('Maximum number of embeddings:')
        self.printLog(str(maximum))
        
        return [clusters, centers, res_lengths, res_infos, maximum]

    def runSamplingPhiTheta(self, starting_lengths, uvwpc):
        [clusters, centers, res_lengths, res_infos, maximum] = self.sampleToGetMoreEmbd(starting_lengths, uvwpc)        
        if self._window:
            self._window.showClusters(clusters, centers)
            self._window.setGraphSequence([GraphCouplerCurve(lengths, window=self._window) for lengths in res_lengths], res_infos)
        return [clusters, centers, res_lengths, res_infos, maximum]

    def findMoreEmbeddings_tree(self, starting_lengths, required_num=None, onlyOne=True):
        self._max_found = False

        if required_num==None:
            self._required_num = self._numAllSol
        else:
            self._required_num = required_num
#        self._reachedMaxs = []
        self._onlyOne = onlyOne
        
        G = GraphEmbedding(starting_lengths, self._graph_type, tmpFileName=self._fileNamePref)
        sols = G.findEmbeddings()
        self._actMaximum = len(sols['real'])
        fromStr = 'from_'+str(self._actMaximum)
        self.printLog('Finding more embeddings - tree search')
        self.printLog('File name: '+self._fileNamePref+'\nStarting lengths:')
        self.printLog(starting_lengths)
        self.printLog(str(self._actMaximum)+ ' embeddings\n\n')

        previous_steps = [['', 0]]
        previous_lengths = [copy.copy(starting_lengths)]
        self.findMoreEmbeddings_recursion(starting_lengths, previous_steps, previous_lengths, self._actMaximum)
        
        outputFilename = './res/'+str(self._actMaximum)+'_embd_'+fromStr+'_'+self._fileNamePref+'.txt'
        os.rename('./res/'+self._fileNamePref+'_intermediateResults.txt', outputFilename)
        
        self.printLog('Reached maximums are in:')
        self.printLog(outputFilename)

    def findMoreEmbeddings_recursion(self, starting_lengths, previous_steps, previous_lengths, prev_max):
        try:
            ind = self._combinations.index(previous_steps[-1][0])
            comb = (self._combinations[ind:]+self._combinations[:ind])[1:]
        except:
            comb = self._combinations
        for uvwpc in comb:
            if not self._max_found:
                self.printLog('Working file name: '+self._fileNamePref+'_intermediateResults.txt')
                self.printLog('Reached maximum: '+str(self._actMaximum))
                self.printLog('Previous steps: '+str(previous_steps[1:]))
                self.printLog('Actual step: ' + str(uvwpc))

                
                lengths_tmp = copy.copy(starting_lengths)
                [clusters, centers, res_lengths, res_infos, maximum] = self.sampleToGetMoreEmbd(lengths_tmp, uvwpc)

                cl = 0                
                for lengths in res_lengths:
                    if not self._max_found:
                        cl +=1
                        steps = previous_steps+[[uvwpc, cl]]
                        all_lengths = previous_lengths + [copy.copy(lengths)]
                        
                        if maximum>self._actMaximum:
                            self._actMaximum = maximum
                            with open('./res/'+self._fileNamePref+'_intermediateResults.txt', 'w') as f:
                                f.write('')
                        
                        if maximum==self._actMaximum:
                            with open('./res/'+self._fileNamePref+'_intermediateResults.txt', 'a') as f:
                                report  = [
                                        str(maximum)+'\n', 
                                        str(steps[1:])+'\n',
                                        str(all_lengths)+'\n\n'
                                         ]
                                f.writelines(report)
                        
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
                            if self._onlyOne:
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

    def findAllMaximumByOneIteration(self, starting_lengths):
        fileName = './res/'+self._fileNamePref+'.txt'
        self.printLog('Applying all possible subgraphs parametrizations to')
        starting_graph = GraphEmbedding(starting_lengths, self._graph_type, window=self._window, tmpFileName=self._fileNamePref)
        sols = starting_graph.findEmbeddings()
        reached_max = len(sols['real'])       
        self.printLog(starting_lengths)
        self.printLog('with '+str(reached_max)+ ' embeddings\n')
        self.printLog('Working file:')
        self.printLog(fileName)
        
        for uvwpc in self._combinations:
            starting_graph.setLengths(starting_lengths)
            [argmax, maximum] = self.computeSamplingPhiTheta(starting_graph, uvwpc)
            if maximum>reached_max:
                reached_max = maximum
                with open(fileName, 'w') as f:
                    f.write('')
                    
            if maximum==reached_max:
                with open(fileName, 'a') as f:
                    for phi, theta in argmax:
                        starting_graph.setPhiTheta(uvwpc, phi, theta)
                        f.write(str(starting_graph.getLengths())+'\n')
        
        os.rename(fileName, './res/generated_'+str(reached_max)+'_embd_'+self._fileNamePref+'.txt')
        self.printLog('Result saved to:')
        self.printLog('./res/generated_'+str(reached_max)+'_embd_'+self._fileNamePref+'.txt')

    def findMoreEmbeddings(self, starting_lengths, required_num=None):
        if required_num==None:
            required_num = self._numAllSol
        
        G = GraphEmbedding(starting_lengths, self._graph_type, tmpFileName=self._fileNamePref)
        sols = G.findEmbeddings()
        actMaximum = len(sols['real'])
        fromStr = 'from_'+str(actMaximum)
        
        
        self.printLog('Finding more embeddings - linear')
        self.printLog('File name: '+self._fileNamePref+'\nStarting lengths:')
        self.printLog(starting_lengths)
        self.printLog(str(actMaximum)+ ' embeddings\n\n')
        
        N = len(self._combinations)
       
        lens_to_check = [[[copy.copy(starting_lengths)], 0, actMaximum]]
        maxSaved = actMaximum

        found = False
        while lens_to_check and not found:
            lengths_seq, comb_counter, embds = lens_to_check.pop()
            lengths = copy.copy(lengths_seq[-1])
            num_not_changed = 0
            prev_max = embds
            while num_not_changed<N and not found:
                uvwpc = self._combinations[comb_counter % N]
                self.printLog('\nWorking file name: '+self._fileNamePref+'_intermediateResults.txt')
                self.printLog('Reached maximum: '+str(actMaximum))
                self.printLog('Iteration: '+str(comb_counter))
                self.printLog('Actual step: ' + str(uvwpc))
                [_, _, res_lengths, _, maximum] = self.sampleToGetMoreEmbd(lengths, uvwpc)
                lengths = copy.copy(res_lengths[0])
                lengths_seq.append(lengths)
                
                if maximum>actMaximum:
                    actMaximum = maximum
                    self.printLog('MAXIMUM INCREASED TO %d' % maximum)
                    with open('./res/'+self._fileNamePref+'_intermediateResults.txt', 'w') as f:
                        f.write('')
                        
                if maximum==actMaximum:
                    with open('./res/'+self._fileNamePref+'_intermediateResults.txt', 'a') as f:
                        report  = [
                                str(maximum)+'\n', 
                                str(comb_counter)+'iterations \n'
                                 ]
                        f.writelines(report)
                        f.writelines([str(L)+'\n' for L in lengths_seq])
                        f.write('\n')
                        
                if maximum == required_num:
                    report  = [
                                str(required_num)+' EMBEDDINGS FOUND:', 
                                'Iterations:',
                                comb_counter,
                                'Lengths:',
                                lengths
                                 ]
                    for r in report:
                        self.printLog(r)                                                
                    if self._window:
                        self._window.showDialog(report)
                    found = True
                
                if prev_max<maximum:
                    num_not_changed = 0
                else:
                    num_not_changed += 1
                prev_max = maximum
                
                comb_counter += 1
#                for lens in res_lengths[1:]:
#                    addLengths(maximum, lengths_seq[:-1], comb_counter)
                for lens in res_lengths[1:]:
                    if maximum>maxSaved:
                        lens_to_check = []
                        maxSaved = maximum
                    if maximum==maxSaved:
                        lens_to_check.append([copy.copy(lengths_seq[:-1]+[lens]), comb_counter, maximum]) 
                self.printLog('Number of lengths in stack: %d' % len(lens_to_check))

        
        outputFilename = './res/'+str(actMaximum)+'_embd_'+fromStr+'_'+self._fileNamePref+'.txt'
        os.rename('./res/'+self._fileNamePref+'_intermediateResults.txt', outputFilename)
        
        self.printLog('Reached maximum is in:')
        self.printLog(outputFilename)



    def printLog(self, s, verbose=0, newLine=True):
        if verbose<=self.verbose:
            if self._window:
                self._window.printLog(str(s), verbose=verbose, newLine=newLine)
            else:
                if newLine:
                    print s
                else:
                    print s, 
                sys.stdout.flush()