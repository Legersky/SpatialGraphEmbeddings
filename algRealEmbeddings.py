import time
#import copy
from sklearn.cluster import DBSCAN

class AlgRealEmbeddings(object):
    def __init__(self, graph, window=None):
        self._window = window

    def getSamplingIterator(self, lengths, num_phi, num_theta,  l_phi,  r_phi, l_theta,  r_theta):
        step_phi = (r_phi-l_phi)/float(num_phi)
        step_theta = (r_theta-l_theta)/float(num_theta)
        phi = l_phi
        while (phi<r_phi+step_phi):
            theta = l_theta
            while theta<r_theta+step_theta:
                yield self.getLengthsForPhiTheta(lengths, phi, theta),  phi, theta
                theta += step_theta
            phi += step_phi

    def runSamplingPhiThetaWithMargins(self, lengths, num_phi, num_theta,  l_phi,  r_phi, l_theta,  r_theta, treshold=0):
        return self.runSampling(self.getSamplingIterator(lengths, num_phi, num_theta,  l_phi,  r_phi, l_theta,  r_theta), treshold=treshold)

    def runSampling(self, iterator, treshold=0, animate=False):
        solutions = []
        for lengths, phi, theta in iterator:
            num_real = len(self.findEmbeddings(lengths, errorMsg=False)['real'])
            
            self.printLog(str([phi, theta, num_real]), verbose=1)
            
#            self.setLengthsAndUpdateFixedTriangle(lengths, setRequiresRecomputing=False)
#            self._window.update_graph2phi()
#            self._window.update_graph2theta()
#            self.computeIntersections()
#            self._window.plotScene()
            if num_real>=treshold:
                solutions.append([phi, theta, num_real])              
        return solutions

    def computeSamplingPhiTheta(self, lengths, num_phi, num_theta):
        start = time.time()
        act_num = len(self.findEmbeddings(lengths)['real'])
        act_phi, act_theta = self.getPhiTheta(lengths)
       
        margin_degree = 5
        margin = margin_degree*math.pi/180.0
        l_phi,  r_phi = -math.pi/2.0 + margin, math.pi/2.0 - margin
        l_theta, r_theta = 0 + margin/2.0, math.pi - margin/2.0
        
        sols = self.runSamplingPhiThetaWithMargins(lengths, num_phi, num_theta,  l_phi,  r_phi, l_theta,  r_theta,  treshold=act_num)
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
            for phi2, theta2, num2 in self.runSamplingPhiThetaWithMargins(lengths, num_phi/4, num_theta/4,
                                                                          phi - step_phi,  phi + step_phi,
                                                                          theta - step_theta,  theta + step_theta,
                                                                          treshold=maximum):
                if num2>maximum:
                    self.printLog(str(num2)+',trying increase max', verbose=2)
                    num_check = len(self.findEmbeddings(self.getLengthsForPhiTheta(lengths, phi2, theta2), usePrev=False)['real'])
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
            n = len(self.findEmbeddings(self.getLengthsForPhiTheta(lengths, phi_c, theta_c))['real'])
            if n < maximum:
                min_dist = dist_angle([phi_c, theta_c], cluster[0])
                phi_c = cluster[0][0]
                theta_c = cluster[0][1]
                for x, y in cluster:
                    d = dist_angle([phi_c, theta_c], [x, y])
                    if d < min_dist:
                        min_dist = d
                        phi_c = x
                        theta_c = y
                self.printLog('Center of cluster does not have maximum number of solutions \n -> nearest point chosen instead.')
            
            centers.append([phi_c, theta_c])
            
            res_lengths.append(self.getLengthsForPhiTheta(lengths, phi_c, theta_c))
            res_infos.append(str([self.getPhiRadian(), self.getThetaRadian()]))
        
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
        
        return [clusters, centers, res_lengths, res_infos]
        
        

    def runSamplingPhiTheta(self, num_phi, num_theta):
        [clusters, centers, res_lengths, res_infos] = self.computeSamplingPhiTheta(self.getLengths(), num_phi, num_theta)        
        if self._window:
            self._window.showClusters(clusters, centers)
            self._window.setGraphSequence([GraphCouplerCurve(lengths, window=self._window) for lengths in res_lengths], res_infos)
