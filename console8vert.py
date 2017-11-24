from algRealEmbeddings import *
import time
from random import uniform
import copy


#v1 = [uniform(a, b), uniform(a, b), uniform(a, b)]
#v7 = [uniform(a, b), uniform(a, b), uniform(a, b)]
#
#v2 = [1.0+uniform(a, b), uniform(a, b), uniform(a, b)]
#v5 = [1.0+uniform(a, b), uniform(a, b), uniform(a, b)]
#
#v3 = [uniform(a, b), 1.0+uniform(a, b), uniform(a, b)]
#
#v6 = [uniform(a, b), uniform(a, b), 1.0+uniform(a, b)]
#
#v4 = [1.0+uniform(a, b), 1.0+uniform(a, b), uniform(a, b)]
#
#v8 = [ 1.0 + uniform(a, b),uniform(a, b), 1.0+uniform(a, b)]

a = -5.0
b = 5.0
v1 = [uniform(a, b), uniform(a, b), uniform(a, b)]
v7 = [uniform(a, b), uniform(a, b), uniform(a, b)]

v2 = [uniform(a, b), uniform(a, b), uniform(a, b)]
v5 = [uniform(a, b), uniform(a, b), uniform(a, b)]

v3 = [uniform(a, b), uniform(a, b), uniform(a, b)]

v6 = [uniform(a, b), uniform(a, b), uniform(a, b)]

v4 = [uniform(a, b), uniform(a, b), uniform(a, b)]

v8 = [uniform(a, b), uniform(a, b), uniform(a, b)]

#v1 = [uniform(a, b), 1.0 + uniform(a, b), uniform(a, b)]
#v8 = [uniform(a, b), 1.0 + uniform(a, b), uniform(a, b)]

#v4 = [-1.0+uniform(a, b), uniform(a, b), uniform(a, b)]
#v2 = [-1.0+uniform(a, b), uniform(a, b), uniform(a, b)]

#v3 = [1.0+uniform(a, b), uniform(a, b), uniform(a, b)]
#v5 = [1.0+uniform(a, b), uniform(a, b), uniform(a, b)]

#v6 = [uniform(a, b), 0.5 + uniform(a, b), 1.0+uniform(a, b)]

#v7 = [uniform(a, b), 0.5 + uniform(a, b), -1.0+uniform(a, b)]

def dist( u, v):
    return float(np.sqrt( (u[0]-v[0])**2 + (u[1]-v[1])**2 + (u[2]-v[2])**2))
#
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
            (2, 8) : dist(v2,v8),
           }


#lengths = {(2, 7): 1.4288466263977495, (4, 7): 17.88035095820198, (2, 6): 1.0237800945583875, (4, 5): 17.84874384133012, (2, 8): 1.0222718473388588,
#           (1, 4): 18.149011451612797, (1, 3): 14.357024144718817, (1, 6): 3.032569805432945, (3, 7): 14.26976782945725, (5, 8): 0.09782324287506733,
#           (1, 2): 3.214648729657953, (6, 8): 0.03813348403581635, (1, 5): 3.032338165266229, (3, 4): 22.24815495162433, (5, 7): 1.0123496857434562, 
#           (5, 6): 0.09495461781471154, (2, 3): 13.512023838396058, (7, 8): 1.0091111235946963}
#           
#lengths = {(2, 7): 1.4288466263977495, (4, 7): 13.177186219752265, (2, 6): 1.0237800945584088, (4, 5): 13.095230370891187, (2, 8): 1.0222587379875785, 
#           (1, 4): 13.365253277214508, (1, 3): 7.715718492411297, (1, 6): 3.032569805432938, (3, 7): 7.296845203770255, (5, 8): 0.09958757549720172, 
#           (1, 2): 3.214648729657953, (6, 8): 0.03293162305363898, (1, 5): 3.032338165266229, (3, 4): 14.352417603885499, (5, 7): 1.0123496857434562, 
#           (5, 6): 0.09495461781465638, (2, 3): 6.567028630020083, (7, 8): 0.9888906344533466}
#           
#lengths = {(2, 7): 1.4288466263977495, (4, 7): 4.842624306865316, (2, 6): 1.0238654856293545, (4, 5): 4.5546452094823247, (2, 8): 1.022390122314935, 
#           (1, 4): 5.501616644087601, (1, 3): 8.51207449423649, (1, 6): 3.0275061849268776, (3, 7): 8.16277576913792, (5, 8): 0.09814257207967064, 
#           (1, 2): 3.214648729657953, (6, 8): 0.033081652556535536, (1, 5): 3.032338165266229, (3, 4): 7.634708567878329, (5, 7): 1.0123496857434562, 
#           (5, 6): 0.09499799392440313, (2, 3): 7.420903583333997, (7, 8): 0.9940938333563378}
           
           
#lengths = {(2, 7): 1.4288466263977495, (4, 7): 3.646175820991137, (2, 6): 1.0238654856293283, (4, 5): 3.3429830998046732, (2, 8): 1.0229346087847153, 
#           (1, 4): 4.49321454664234, (1, 3): 12.291475166368201, (1, 6): 3.027506184926906, (3, 7): 12.213160380343274, (5, 8): 0.09529692624044482, 
#           (1, 2): 3.214648729657953, (6, 8): 0.0179005231972465, (1, 5): 3.032338165266229, (3, 4): 10.92113783226037, (5, 7): 1.0123496857434562, 
#           (5, 6): 0.0949979939254927, (2, 3): 11.53333847125022, (7, 8): 1.0019027335371273}


start = time.time()
print 'Magnitude system:'
alg = AlgRealEmbeddings('Max8vertices', name='test_magn_syst')
alg.sampleToGetMoreEmbd(lengths, [8, 2, 7, 6, 5],  0)
end = time.time()
print 'Final time: '+str(end - start)


start = time.time()
print 'Distance system:'
alg = AlgRealEmbeddings('Max8vertices_distSyst', name='test_dist_syst')
alg.sampleToGetMoreEmbd(lengths, [8, 2, 7, 6, 5],  0)
end = time.time()
print 'Final time: '+str(end - start)



#start = time.time()
#print '\nDistance system'
#G = GraphEmbedding(lengths, 'Max8vertices_distSyst')
#sols = G.findEmbeddings()
#print '# embeddable (dist syst):'
#print len(sols['real'])
#end = time.time()
#print 'time: '+str(end - start)

#if m==2*len(sols['real']):
#    print 'OK'
#else:
#    print 'PROBLEM!!!'

#print len(sols['complex'])
#lengths = {(1, 2): 3.2232294709034446, (4, 7): 3.8174034354261854, (2, 6): 1.023865485629374, (4, 5): 3.4734613814676014, (2, 8): 1.0229346087848081,
#           (1, 4): 4.5144979036144388, (1, 3): 18.211680156322096, (1, 6): 3.0275061849268874, (3, 7): 18.412602590272094, (5, 8): 0.09529692623985034,
#           (2, 7): 1.4047534656392475, (6, 8): 0.018779861626871987, (1, 5): 3.032338165266229, (7, 8): 1.0066546644867307, (5, 7): 1.0123496857434562,
#           (5, 6): 0.0949979939274194, (2, 3): 17.74254212711772, (3, 4): 16.70963254790826}
#   
#
#L = [{(1, 2): 3.2232294709034446, (4, 7): 3.8174034354261854, (2, 6): 1.023865485629374, (4, 5): 3.4734613814676014, (2, 8): 1.0229346087848081, (1, 4): 4.514497903614439, (1, 3): 11.960918976585539, (1, 6): 3.0275061849268874, (3, 7): 12.043848163688027, (5, 8): 0.09529692623985034, (2, 7): 1.4047534656392475, (6, 8): 0.018779861626871987, (1, 5): 3.032338165266229, (7, 8): 1.0066546644867307, (5, 7): 1.0123496857434562, (5, 6): 0.0949979939274194, (2, 3): 11.396942840889189, (3, 4): 10.534284999704216}, 
#    {(1, 2): 3.2232294709034446, (4, 7): 3.8174034354261854, (2, 6): 1.023865485629374, (4, 5): 3.4734613814676014, (2, 8): 1.0229346087848081, (1, 4): 4.514497903614439, (1, 3): 11.448361258356769, (1, 6): 3.0275061849268874, (3, 7): 11.523580083328511, (5, 8): 0.09529692623985034, (2, 7): 1.4047534656392475, (6, 8): 0.018779861626871987, (1, 5): 3.032338165266229, (7, 8): 1.0066546644867307, (5, 7): 1.0123496857434562, (5, 6): 0.0949979939274194, (2, 3): 10.879761937089897, (3, 4): 10.040497985715902}, 
#    {(1, 2): 3.2232294709034446, (4, 7): 4.160009183178643, (2, 6): 1.023865485629374, (4, 5): 3.827147328563408, (2, 8): 1.0229346087848081, (1, 4): 4.782601524915176, (1, 3): 18.211680156322096, (1, 6): 3.0275061849268874, (3, 7): 18.412602590272094, (5, 8): 0.09529692623985034, (2, 7): 1.4047534656392475, (6, 8): 0.018779861626871987, (1, 5): 3.032338165266229, (7, 8): 1.0066546644867307, (5, 7): 1.0123496857434562, (5, 6): 0.0949979939274194, (2, 3): 17.74254212711772, (3, 4): 16.748943042533085}, 
#    {(1, 2): 3.2232294709034446, (4, 7): 3.9218264570026635, (2, 6): 1.023865485629374, (4, 5): 3.581476938942793, (2, 8): 1.0229346087848081, (1, 4): 4.595148522944032, (1, 3): 18.211680156322096, (1, 6): 3.0275061849268874, (3, 7): 18.412602590272094, (5, 8): 0.09529692623985034, (2, 7): 1.4047534656392475, (6, 8): 0.018779861626871987, (1, 5): 3.032338165266229, (7, 8): 1.0066546644867307, (5, 7): 1.0123496857434562, (5, 6): 0.0949979939274194, (2, 3): 17.74254212711772, (3, 4): 16.790694351414444}, 
#    {(1, 2): 3.2232294709034446, (4, 7): 3.8174034354261854, (2, 6): 1.023865485629374, (4, 5): 3.4734613814676014, (2, 8): 1.0229346087848081, (1, 4): 4.514497903614439, (1, 3): 13.927210641812112, (1, 6): 3.0275061849268874, (3, 7): 14.051203234509808, (5, 8): 0.09529692623985034, (2, 7): 1.4047534656392475, (6, 8): 0.018779861626871987, (1, 5): 3.032338165266229, (7, 8): 1.0066546644867307, (5, 7): 1.0123496857434562, (5, 6): 0.0949979939274194, (2, 3): 13.36776176442962, (3, 4): 12.591946342294957}, 
#    {(1, 2): 3.2232294709034446, (4, 7): 3.8174034354261854, (2, 6): 1.023865485629374, (4, 5): 3.4734613814676014, (2, 8): 1.0229346087848081, (1, 4): 4.514497903614439, (1, 3): 13.227901143687951, (1, 6): 3.0275061849268874, (3, 7): 13.334410186799714, (5, 8): 0.09529692623985034, (2, 7): 1.4047534656392475, (6, 8): 0.018779861626871987, (1, 5): 3.032338165266229, (7, 8): 1.0066546644867307, (5, 7): 1.0123496857434562, (5, 6): 0.0949979939274194, (2, 3): 12.647927033216305, (3, 4): 11.923017253203193}, 
#    {(1, 2): 3.2232294709034446, (4, 7): 3.8174034354261854, (2, 6): 1.023865485629374, (4, 5): 3.4734613814676014, (2, 8): 1.0229346087848081, (1, 4): 4.514497903614439, (1, 3): 17.130381998221665, (1, 6): 3.0275061849268874, (3, 7): 17.346872128951311, (5, 8): 0.09529692623985034, (2, 7): 1.4047534656392475, (6, 8): 0.018779861626871987, (1, 5): 3.032338165266229, (7, 8): 1.0066546644867307, (5, 7): 1.0123496857434562, (5, 6): 0.0949979939274194, (2, 3): 16.679353422619492, (3, 4): 15.662725189566396}, 
#    {(1, 2): 3.2232294709034446, (4, 7): 3.8174034354261854, (2, 6): 1.023865485629374, (4, 5): 3.4734613814676014, (2, 8): 1.0229346087848081, (1, 4): 4.514497903614439, (1, 3): 34.741811568382815, (1, 6): 3.0275061849268874, (3, 7): 35.036883454355475, (5, 8): 0.09529692623985034, (2, 7): 1.4047534656392475, (6, 8): 0.018779861626871987, (1, 5): 3.032338165266229, (7, 8): 1.0066546644867307, (5, 7): 1.0123496857434562, (5, 6): 0.0949979939274194, (2, 3): 34.320764213216783, (3, 4): 33.50090321431538}, 
#    {(1, 2): 3.2232294709034446, (4, 7): 3.8174034354261854, (2, 6): 1.023865485629374, (4, 5): 3.4734613814676014, (2, 8): 1.0229346087848081, (1, 4): 4.514497903614439, (1, 3): 30.917431458441378, (1, 6): 3.0275061849268874, (3, 7): 31.199937656313235, (5, 8): 0.09529692623985034, (2, 7): 1.4047534656392475, (6, 8): 0.018779861626871987, (1, 5): 3.032338165266229, (7, 8): 1.0066546644867307, (5, 7): 1.0123496857434562, (5, 6): 0.0949979939274194, (2, 3): 30.484198141245734, (3, 4): 29.641400362571787}]
#
#for lengths in L:











