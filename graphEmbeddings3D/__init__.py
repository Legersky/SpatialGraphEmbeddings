import os

if not os.path.isdir(os.path.dirname(os.path.realpath(__file__))+'/../tmp'):
    os.makedirs(os.path.dirname(os.path.realpath(__file__))+'/../tmp')
if not os.path.isdir(os.path.dirname(os.path.realpath(__file__))+'/../results'):
    os.makedirs(os.path.dirname(os.path.realpath(__file__))+'/../results')

import sys
for i in sys.path:
    print i
print '-----------------'

sys.path.insert(0,"/scratch/jlegersk/lib/python2.7/site-packages/sklearn/")
sys.path.insert(0,"/scratch/jlegersk/lib/python2.7/site-packages/sklearn/")
for i in sys.path:
    print i
print '-----------------'
