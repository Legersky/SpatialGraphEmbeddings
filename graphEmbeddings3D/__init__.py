import os

if not os.path.isdir(os.path.dirname(os.path.realpath(__file__))+'/../tmp'):
    os.makedirs(os.path.dirname(os.path.realpath(__file__))+'/../tmp')
if not os.path.isdir(os.path.dirname(os.path.realpath(__file__))+'/../results'):
    os.makedirs(os.path.dirname(os.path.realpath(__file__))+'/../results')

import sys

sys.path.insert(0,"/scratch/jlegersk/lib/python2.7/site-packages")
sys.path.insert(0,"/scratch/jlegersk/lib/python2.7/site-packages/sklearn")
