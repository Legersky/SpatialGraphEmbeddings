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




import os

if not os.path.isdir(os.path.dirname(os.path.realpath(__file__))+'/../tmp'):
    os.makedirs(os.path.dirname(os.path.realpath(__file__))+'/../tmp')
if not os.path.isdir(os.path.dirname(os.path.realpath(__file__))+'/../results'):
    os.makedirs(os.path.dirname(os.path.realpath(__file__))+'/../results')

import sys

sys.path.insert(0,"/scratch/jlegersk/lib/python2.7/site-packages")
sys.path.insert(0,"/scratch/jlegersk/lib/python2.7/site-packages/sklearn")
