#  This file is based on code written by Alvaro Javier Fuentes Suarez

import subprocess
import numpy as np
import time

class Visualization(object):

    def __init__(self):

        self.subdivision_number = 0

        self.meshes = {} # name -> points,facets,color
        self.normals = {} #name -> normals
        self.polylines = {} # name -> points,lens,color
        self.points = {} # name -> points,color,size

    def set_subdivision_number(self,n):
        self.subdivision_number = n

    def add_mesh(self, mesh_points, mesh_facets, color=[0,0,255], name=None):
        if name is None:
            name = "Mesh%d" % len(self.meshes)
        if name in self.meshes.keys():
            n = len(self.meshes[name][0])
            self.meshes[name][0].extend(mesh_points)
            fs = list(map(list,mesh_facets))
            for f in fs:
                for i in xrange(len(f)):
                    f[i] += n
            self.meshes[name][1].extend(fs)
        else:
            self.meshes[name] = (list(mesh_points), list(mesh_facets), list(color) )

    def add_normals(self, normals, name):
        if self.normals.has_key(name):
            self.normals[name].extend(normals)
        else:
            self.normals[name] = list(normals)

    def add_polyline(self, poly_points, color=[255,0,0], name=None):
        if name is None:
            name = "Polyline%d" % len(self.polylines)
        if name in self.polylines.keys():
            self.polylines[name][0].extend(list(poly_points))
            self.polylines[name][1].append(len(poly_points))
        else:
            self.polylines[name] = (list(poly_points), [len(poly_points)], list(color) )

    def add_points(self, points, color=[0,255,0], name=None,  size=0.6):
        if name is None:
            name = "Points%d" % len(self.points)

        if name in self.points.keys():
            self.points[name][0].extend(list(points))
        else:
            self.points[name] = (list(points), list(color) , size)

    def show(self):
        raise NotImplementedError()

class VisualizationAxel(Visualization):
    """docstring for Visualization2"""
    def __init__(self, axel_path,output_file="vis.axl"):
        super(VisualizationAxel, self).__init__()

        self.output_file = output_file

        self.axel_path = axel_path


    def set_mesh_to_subdivide(self,name):
        self.mesh_to_subdivide = name

    def write_mesh(self, f, name):
        ps, fs, color = self.meshes[name]
        f.write('<mesh color="%d %d %d 1" shader="" name="%s" size="0.05">\n' % tuple(color + [name]) )
        f.write('\t<count>%d 0 %d</count>\n' % (len(ps),len(fs)) ) #points edges facets
        f.write('\t<points>\n')
        for p in ps:
            f.write('\t\t%f %f %f\n' % tuple(p))
        f.write('\t</points>\n')

        if self.normals.has_key(name):
            ns = self.normals[name]
            if len(ns) != len(ps):
                print "Error: not the same number of normals and points in the mesh: %" % name
            else:
                f.write('\t<normals>\n')
                for n in ns:
                    f.write('\t\t%f %f %f\n' % tuple(n))
                f.write('\t</normals>\n')

        #f.write('\t<edges></edges>\n')
        f.write('\t<faces>\n')
        for fa in fs:
            f.write('\t\t%d %s\n' % (len(fa),' '.join(map(str,fa))))
        f.write('\t</faces>\n')
        f.write('</mesh>\n')

    def write_polyline(self, f, name):
        ps, ls, color = self.polylines[name]
        f.write('<mesh color="%d %d %d 1" shader="" name="%s" size="0.2">\n' % tuple(color + [name]) )
        f.write('\t<count>%d %d 0</count>\n' % (len(ps),len(ls)) ) #points edges facets
        f.write('\t<points>\n')
        for p in ps:
            f.write('\t\t%f %f %f\n' % tuple(p))
        f.write('\t</points>\n')
        f.write('\t<edges>\n')
        n = 0
        for l in ls:
            f.write('\t\t%d %s\n' % (l,' '.join(map(str,range(n,n+l))) ) )
            n += l
        f.write('\t</edges>\n')
        f.write('</mesh>\n')

    def write_points(self, f, name):
        ps, color, size = self.points[name]
        f.write('<mesh color="%d %d %d 1" shader="" name="%s" size="%f">\n' % tuple(color + [name]+[size]) )
        f.write('\t<count>%d 0 0</count>\n' % len(ps) ) #points edges facets
        f.write('\t<points>\n')
        for p in ps:
            f.write('\t\t%f %f %f\n' % tuple(p))
        f.write('\t</points>\n')
        f.write('</mesh>\n')

    def write_subdiv_off(self, f):
        ps, fs, _ = self.meshes[self.mesh_to_subdivide]
        f.write('OFF\n')
        f.write('%d %d 0\n' % (len(ps),len(fs)))
        for p in ps:
            f.write('%f %f %f\n' % tuple(p))
        for fa in fs:
            f.write('%d %s\n' % (len(fa),' '.join(map(str,fa))))
        f.write('\n')

    def read_subdiv_off(self,f):
        f.readline() #discard the first line

        line = f.readline()
        while line[0]=="#":
            line = f.readline() #discard comment lines

        #read teh number of points and faces
        nump,numf,_ = map(int,line.split(' '))
        line = f.readline()

        ps = []
        fs = []
        while nump: #read the points
            #discard comment and empty lines
            if line[0]=="#" or len(line)<3:
                line = f.readline()
                continue

            ps.append(np.array(map(float,line.replace('\n', '').split(' '))))

            line = f.readline()
            nump -= 1


        while numf: #read the faces
            #discard comment and empty lines
            if line[0]=="#" or len(line)<3:
                line = f.readline()
                continue

            fs.append(map(int,line.replace('\n', '').split(' ')[2:]))

            line = f.readline()
            numf -= 1

        color = self.meshes[self.mesh_to_subdivide][2]
        del self.meshes[self.mesh_to_subdivide]
        self.add_mesh(ps,fs,color,name=self.mesh_to_subdivide)

    def show(self,save_and_show=True):

        if self.subdivision_number and self.mesh_to_subdivide:
            with open('subdiv_in.off','wt') as fin:
                self.write_subdiv_off(fin)
            with open('subdiv_in.off','rt') as fin, open('subdiv_out.off','wt') as fout:
                subprocess.call([self.subdiv_path, '%d' % self.subdivision_number],stdin=fin,stdout=fout)
                fout.flush()
            with open('subdiv_out.off','rt') as fout:
                self.read_subdiv_off(fout)

            import os
            os.remove('subdiv_out.off')
            os.remove('subdiv_in.off')

        with open(self.output_file,'wt') as f:
            f.write('<axl>\n')
            for name in self.meshes.keys():
                self.write_mesh(f,name)
            for name in self.polylines.keys():
                self.write_polyline(f,name)
            for name in self.points.keys():
                self.write_points(f,name)
            f.write('</axl>\n')
        if save_and_show:
            #subprocess.call([self.axel_path, 'vis.axl'])
            subprocess.Popen([self.axel_path, self.output_file])
            time.sleep(3)

    def save(self):
        self.show(False)


