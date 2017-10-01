from OpenGL.GLUT import *
from OpenGL.GLU import *
from OpenGL.GL import *
import sys
import math
import traceback
from topology import *
from pyshader import Editor

class Viz(Editor):
    def __init__(self, resolution, chain):
        Editor.__init__(self, resolution)
        self.chain = chain
        self.shader('render', 'test.frag', 'project.vert')

    def draw(self):
        try:
            glClearColor(0.,0.,0.,1.)
            glShadeModel(GL_SMOOTH)
            glEnable(GL_CULL_FACE)
            glEnable(GL_DEPTH_TEST)
            glPushMatrix()
            self.setup_scene()
            color = [1.0,0.,0.,1.]
            glutSolidSphere(10.0, 20, 20)
            glMaterialfv(GL_FRONT,GL_DIFFUSE,color)
            self.shader('render').draw()
            self.chain.render()
            glPopMatrix()
        except:
            traceback.print_exc()
            raw_input("?")


f = ForceField()
villin = """MTKLSAQVKGSLNITTPG"""
c = f.create_chain(villin) # todo HISTIDINE
r = Viz(Viz.RES1080P, c)
r.run()
