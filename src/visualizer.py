from OpenGL.GLUT import *
from OpenGL.GLU import *
from OpenGL.GL import *
from pyshader import Editor


class Visualizer(Editor):
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