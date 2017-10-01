from OpenGL.GLUT import *
from OpenGL.GLU import *
from OpenGL.GL import *
import sys
import math
import traceback
from topology import *

class Renderer:
    RES8K = [7680, 4320]
    RES4K = [3840, 2160]
    RES720P = [1280, 720]
    RES1440P = [2560, 1440]
    RES1080P = [1920, 1080]
    RES1024 = [1024, 768]
    RES2880 = [2880, 1800]
    RES1280 = [1280, 800]
    RES1440 = [1440, 900]
    RES1920 = [1920, 1200]

    def __init__(self, resolution, chain):
        self.chain = chain
        self.eye_pos = [-4.95, 13.65, -10.43]
        self.target_dist = 10.
        self.angle = 0.0
        self.tilt_angle = 0.0
        self.state = None
        self.scale = 1.0
        self.angle = 0.0
        self.width = resolution[0]
        self.height = resolution[1]
        self.init_gl()

    def move(self, key, x, y):
        """
            Handle keyboard controls
        """
        if key == 'w':
            self.scale += 0.1
        elif key == 'a':
            self.angle += 1.0
        elif key == 's':
            self.scale -= 0.1
        elif key == 'd':
            self.angle -= 1.

    def init_gl(self):
        glutInit(sys.argv)
        glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH)
        glutInitWindowSize(self.width, self.height)
        glutCreateWindow('Visualization')
        glClearColor(0.,0.,0.,1.)
        glShadeModel(GL_SMOOTH)
        glEnable(GL_CULL_FACE)
        glEnable(GL_DEPTH_TEST)
        glutDisplayFunc(self.display)
        glutKeyboardFunc(self.move)
        glMatrixMode(GL_PROJECTION)
        gluPerspective(45, 1.5, 0.1, 1000.)
        glMatrixMode(GL_MODELVIEW)
        gluLookAt(0,0,10,
                  0,0,0,
                  0,1,0)
        glPushMatrix()
        glutIdleFunc(self.idle)
        glutMainLoop()

    def display(self):
        glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT)
        glPushMatrix()
        color = [1.0,0.,0.,1.]
        glMaterialfv(GL_FRONT,GL_DIFFUSE,color)
        glPushMatrix()
        glRotatef(self.angle, 0, 1.0, 0)
        glScalef(self.scale, self.scale, self.scale)
        try:
            c.render()
        except:
            traceback.print_exc()
            raw_input("?")
        glPopMatrix()
        glPopMatrix()
        glutSwapBuffers()
    
    def idle(self):
        glutPostRedisplay()


f = ForceField()
c = f.create_chain('ARNDCQEGILKMPFSTWYV') # todo HISTIDINE does not work!
r = Renderer(Renderer.RES1080P, c)
