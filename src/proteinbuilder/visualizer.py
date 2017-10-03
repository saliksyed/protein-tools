from OpenGL.GLUT import *
from OpenGL.GLU import *
from OpenGL.GL import *
import pyrr
import numpy as np
import math
import traceback

class Visualizer:
    def __init__(self, width, height, chain):
        glutInit(sys.argv)
        glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH)
        self.width = width
        self.height = height
        glutInitWindowSize(width, height)
        glutCreateWindow('Visualization')
        glClearColor(0.,0.,0.,0.)
        glClampColor(GL_CLAMP_READ_COLOR, GL_FALSE)
        glEnable(GL_CULL_FACE)
        glEnable(GL_DEPTH_TEST)
        glEnable( GL_PROGRAM_POINT_SIZE )
        glMatrixMode(GL_MODELVIEW)
        glutIdleFunc(self.idle)
        glutDisplayFunc(self.display)
        self.chain = chain
        self.t = 0
        glutMainLoop()

    def idle(self):
        glutPostRedisplay()

    def get_bounds(self, chain):
        points = []
        for i, atom in enumerate(chain.atoms):
            if atom:
                points.append(atom.get_transformed_position())
        return pyrr.aabb.create_from_points(np.array(points))

    def render_chain(self, chain):
        # TODO: rendering code needs to be cleaner and less hacky
        # render each of the atoms
        while True:
            for i, atom in enumerate(chain.atoms):
                if atom != None:
                    r = 0.1
                    pos = atom.get_transformed_position()
                    glPushMatrix()
                    glTranslatef(pos[0], pos[1], pos[2])
                    if i == chain.child_atom_idx:
                        glColor3f(1.0, 0.0, 0.0)
                        r = 0.2
                        glutSolidSphere(r, 20, 20)
                    elif i == chain.parent_atom_idx:
                        glColor3f(0.0, 0.0, 1.0)
                        r = 0.2
                        glutSolidSphere(r, 20, 20)
                    else:
                        glColor3f(0.5, 0.5, 0.5)
                        glutSolidSphere(0.1, 20, 20)
                    glPopMatrix()

            glLineWidth(3.0)
            glBegin(GL_LINES)
            for bond in chain.bonds:
                a1 = bond.get_atom_1()
                a2 = bond.get_atom_2()
                if a1 == None or a2 == None:
                    continue
                pos1 = a1.get_transformed_position()
                pos2 = a2.get_transformed_position()
                glColor3f(chain.color[0], chain.color[1], chain.color[2])
                glVertex3f(pos1[0], pos1[1], pos1[2])
                glVertex3f(pos2[0], pos2[1], pos2[2])
            glEnd()
            # render the child
            if chain.child_peptide:
                glLineWidth(4.0)
                glBegin(GL_LINES)
                glColor3f(1.0, 0.0, 1.0)
                loc_from = chain.atoms[chain.parent_atom_idx].get_transformed_position()
                loc_to = chain.child_peptide.atoms[chain.child_peptide.child_atom_idx].get_transformed_position()
                glVertex3f(loc_from[0], loc_from[1], loc_from[2])
                glVertex3f(loc_to[0], loc_to[1], loc_to[2])
                glEnd()
                chain = chain.child_peptide
            else:
                break


    def display(self):
        glViewport(0, 0, int(self.width), int(self.height))
        self.draw()
        glutSwapBuffers()

    def draw(self):
        try:
            self.t += 0.1
            glClearColor(0.,0.,0.,1.)
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
            glShadeModel(GL_SMOOTH)
            glEnable(GL_CULL_FACE)
            glEnable(GL_DEPTH_TEST)
            glPushMatrix()
            glMatrixMode(GL_PROJECTION)
            glLoadIdentity()
            gluPerspective(45.0, 1.5, 0.1, 1000.0)
            glMatrixMode(GL_MODELVIEW)
            glLoadIdentity()
            bounds = self.get_bounds(self.chain)
            center = 0.5 * (pyrr.aabb.minimum(bounds) + pyrr.aabb.maximum(bounds))
            gluLookAt(center[0], center[1], center[2]-50.0, center[0], center[1], center[2], 0, 1., 0)
            glPushMatrix()
            glRotatef(self.t, 0.0, 1.0, 0.0)
            self.render_chain(self.chain)
            glPopMatrix()
            glPopMatrix()
        except:
            traceback.print_exc()
            exit()


