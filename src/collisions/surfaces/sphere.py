import trep
import numpy as np
from surface import Surface

# Import the visualization module. If this fails, the constraint will not
#   be drawn in the simulation.
try:
    from OpenGL.GL import *
    from OpenGL.GLU import *
    from OpenGL.GLUT import *
    _opengl = True
except:
    _opengl = False


class Sphere(Surface):
    def __init__(self, system, X, Y, Z, R, frame, tolerance=1e-10, name=None):
        Surface.__init__(self, system, name, tolerance)
        self.X = X
        self.Y = Y
        self.Z = Z
        self.R = R
        self.frame = frame
        self.sign = -1.0

    def phi(self):
        x, y, z = self.frame.p()[:-1]
        return (x-self.X)**2 + (y-self.Y)**2 + (z-self.Z)**2 - self.R**2

    def phi_dq(self, q1):
        x, y, z = self.frame.p()[:-1] 
        dx, dy, dz = self.frame.p_dq(q1)[:-1]
        return 2*(x-self.X)*dx + 2*(y-self.Y)*dy + 2*(z-self.Z)*dz

    def phi_dqdq(self, q1, q2):
        return 1.0
        x, y, z = self.frame.p()[:-1]
        dx01, dy01, dz01 = self.frame.p_dq(q2)[:-1]
        dx10, dy10, dz10 = self.frame.p_dq(q1)[:-1]
        dx11, dy11, dz11 = self.frame.p_dqdq(q1,q2)[:-1]
        return 2*(dx01*dx10+dy01*dy10+dz01*dz10+
                (x-self.X)*dx11+
                (y-self.Y)*dy11+
                (z-self.Z)*dz11)
        
    def __repr__(self):
            return "<Sphere surface frame=%s>" %self.frame

    if _opengl:
        def opengl_draw(self):
            """ Draw a simple sphere. """
            glPushAttrib(GL_CURRENT_BIT | GL_LIGHTING_BIT | GL_LINE_BIT)
            glColor3f(255/255.0, 0/255.0, 0/255.0)
            subdivisions = 50
            color = [1,0,0,1]
            quad = gluNewQuadric()
            gluQuadricNormals(quad, GLU_SMOOTH)
            glPushMatrix()
            glTranslatef(self.X, self.Y, self.Z)
            glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, color)
            glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, .1)
            gluSphere(quad, self.R, subdivisions, subdivisions)
            glPopMatrix()
            glPopAttrib()
            gluDeleteQuadric(quad)

