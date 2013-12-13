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


class Plane(Surface):
    def __init__(self, system, point, normal, frame, tolerance=1e-10, name=None):
        Surface.__init__(self, system, name, tolerance)
        self.normal = normal
        self.point = point
        self.frame = frame
        self.plane_d = -np.dot(point, normal)
        self.plane_r = np.sqrt(np.dot(normal,normal))
        self.sign = -1.0

    def phi(self):
        p = self.frame.p()[:-1]
        return (np.dot(p,self.normal)+self.plane_d)/self.plane_r

    def phi_dq(self, q1):
        x, y, z = self.frame.p()[:-1] 
        dq = self.frame.p_dq(q1)[:-1]
        return np.dot(p,self.normal)/self.plane_r 

    def __repr__(self):
            return "<Plane surface frame=%s>" %self.frame

    if _opengl:
        def opengl_draw(self):
            """ Draw a simple sphere. """
            glPushAttrib(GL_CURRENT_BIT | GL_LIGHTING_BIT | GL_LINE_BIT)
            glColor3f(96/255.0, 96/255.0, 96/255.0)
            glDisable(GL_LIGHTING)
            glLineWidth(2.0)

            glBegin(GL_QUADS)
            glVertex4f(x_min, y_min, 0.0)
            glVertex4f(x_min, y_max, 0.0)
            glVertex4f(x_max, y_max, 0.0)
            glVertex4f(x_max, y_min, 0.0)
            glEnd()

            glPopAttrib() 
