import trep
import numpy as np
from surface import Surface

try:
    from OpenGL.GL import *
    from OpenGL.GLU import *
    from OpenGL.GLUT import *
    _opengl = True
except:
    _opengl = False


def find_closest_frame(whisker, Y, Z):
    dist = []
    for i in range(whisker.num_links-1):
        frame = whisker.get_frame('Link-%d' %i)
        if frame.p()[1] > Y:
            return frame

class Peg(Surface):
    def __init__(self, system, Y, Z, frame, tolerance=1e-10, name=None):
        Surface.__init__(self, system, name, tolerance)
        self.Y = Y
        self.Z = Z
        self.frame = frame
        self.sign = -1.0
        self.pegrw = np.array((0,Y,Z,1))
    
    def phi(self):
        closest_frame = find_closest_frame(self.system, self.Y, self.Z)
        if self.frame != closest_frame:
            return 1.0
        return np.dot(closest_frame.g_inv(), self.pegrw)[2]

    def phi_dq(self, q1):
        closest_frame = find_closest_frame(self.system, self.Y, self.Z)
        return np.dot(closest_frame.g_inv_dq(q1), self.pegrw)[2]
    
    if _opengl:
        def opengl_draw(self):
            glPushAttrib(GL_CURRENT_BIT | GL_LIGHTING_BIT | GL_LINE_BIT)
            glColor3f(1.0, 0.0, 0.0)
            glDisable(GL_LIGHTING)
            
            glPointSize(5.0)
            glBegin(GL_POINTS)
            glVertex3f(0,self.Y,self.Z)
            glEnd()
                        
            glPopAttrib()

