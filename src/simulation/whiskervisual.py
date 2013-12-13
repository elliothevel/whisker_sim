import trep
from trep import tx, ty, tz, rx, ry, rz
import numpy as np
import os
import pickle
import trep_collisions as tc
from OpenGL.GL import *
from OpenGL.GLUT import *
from OpenGL.GLU import *

DEFAULT_CAMERA_POS = np.array([0.0205725, 0.06025617, 0.04004043])
DEFAULT_CAMERA_ANG = np.array([4.49659265, 0.59, 0.])

class WhiskerVisual(trep.visual.VisualItem3D):
    """Class for drawing rat's head and whisker."""
    def __init__(self, *args, **kwds):
        super(WhiskerVisual, self).__init__(*args, **kwds)
        dirname, filename = os.path.split(os.path.abspath(__file__))
        self.attachDrawing('Head', trep.visual.stlmodel(dirname+'/rathead.stl').draw)
        self.attachDrawing(None, self.draw_constraints)
        self.attachDrawing(None, self.draw_whisker) 

    def draw_constraints(self):
        for x in self.system.constraints:
            x.opengl_draw()

    def draw_whisker(self):
        quad = gluNewQuadric()
        gluQuadricNormals(quad, GLU_SMOOTH)
        mass_list = self.system.masses

        draw_list = [self.system.get_frame('Base Point')]
        for j in range(self.system.num_links):
            draw_list.append(self.system.get_frame('Link-%d' %j))

        for j in range(len(draw_list)-1):    
            frame1 = draw_list[j]
            frame2 = draw_list[j+1]
            frame1_g = frame1.g()
            frame2_g = frame2.g()
            
            glPushAttrib(GL_LIGHTING_BIT )
            glDisable(GL_LIGHTING)
            glBegin(GL_LINES)
    
            glVertex3d(frame1_g[0][3], frame1_g[1][3], frame1_g[2][3])
            glVertex3d(frame2_g[0][3], frame2_g[1][3], frame2_g[2][3]) 

            glEnd()
            glPopAttrib()
            
        gluDeleteQuadric(quad)  

class SimpleWhisker3D(trep.System):
    """
    This defines a series of point masses that lie along a whisker. It is
    used to visualize a full whisker away since it can represent any whisker
    system just by inputing the correct system states. Moreover, only one
    system needs to be created to visualize the entire array.
    """
    def __init__(self, num_links=None, masses=None):
        trep.System.__init__(self)

        # Note that by the time this system is created, all of the points'
        # trajectories have been calculated, so there is no need to add in the
        # damping, stiffness, etc. This greatly simplifies the system that
        # actually needs to be visualized.
        if masses is None: masses = [0.0 for j in range(num_links)]
        self.import_frames(self.point_mass_frames(num_links, masses))
        self.num_links = num_links

    def point_mass_frames(self, num_links, masses):
        frames = [tx('x_base'), [ty('y_base'), [tz('z_base', name='Base Point')]] ]
        for j in range(num_links):
            frames += [tx('x-%d' %j), [ty('y-%d' %j), [tz('z-%d' %j,
                mass=masses[j], name='Link-%d' %j)]] ]
        frames += [rz(0.0, name='Head')]    
        return frames

def convert_trajectory_for_simple_system(whisker, q):
    """ 
    Converts the trajectory of a real whisker to one that can be represented
    with the simple whisker system. 
    """
    # Can probably hard code some of these indices.

    simple_whisker = SimpleWhisker3D(whisker.num_links)        
    Q = []
    #head_index   = simple_sys.get_config('head_angle').index
    x_base_index = simple_whisker.get_config('x_base').index
    y_base_index = simple_whisker.get_config('y_base').index
    z_base_index = simple_whisker.get_config('z_base').index
    for q_point in q:
        q_i = np.zeros(simple_whisker.nQ)
        # Set the system state to the current point in the trajectory.
        whisker.set_state(q=q_point)
        
        # Add base position. This is a constant for every time step. (Not if
        # head is rotating!)
        #q_i[head_index] = real_sys.get_config('head_angle').q
        pos = whisker.get_frame('Base Point').p()
        q_i[x_base_index] = pos[0]
        q_i[y_base_index] = pos[1]
        q_i[z_base_index] = pos[2]

        # Add the positions of each node of the whisker.
        for j in range(whisker.num_links):
            pos = whisker.get_frame('Link-%d' %j).p()
            q_i[simple_whisker.get_config('x-%d' %j).index] = pos[0]
            q_i[simple_whisker.get_config('y-%d' %j).index] = pos[1]
            q_i[simple_whisker.get_config('z-%d' %j).index] = pos[2]

        # Save the current trajectory point to Q.         
        Q.append(q_i)

    return Q  

def visualize_array(surface=None, simulation_data=None, time_scale=1.0):
    """ 
    Draws the full whisker array for a simulation. The argument is a
    dictionary containing the simple configurations for each whisker over time.
    These trajectories are simply the locations of each mass of the whisker.
    """
    if simulation_data is None:
        try:
            f = open(os.getcwd()+'/outputdata/simple_data.p', 'r')
            simulation_data = pickle.load(f)
            f.close()
        except:
            raise Exception("Could not find saved data.")    

    simple_whisker = SimpleWhisker3D(simulation_data['N'])

    if surface:
        if surface['type'] == 'sphere':
            (X, Y, Z) = surface['center']
            R = surface['R']
            surf = tc.surfaces.Sphere(simple_whisker, X, Y, Z,
                    R,frame=simple_whisker.world_frame)
        else:
            pass
        surf.activate_constraint()
        #add_constraint_to_system(simple_whisker.get_frame('World'))

    t = simulation_data['T']
    visual_list = []
    for name in simulation_data.keys():
        if (name == 'T') or (name == 'N'):
            pass
        else:
            q = simulation_data[name]['q']
            if len(q) == len(t):
                visual_list.append(WhiskerVisual(simple_whisker, np.array(t)/time_scale, q))        

    trep.visual.visualize_3d(visual_list, camera_pos=DEFAULT_CAMERA_POS,
            camera_ang=DEFAULT_CAMERA_ANG)        
