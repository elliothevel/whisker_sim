import trep
from trep import rx, ry, rz, tx, ty, tz, const_txyz
import numpy as np
from parameters import *


class FixedAngle(trep.Constraint):
    """ Defines a constraint on an angle that should remain at a fixed value. """
    def __init__(self, system, angle_name, reference, name):
        trep.Constraint.__init__(self, system, name=name)
        self.angle = system.get_config(angle_name)
        self.reference = reference

    def h(self):
        return self.angle.q - self.reference

    def h_dq(self, q1):
        if q1 == self.angle:
            return 1.0
        return 0.0

    def h_dqdq(self, q1, q2):
        return 0.0
 
# Define the constraint on the base angle.
class BaseAngle(trep.Constraint):
    """ This defines a constraint on the base motion of the whisker.
        It is used to drive the whisking motion. """
    def __init__(self, system, angle_name, h_func, hdt_func, name=None):
        trep.Constraint.__init__(self, system, name=name)
        self.base_angle = system.get_config(angle_name)
        self.h_func = h_func
        self.hdt_func = hdt_func

    def h(self):
        return self.base_angle.q - self.h_func(self.system.t)

    def h_dq(self, q1):
        if self.base_angle == q1:
            return 1.0
        return 0.0

    def h_dqdq(self, q1, q2):
        return 0.0

    def h_dt(self):
        return self.hdt_func(self.system.t)


class Whisker2D(trep.System):
    """ Defines a single, two dimensional whisker. """
    def __init__(self, parameters, base_motion):
        trep.System.__init__(self)

        self.num_links = parameters['N']
        self.lengths = parameters['L'] 
        self.import_frames(self.whisker_frames(parameters['a']))
        self.add_base_constraints(base_motion)

        self.add_node_springs(parameters['k'])
        self.add_node_damping(parameters['c'])
        self.set_masses(parameters['m'])
        self.satisfy_constraints()

    def whisker_frames(self, curvature):
        """ Creates a list of frames that define the whisker. """
        if curvature is None:
            ref_angles = np.zeros(self.num_links)
        else:
            ref_angles = get_angles_from_curvature(self.lengths, curvature)

        frames = []
        for j in reversed(range(1, self.num_links)):
            frames = [ty(self.lengths[j], name='Link-%d' %j), frames]
            frames = [rx('theta-%d' %j), frames]
            frames = [rx(ref_angles[j]), frames]
        frames = [ty(self.lengths[0], name='Link-0'), frames]
        frames = [tz('z'), [ty('y', name='Base_Point'), frames]]
        frames = [rx('theta-0'), frames]
        frames = [rx(ref_angles[0]), frames]
        return frames

    def set_masses(self, m):
        """ Sets masses at each of the whisker's nodes. """
        for j in range(self.num_links):
            self.get_frame('Link-%d' %j).set_mass((m[j],0,0,0))
        self.get_frame('Base_Point').set_mass((m[0],0,0,0))
            
    def add_node_springs(self, k):
        """ Adds configuration (torsional) springs at each node. """
        for j in range(self.num_links):
            trep.potentials.ConfigSpring(self, 'theta-%d' %j, k=k[j], name='spring-%d' %j)

    def add_node_damping(self, c, default=0.0):
        """ Adds damping at each node. """
        damp_dict = {}
        for j in range(self.num_links):
            damp_dict.update({'theta-%d' %j: c[j]})   
        trep.forces.Damping(self, default, damp_dict)

    def add_base_constraints(self, base_motion):
        """ Adds constraints on the rotation and translation at whisker's base. """
        BaseAngle(self, 'theta-0', base_motion['h'], base_motion['hdt'], name='M')
        trep.constraints.PointOnPlane(self, 'Base_Point', (0,0,1), 'World', name='FZ') 
        trep.constraints.PointOnPlane(self, 'Base_Point', (0,1,0), 'World', name='FY')


class Whisker3D(trep.System):
    """ Defines a single, three dimensional whisker. """
    def __init__(self, parameters, base_motion, name=None):
        trep.System.__init__(self)

        self.name = name
        self.num_links = parameters['N']
        self.lengths = parameters['L'] 
        self.import_frames(self.whisker_frames(parameters['a'],
                           parameters['base_pos'], parameters['base_rot']))
        self.add_base_constraints(base_motion)

        self.add_node_springs(parameters['k'])
        self.add_node_damping(parameters['c'])
        self.set_masses(parameters['m'])
        self.satisfy_constraints()

    def whisker_frames(self, curvature, base_pos, base_rot):
        """ Creates a list of frames that define the whisker. """
        if curvature is None:
            ref_angles = np.zeros(self.num_links)
        else:
            ref_angles = get_angles_from_curvature(self.lengths, curvature)

        frames = []
        for j in reversed(range(1, self.num_links)):
            frames = [tx(self.lengths[j], name='Link-%d' %j), frames]
            frames = [rz('theta-%d_z' %j), frames]
            frames = [ry('theta-%d_y' %j), frames]
            #frames = [rx('theta-%d_x' %j), frames]
            frames = [rz(-ref_angles[j]), frames]
        frames = [tx(self.lengths[0], name='Link-0'), frames]
        frames = [rz('theta-0_z', name='Moving Base Point'), frames]
        frames = [ry('theta-0_y'), frames]
        frames = [rx('theta-0_x'), frames]
        frames = [tz('z'), [ty('y'), [tx('x'), frames]]]

        (X, Y, Z)          = base_pos
        (theta, phi, zeta) = base_rot

        # Rotate to the correct position.
        frames = [rz(theta), [ry(phi), [rx(zeta), frames]]]

        # Place the whisker at the correct spot on the mystacial pad and add an angle
        #   that drives the whisking motion.
        frames = [tx(X), [ty(Y),[ tz(Z, name='Base Point'),
                 [rz('Driving Angle', name="Driving Angle"), frames]]]]

        frames = [ rz(0.0, name='Head'), frames]

        return frames

    def set_masses(self, m):
        """ Sets masses at each of the whisker's nodes. """
        for j in range(self.num_links):
            self.get_frame('Link-%d' %j).set_mass((m[j],0,0,0))
        self.get_frame('Base Point').set_mass((m[0],0,0,0))
            
    def add_node_springs(self, k):
        """ Adds configuration (torsional) springs at each node. """
        for j in range(self.num_links):
            trep.potentials.ConfigSpring(self, 'theta-%d_z' %j, k=k[j],
                    name='spring-%d_z' %j)
            trep.potentials.ConfigSpring(self, 'theta-%d_y' %j, k=k[j],
                    name='spring-%d_y' %j)
            if j == 0:
                trep.potentials.ConfigSpring(self, 'theta-%d_x' %j, k=k[j],
                        name='spring-%d_x' %j)

    def add_node_damping(self, c, default=0.0):
        """ Adds damping at each node. """
        damp_dict = {}
        for j in range(self.num_links):
            damp_dict.update({'theta-%d_y' %j: c[j],
                'theta-%d_z' %j: c[j]})  
            if j==0:
                damp_dict.update({'theta-%d_x' %j: c[j]})
        trep.forces.Damping(self, default, damp_dict)

    def add_base_constraints(self, base_motion):
        """ Adds constraints on the rotation and translation at whisker's base. """
        BaseAngle(self, 'Driving Angle', base_motion['h'], base_motion['hdt'])
        trep.constraints.PointOnPlane(self, 'Moving Base Point', (1,0,0), 'Base Point',
                name='FX')
        trep.constraints.PointOnPlane(self, 'Moving Base Point', (0,1,0), 'Base Point',
                name='FY') 
        trep.constraints.PointOnPlane(self, 'Moving Base Point', (0,0,1), 'Base Point',
                name='FZ')
        FixedAngle(self, 'theta-0_x', 0.0, name='MX')
        FixedAngle(self, 'theta-0_y', 0.0, name='MY')
        FixedAngle(self, 'theta-0_z', 0.0, name='MZ')


        
            
    
        


        

