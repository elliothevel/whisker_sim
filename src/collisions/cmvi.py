import trep
import numpy as np
from detect import *
from impact import *
from release import *
from util import *


class CollisionMVI(trep.MidpointVI):
    """ Collision midpoint variational integrator. """

    def __init__(self, system, surfaces, kin_func=None, tolerance=1e-10, 
                 impact_method='root', release_method='endpoint'):
        trep.MidpointVI.__init__(self, system, tolerance=tolerance)

        # Impact surfaces.
        self._surfaces = surfaces
        self._surface = None
        self._state1_impacts = []
        self._state1_releases = []
                  
        # Method for solving release state.          
        assert release_method in ['root', 'interp', 'endpoint'],\
                "Invalid release method"
        self._release_solve_method = release_method
        self._releases_off = False

        # Method for solving impact state.
        assert impact_method in ['root', 'interp', 'endpoint'],\
                "Invalid impact method"
        self._impact_solve_method = impact_method

        # Function for calculating kinematic configuration.
        self.kin_func = kin_func

        # Velocity at state 1.
        self._dq1 = 0.0

        # Multiplier for the previous time-step.
        self.lambda0 = 0.0

    @property
    def surfaces(self):
        """Surfaces the system can contact."""
        return self._surfaces

    @property
    def active_surfaces(self):
        """Constraints currently imposed on the system."""
        return [surface for surface in self._surfaces if surface.active]

    @property
    def inactive_surfaces(self):
        """Constraints not currently imposed on the system."""
        return [surface for surface in self._surfaces if not surface.active]  

    @property
    def _s0(self):
        """State 0 of the integrator."""
        return (self.t0, self.q0, self.p0)

    @_s0.setter
    def _s0(self, value):
        self.t0, self.q0, self.p0 = value
    
    @property
    def _s1(self):
        """State 1 of the integrator."""
        return (self.t1, self.q1, self.p1)

    @_s1.setter
    def _s1(self, value):
        self.t1, self.q1, self.p1 = value

    @property
    def _s2(self):
        """State 2 of the integrator."""
        return (self.t2, self.q2, self.p2)

    @_s2.setter
    def _s2(self, value):
        self.t2, self.q2, self.p2 = value

    def get_surface(self, name):
        """Finds the surface with a specified name."""
        for surface in self._surfaces:
            if surface.name == name:
                return surface
        return None 

    def kin_config(self, t):
        """Returns the kinematic configuration as a function of time."""
        if self.kin_func is None:
            return tuple()
        return self.kin_func(t)

    def cstep(self, t2, k2, *args, **kwargs):
        """
        Steps the system to time t2. This satisfies the dynamics and resolves
        collisions.
        """
        # First have trep step the normal mvi.
        self._s0 = self._s1
        self.lambda0 = self.lambda1
        self.step(t2, k2=self.kin_config(t2), *args, **kwargs)
        self._state1_impacts = []
        self._state1_releases = []

        # Then solve for a new state (t2, q2, p2) consistent with the impact
        # surfaces.
        self.solve_collisions(self.t2)
        
    def solve_collisions(self, t2):
        """Solves impacts & releases to find an acceptable state 2."""
        self.t_start, self.t_end = self.t1, t2
        self.current_dt = self.t_end-self.t_start

        # Loop until an acceptable state is found.
        while True:

            # Check if any events have occurred.
            impact_surfaces  = detect_impacts(self)
            release_surfaces = detect_releases(self)

            # Case 1: no events, return with the current state.
            if not any((release_surfaces, impact_surfaces)):
                return    

            # Case 2: impact but no releases, solve for impacts.
            elif (impact_surfaces and (not release_surfaces)):
                surfaces, impact_state = find_first_impact(self,
                        impact_surfaces, self._impact_solve_method)
                self.impact_update(surfaces, impact_state) 
                self._state1_impacts = surfaces
                #print '\tImpact', impact_state['t'], surfaces
                                
            # Case 3: releases but no impacts, solve for releases.    
            elif (release_surfaces and (not impact_surfaces)):
                surfaces, release_state = find_first_release(self,
                        release_surfaces, self._release_solve_method)
                self.release_update(surfaces, release_state) 
                self._state1_releases = []
                #print '\tRelease', release_state['t'], surfaces

            # Case 4: both impacts and releases, find first event.
            else:
                impact_surfaces, impact_state = find_first_impact(self,
                        impact_surfaces, self._impact_solve_method)
                release_surfaces, release_state = find_first_release(self,
                        release_surfaces, self._release_solve_method)

                # Times are very close to one another - transfer of contact.
                if abs(impact_state['t']-release_state['t'])<1e-3:
                    self.constraint_transition_update(impact_surfaces,
                            release_surfaces, impact_state)
                    #print '\tTransition of constraints'

                # Impact occurred first.
                elif impact_state['t'] < release_state['t']:
                    self.impact_update(impact_surfaces, impact_state)
                    self._state1_impacts = impact_surfaces

                # Release occurred first.
                else:
                    self.release_update(release_surfaces, release_state)


    def impact_update(self, impact_surfaces, impact_state):
        """
        Given a state at the time of impact, this function solves the impact
        map and adds the new constraint to the system. When passed to
        impact_update, state1 is last state before impact. 
        """   
        # Set MVI state 2 to impact state.
        self._s2 = (impact_state['t'], impact_state['q'], impact_state['p'])
        self._tau1 = TAU(self)

        # Apply the discrete impact map.
        dis = discrete_plastic_impact(self, impact_surfaces)
        
        # Add constraints to the system.
        add_constraints(self, impact_surfaces)

    def release_update(self, surfaces, release_state):
        """
        Given a state at the point of release, this function removes the
        constraint(s) from the system and intagrates to the next time step.
        """ 
        # Set the integrator to the release state.
        self._s2 = (release_state['t'], release_state['q'], release_state['p'])

        # Remove the constraints from the system.
        remove_constraints(self, surfaces)
        #remove_list = []
        #for con in self._state1_impacts:
        #    if con in self.system.constraints:
        #        remove_list.append(con)
        #remove_constraints(self, remove_list)        

        self._tau1 = TAU(self)

        # Step to the end of the current timestep.
        if self.t_end != self.t2:
            self.lambda0 = self.lambda1
            self.step(self.t_end, k2=self.kin_config(self.t_end))


    def constraint_transition_update(self, impact_surfaces, release_surfaces,
                                     transfer_state):
        """
        This solves the problem of constraints being added and removed at the
        same time.       
        """
        self._s2 = (transfer_state['t'], transfer_state['q'], transfer_state['p'])

        # Remove the constraints from the system.
        remove_constraints(self, release_surfaces)

        # Apply the impact map to the end of the interval.
        discrete_plastic_impact(self, impact_surfaces) 

        self._state1_impacts = impact_surfaces
        self._state1_releases = release_surfaces
                
    def prepare_to_visualize(self):
        """ 
        Should be called before visualizing system. It will draw the constraint(s)
        even if there are no active surfaces at the end of the simulation.
        """
        self.system.hold_structure_changes()
        for surface in self.inactive_surfaces:
            surface.activate_constraint()
        self.system.resume_structure_changes()  

