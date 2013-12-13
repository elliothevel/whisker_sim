"""
This file simuates free (non-contact) whisking in 2D. The motion at the 
base can be any function of time.
"""

import numpy as np
import trep
from trep.visual import *
import whisker_sim as ws

# Import the whisker parameters and specify number of links
N = 15
parameters = ws.whisker.parameters.parameters_from_name('RC1', N)
# Define the motion of the base. We will use a simple sinusoid.
amp  = np.radians(40)   # amplitude in radians
freq = 1                # frequency in Hz.
def h_func(t):
    return amp*np.sin(2*np.pi*freq*t)

def hdt_func(t):
    return 2*np.pi*freq*amp*np.cos(2*np.pi*freq*t)

base_motion = {'h': h_func, 'hdt': hdt_func}

# Create the whisker and the integrator.
whisker = ws.whisker.Whisker3D(parameters, base_motion)
mvi = trep.MidpointVI(whisker)

# Simulate and visualize the whisking.
tf = 1.
dt = 0.01
(t, q, lam) = ws.collisions.sim.simulate(mvi, tf, dt)



simple_whisker = ws.simulation.whiskervisual.SimpleWhisker3D(whisker.num_links)
Q = ws.simulation.whiskervisual.convert_trajectory_for_simple_system(whisker, q)
ws.simulation.visualize_array(simulation_data={'RC1': {'q': Q, 'lam': lam }, 'T': t, 'N': whisker.num_links})
#forces = w.simulation.get_base_forces_3d(whisker, lam, dt)
#print whisker.get_frame('Base Point').p()
#visualize_3d([w.simulation.WhiskerVisual(whisker, t, q)])
#visualize_3d([VisualItem3D(whisker, t, q)])
#visualize_3d([VisualItem3D(whisker, [0,1], [whisker.q, whisker.q])])

