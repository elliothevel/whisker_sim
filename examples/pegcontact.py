import numpy as np
import trep
from trep.visual import *
import matplotlib.pyplot as plt
import whisker_sim as ws

# Import the single whisker parameters. (or use different parameters)
parameters = ws.whisker.parameters.SINGLE_WHISKER_DEFAULT_PARAMETERS

# Define the motion of the base. We will use a simple sinusoid.
amp  = np.radians(10)   # amplitude in radians
freq = 8                # frequency in Hz.
def h_func(t):
    return amp*np.sin(2*np.pi*freq*t)

def hdt_func(t):
    return 2*np.pi*freq*amp*np.cos(2*np.pi*freq*t)

base_motion = {'h': h_func, 'hdt': hdt_func}

# Create the whisker.
whisker = ws.whisker.Whisker2D(parameters, base_motion)

# Creat the peg surface.
L = sum(whisker.lengths)
Y, Z = 0.8*L, 0.05*L
surfaces = ws.collisions.surfaces.global_surfaces(
        ws.collisions.surfaces.Peg, whisker, kwargs={'Y': Y, 'Z': Z})

# Create the integrator.
mvi = ws.collisions.cmvi.CollisionMVI(whisker, surfaces)

# Simulate and visualize the whisking.
tf = 1.0
dt = 0.001
(t, q, lam) = ws.collisions.sim.simulate(mvi, tf, dt)

mvi.prepare_to_visualize()   # called to draw the peg
t = np.array(t)/0.05         # scaling the time will slow down the animation

visualize_3d([VisualItem3D(whisker, t, q)], 
             camera_pos=[.12, 5.97e-16, 0.], 
             camera_ang=[3.14159, 0., 0.])

# Uncomment the following to get the base forces from the simulation. Do not
# use trep's visualization and a matplotlib figure at once though, due to some
# bug in the graphics backend.
base_forces = ws.simulation.simulate.get_base_forces_2d(whisker, lam, dt)
plt.plot(t, base_forces['Fy']*1e6, t, base_forces['Fz']*1e6)
plt.legend(['Axial Force','Transverse Force'])
plt.xlabel('time (sec)')
plt.ylabel('Force (micro-N)')
#plt.show()
