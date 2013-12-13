"""
This file simuates free-air (non-contact) whisking in 2D. The motion at the 
base can be any function of time. The simulation is run with different
intrinsic curvatures and the resulting axial forces are plotted.
"""
import numpy as np
import trep
from trep.visual import *
import matplotlib.pyplot as plt
import matplotlib
import whisker_sim as ws

# Import the single whisker parameters.
parameters = ws.whisker.parameters.SINGLE_WHISKER_DEFAULT_PARAMETERS

# Define the motion of the base. We will use a simple sinusoid.
amp  = np.radians(10)   # amplitude in radians
freq = 8                # frequency in Hz.
def h_func(t):
    return amp*np.sin(2*np.pi*freq*t)

def hdt_func(t):
    return 2*np.pi*freq*amp*np.cos(2*np.pi*freq*t)

base_motion = {'h': h_func, 'hdt': hdt_func}

# Create the whisker and the integrator.
fx = []
for a in [0, 10, 20, 30]:
    parameters['a'] = a
    whisker = ws.whisker.Whisker2D(parameters, base_motion)
    mvi = trep.MidpointVI(whisker)

    # Simulate and visualize the whisking.
    tf = 2.0
    dt = 0.005
    (t, q, lam) = ws.collisions.sim.simulate(mvi, tf, dt)

    base_forces = ws.simulation.simulate.get_base_forces_2d(whisker, lam, dt)
    fx.append(base_forces['Fy']*1e6)

# Plot the results for each curvature.
N = 200
plt.rc('text', usetex=True)
matplotlib.rcParams.update({'font.size':18})
fig = plt.figure(facecolor='white',figsize=(12,5))
ax = plt.axes()
ax.plot(t[:N], fx[0][:N], t[:N], fx[1][:N], t[:N], fx[2][:N], t[:N], fx[3][:N])
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
ax.set_xlabel('$t\;(\mathrm{sec})$')
ax.set_ylabel('$F_{\mathrm{axial}} \; (\mu\mathrm{N})$')
plt.legend(['$a=0$','$a=10$','$a=20$','$a=30$'])
plt.show()
