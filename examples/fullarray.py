"""
This is a simulation of whiskers contacting a sphere.
"""

import whisker_sim as ws

# Whisker parameters - number of links and names.
N = 10
whisker_names = ['RA1', 'RC1', 'RB2']#['RB4','RB3','RB2','RC2','RC1','RA2','RB1','RA1','RC3','RA3']

# Base motion parameters.
amp    = 1.2      # Amplitude of whisking motion.
freq   = 1.0      # Frequency of whisking motion.
offset = 0.20     # Offset of whisking motion.
base_motion = (amp, freq, offset)

# Surface parameters. (Sphere)
R = .015                           # Radius 
(X, Y, Z) = (.015, .020+1.00, -.005)    # Center position
surface_parameters = {'type': 'sphere', 'center': (X,Y,Z), 'R': R}

# Simulation parameters - final time and timestep.
tf = .11
dt = 0.001

# This function will simulate the dynamics for each whisker, save the
# results to a file, and visualize the results when all the whiskers are
# finished.
ws.simulation.simulate_array(whisker_names, N, base_motion,
                            surface_parameters, tf, dt, visual=True)
