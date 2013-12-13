import argparse
import pickle
import numpy as np
import whisker_sim as ws
from whisker_sim.whisker.parameters import  parameters_from_name
from whisker_sim.whisker import Whisker3D
from whiskervisual import convert_trajectory_for_simple_system


def build_and_simulate_whisker(name, N, tf, dt, surface, base_motion):
    """ Builds the whisker and simulates. """
    parameters = parameters_from_name(name, N)
    whisker = Whisker3D(parameters, base_motion)
    cmvi = ws.collisions.CollisionMVI(whisker, surface)
    cmvi.initialze_state(0.0, whisker.q, np.zeros(cmvi.nd))
    return ws.simulation.simulate(cmvi, tf, dt)


def save_data(whisker, Q, T, LAM, RXN, path):
    """ Saves the force and trajectory data to a file. """
    try:
        f = open(path+'/outputdata/full_data.p','r+')
        full_data = pickle.load(f)
        f.close()
    
    except:
        f = open(path+'/outputdata/full_data.p','wb')
        f.close()
        full_data = {}

    full_data[whisker.name] = {'q': Q, 'lam': LAM, 'rxn': RXN}

    if 'T' not in full_data:
        full_data['T'] = T 

    f = open(path+'/outputdata/full_data.p','wb')
    pickle.dump(full_data, f)
    f.close()
        
    try:
        f = open(path+'/outputdata/simple_data.p','r+')
        simple_data = pickle.load(f)
        f.close()
    
    except:
        f = open(path+'/outputdata/simple_data.p','wb')
        f.close()
        simple_data = {}

    Q_simple = convert_trajectory_for_simple_system(whisker, Q)
    simple_data[whisker.name] = {'q': Q_simple}

    if 'T' not in simple_data:
        simple_data['T'] = T

    if 'N' not in simple_data:
        simple_data['N'] = whisker.num_links

    f = open(path+'/outputdata/simple_data.p','wb')
    pickle.dump(simple_data, f)
    f.close()
    

def get_base_forces_2d(whisker, lam, dt):
    """ Extracts the base forces from lambda. """
    M_index = whisker.get_constraint('M').index
    Fy_index = whisker.get_constraint('FY').index
    Fz_index = whisker.get_constraint('FZ').index
    Fy, Fz = [], []
    M = []
    for j in range(len(lam)):
        Fy.append(lam[j][Fy_index]/dt)
        Fz.append(lam[j][Fz_index]/dt)
        M.append(lam[j][M_index]/dt)

    base_forces = {'Fy': np.array(Fy), 'Fz': np.array(Fz), 'M': np.array(M)}        
    return base_forces


def get_base_forces_3d(whisker, lam, dt):
    """ Extracts base forces from lambda. """
    Mx_index = whisker.get_constraint('MX').index
    My_index = whisker.get_constraint('MY').index
    Mz_index = whisker.get_constraint('MZ').index
    Fx_index = whisker.get_constraint('FX').index
    Fy_index = whisker.get_constraint('FY').index
    Fz_index = whisker.get_constraint('FZ').index
    Fx, Fy, Fz = [], [], []
    Mx, My, Mz = [], [], []
    for j in range(len(lam)):
        Fx.append(lam[j][Fx_index]/dt)
        Fy.append(lam[j][Fy_index]/dt)
        Fz.append(lam[j][Fz_index]/dt)
        Mz.append(lam[j][Mz_index]/dt)
        Mx.append(lam[j][Mx_index]/dt)
        My.append(lam[j][My_index]/dt)

    base_forces = {'Fx': np.array(Fx), 'Fy': np.array(Fy), 'Fz': np.array(Fz),
                   'Mx': np.array(Mx), 'My': np.array(My), 'Mz': np.array(Mz)}        
    return base_forces
    

# When simulating the array, the entire file is called so that the systems are
# destroyed after each run.
if __name__ == '__main__':
    
    # Create a parser to get the inputs to the file.
    parser = argparse.ArgumentParser(description='Gets whisker and simulation parameters')
    parser.add_argument('--tf', type=float)
    parser.add_argument('--dt', type=float)
    parser.add_argument('--N', type=int)
    parser.add_argument('--surface',type=str)
    parser.add_argument('--whisker')
    parser.add_argument('--motion')
    parser.add_argument('--path',type=str)
    args = vars(parser.parse_args())

    # Now make sense of the inputs and create relevant variables.
    tf = args['tf']
    dt = args['dt']
    N  = args['N']
    sp = args['surface'].split(',')
    surface = {'type':sp[0], 'X': float(sp[1]), 'Y': float(sp[2]), 'Z': float(sp[3]), 
            'R': float(sp[4])}
    whisker_name = args['whisker']
    m = args['motion'].split(',')
    motion = [float(m[0]), float(m[1]), float(m[2])]

    amp  = motion[0]
    freq = motion[1]
    offset = motion[2]

    def h_func(t):
        return amp*np.sin(2*np.pi*freq*t)-offset

    def hdt_func(t):
        return 2*np.pi*freq*amp*np.cos(2*np.pi*freq*t)
  
    base_motion = {'h': h_func, 'hdt': hdt_func}

    parameters = parameters_from_name(whisker_name, N)
    whisker = Whisker3D(parameters, base_motion, name=whisker_name)

    # Create the surface and the base motion functions.
    if surface['type'] == 'sphere':
        (X, Y, Z, R) = (surface['X'], surface['Y'],
                        surface['Z'], surface['R'])
    #surf = tc.surfaces.Sphere(X, Y, Z, R, whisker)
    surfaces = ws.collisions.surfaces.global_surfaces(ws.collisions.surfaces.Sphere,
            whisker, kwargs={'X': X, 'Y': Y, 'Z': Z, 'R': R})

    # Next, run the simulation with the given parameters.
    cmvi = ws.collisions.CollisionMVI(whisker, surfaces)
    cmvi.initialize_from_state(0.0, whisker.q, np.zeros(cmvi.nd))
    (T, Q, LAM) = ws.collisions.sim.simulate(cmvi, tf, dt)
    
    RXN = get_base_forces_3d(whisker, LAM, dt)
    save_data(whisker, Q, T, LAM, RXN, args['path'])
