import subprocess
import time
import os

def simulate_array(whisker_names, N,  base_motion, surface_parameters, tf, dt,
        visual=True):
    if not os.path.exists('./outputdata'):
        os.makedirs('./outputdata')
    
    open('./outputdata/full_data.p','w').close()
    open('./outputdata/simple_data.p','w').close()
    dirname, filename = os.path.split(os.path.abspath(__file__))
    save_location = os.getcwd()

    amp, freq, offset = base_motion
    s_type = surface_parameters['type']
    X, Y, Z = surface_parameters['center']
    R = surface_parameters['R']

    i = 1
    print 'Simulating array...'
    for name in whisker_names:
        whisker_input = '--whisker={0} --N={1} '.format(name, N)
        motion_input = '--motion={0},{1},{2} '.format(amp, freq, offset)
        sim_input = '--dt={0} --tf={1} '.format(dt,tf)
        surface_input = '--surface={0},{1},{2},{3},{4} '.format(s_type, X, Y, Z, R)
        path_input = '--path=%s' %os.getcwd()
        start_time = time.time()    
        subprocess.call(['python '+dirname+'/simulate.py '+whisker_input+motion_input
                     +sim_input+surface_input+path_input],shell=True)
    
        print '\tSimulated whisker %d of %d - %s (runtime=%f)' %(i,
                len(whisker_names), name, time.time()-start_time)
        i += 1
    print 'Finished simulating array!'

    if visual:
        from whiskervisual import visualize_array
        print 'Preparing visual...'
        visualize_array(surface_parameters, time_scale=dt/0.1)
    print 'Done!'    
