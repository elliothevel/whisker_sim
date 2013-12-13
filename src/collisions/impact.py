import trep
import numpy as np
import scipy
import scipy.optimize
from detect import *
from util import *


TIME_TOL = 1e-5


def linear_interpolate_impact(mvi):
    """ 
    Finds a time and configuration when phi=0 based on linear
    interpolation. Used to get a guess with which to seed the root-solver.
    """
    set_system_state(mvi, 1)
    phi1 = phi(mvi)
    set_system_state(mvi, 2)
    phi2 = phi(mvi)

    ti = mvi.t1 - (phi1/(phi2-phi1))*(mvi.t2-mvi.t1)
    frac = (ti-mvi.t1)/(mvi.t2-mvi.t1)
    qi = frac*(mvi.q2-mvi.q1)+mvi.q1
    return ti, qi

def solve_impact_time(mvi, surface):
    """ 
    Finds the time and configuration at impact across one surface.
    When passed to this function, mvi should have following states:
        1 - state just before impact
        2 - invalid configuration to be corrected
    This method modifies only state 2 and lambda1.
    """
    mvi._surface = surface
    nd, m = mvi.nd, mvi.nc

    set_system_state(mvi, 1)
    Dh1T = np.transpose(Dh(mvi))
    D2h1 = D2h(mvi)

    # Here, check if the impact occurs exactly at t1. If 
    # | phi | < tolerance, set the impact state at 1. The root-solver
    # can't find this because it will cause a divide by zero error
    # (t2-t1 = 0).
    if abs(phi(mvi))<mvi._surface.tolerance:
        mvi._s2 = (mvi.t1, mvi.q1, mvi.p1)
        mvi.lambda1 = mvi.lambda0   
        return {'t': mvi.t2, 'q': mvi.q2, 'p': mvi.p2, 
                'lambda1': mvi.lambda1}

    t_guess, q_guess = linear_interpolate_impact(mvi)

    # First, if no impacts have occurred at state 1, we use the normal
    # VI equations plus the condition phi(qi) = 0. 
    # f = [DEL(qi, ti), h(qi), phi(qi)].
    if not mvi._state1_impacts:

        x0 = np.hstack((q_guess[:nd], mvi.lambda1, t_guess))

        def func(x):  
            mvi.t2 = x[-1]
            mvi.q2 = np.append(x[:nd], mvi.kin_config(x[-1]))
            mvi.lambda1 = x[nd:(nd+m)]

            mvi.set_midpoint()
            f1 = mvi.p1+D1L2(mvi)+fm2(mvi)-np.dot(Dh1T,mvi.lambda1)

            set_system_state(mvi, 2)
            f2 = h(mvi)
            f3 = phi(mvi)

            return np.hstack((f1, f2, f3))

        def fprime(x):
            mvi.t2 = x[-1]
            mvi.q2 = np.append(x[:nd], mvi.kin_config(x[-1]))
            mvi.lambda1 = x[(nd):(nd+m)]
        
            mvi.set_midpoint()
            Df11 = D2D1L2(mvi) + D2fm2(mvi)
            Df13 = D4D1L2(mvi) + D4fm2(mvi)
        
            set_system_state(mvi, 2)
            Df21 = Dh(mvi)
            Df31 = Dphi(mvi)
            Df23 = D2h(mvi).reshape((mvi.nc, 1))

            Df1 = np.hstack((Df11, -Dh1T, Df13))
            Df2 = np.hstack((Df21, np.zeros((m, m)), Df23))
            Df3 = np.hstack((Df31, np.zeros(m+1)))

            return np.vstack((Df1, Df2, Df3))

    # Otherwise, we need to use the modified impact time solver that 
    # also applies the plastic impact map. Keep in mind that here the 
    # impact surfaces have already been added to the system, i.e., they 
    # are included in h. f = [DEL(qi, ti), h(qi), IMPACT_MAP(qi, ti)]
    else:    
        impact_surfaces = mvi._state1_impacts
        p = len(impact_surfaces)
        tau = mvi._tau1

        x0 = np.hstack((q_guess[:nd], mvi.lambda1, t_guess, 
                        np.dot(D2h1, mvi.lambda1)-tau))

        def func(x):  
            mvi.t2 = x[nd+m]
            mvi.q2 = np.append(x[:nd], mvi.kin_config(mvi.t2))
            mvi.lambda1 = x[nd:(nd+m)]
            Ec = x[-1]

            mvi.set_midpoint()
            f1 = mvi.p1+D1L2(mvi)+fm2(mvi)-np.dot(Dh1T,mvi.lambda1)
            f3 = tau+D3L2(mvi)-np.dot(D2h1, mvi.lambda1)+Ec

            set_system_state(mvi, 2)
            f2 = h(mvi)
            f4 = phi(mvi)
               
            return np.hstack((f1, f2, f3, f4))

        def fprime(x):
            mvi.t2 = x[nd+m]
            mvi.q2 = np.append(x[:nd], mvi.kin_config(mvi.t2))
            mvi.lambda1 = x[nd:(nd+m)]
            Ec = x[-1]
        
            mvi.set_midpoint()
            Df11 = D2D1L2(mvi) + D2fm2(mvi)
            Df13 = D4D1L2(mvi) + D4fm2(mvi)
            Df31 = D2D3L2(mvi)
            Df33 = D4D3L2(mvi)
        
            set_system_state(mvi, 2)
            Df21 = Dh(mvi)
            Df23 = D2h(mvi).reshape((m, 1))
            Df32 = -D2h(mvi)
            Df41 = Dphi(mvi)

            Df1 = np.hstack((Df11, -Dh1T, Df13, np.zeros((nd,1))))
            Df2 = np.hstack((Df21, np.zeros((m,m)), Df23, np.zeros((m,1))))
            Df3 = np.hstack((Df31, Df32, Df33, (1.0,)))
            Df4 = np.hstack((Df41, np.zeros(m+2)))

            return np.vstack((Df1, Df2, Df3, Df4))

    # Call the root-solver to find the time of impact.
    try:
        (x, infodict, ier, mesg) = \
                scipy.optimize.fsolve(func, x0, full_output=True, fprime=fprime)
    except ZeroDivisionError: 
        raise Exception("Divide by zero in impact time solver.")

    if ier != 1:
        print mesg
        raise trep.ConvergenceError("Find impact failed to converge.") 

    mvi.t2 = x[nd+m]
    mvi.q2 = np.append(x[:nd], mvi.kin_config(mvi.t2))
    mvi.lambda1 = x[nd:(nd+m)]
    mvi.set_midpoint()    
    mvi.p2 = D2L2(mvi)

    return {'t': mvi.t2, 'q': mvi.q2, 'p': mvi.p2,
            'lambda1': mvi.lambda1}
  

def discrete_plastic_impact(mvi, impact_surfaces=[]):
    """ 
    Finds the configuration after impact.
    when passed to this function, mvi should have following states
        1 - state just before impact
        2 - impact state
    plastic_impact updates so that state2 is post-impact, state1 is impact     
    """
    p = len(impact_surfaces)
    nd, m = mvi.nd, mvi.nc

    mvi._s2 = (mvi.t2, mvi.q2, mvi.p2)
    mvi.lambda0 = mvi.lambda1
    tau = mvi._tau1

    # It is always possible that an impact happens very close to the 
    # time step t2. If this is the case, then we need to take a whole new
    # step to ensure the impact map is computed. This step is assumed to 
    # be the same size as the previous one. Not ideal, but definitely 
    # needed.
    if abs(mvi.t2-mvi.t_end)<TIME_TOL:
        mvi.t_end = mvi.t_start+2*mvi.current_dt 
    mvi.t2 = mvi.t_end

    set_system_state(mvi, 1)
    Dh1T  = np.transpose(Dh(mvi))
    D2h1  = D2h(mvi) 

    # If there are impact surfaces, apply normal plastic impact update.
    if p != 0:

        Dphi1 = Dphi(mvi, surfaces=impact_surfaces)
        x0 = np.hstack((mvi.q2[:nd], mvi.lambda1, np.dot(Dphi1, mvi.p1),
                        np.dot(D2h1, mvi.lambda1)-tau))

        def func(x):
            mvi.q2 = np.append(x[:nd], mvi.kin_config(mvi.t2))
            mvi.lambda1 = x[nd:(nd+m)]
            lambda_c = x[(nd+m):-1]
            Ec = x[-1]        

            set_system_state(mvi, 2)
            f2 = h(mvi)
            f4 = phi(mvi, surfaces=impact_surfaces)
            
            mvi.set_midpoint()
            f1 = (mvi.p1+D1L2(mvi)+fm2(mvi)-np.dot(Dh1T,mvi.lambda1)-
                  np.dot(np.transpose(Dphi1),lambda_c))
            f3 = tau+D3L2(mvi)-np.dot(D2h1, mvi.lambda1)+Ec

            return np.hstack((f1, f2, f3, f4))

        def fprime(x):
            mvi.q2 = np.append(x[:nd], mvi.kin_config(mvi.t2))
            mvi.lambda1 = x[nd:nd+m]
            lambda_c = x[(nd+m):-1]
            Ec = x[-1]

            mvi.set_midpoint()
            Df11 = D2D1L2(mvi) + D2fm2(mvi)
            Df31 = D2D3L2(mvi)

            set_system_state(mvi, 2)
            Df21 = Dh(mvi)
            Df41 = Dphi(mvi, surfaces=impact_surfaces)
            Df32 = -D2h(mvi)

            Df1 = np.column_stack((Df11, -Dh1T, -np.transpose(Dphi1), 
                                   np.zeros((nd,1))))
            Df2 = np.hstack((Df21, np.zeros((mvi.nc, mvi.nc+p+1))))
            Df3 = np.hstack((Df31, Df32, np.zeros(p), (1.0,)))
            Df4 = np.hstack((Df41,  np.zeros((p,mvi.nc+p+1))))

            return np.vstack((Df1, Df2, Df3, Df4))

    # Because of multiple events, it is possible we need to solve the
    # impact update after the constraints have already been added to the 
    # system. This is the same problem, except lambda_c is a part of lambda
    # already.
    else:
        x0 = np.hstack((mvi.q2[:nd], mvi.lambda1, 
                        np.dot(D2h1, mvi.lambda1)-tau))

        def func(x):
            mvi.q2 = np.append(x[:nd], mvi.kin_config(mvi.t2))
            mvi.lambda1 = x[nd:(nd+m)]
            Ec = x[-1]        

            set_system_state(mvi, 2)
            f2 = h(mvi)
            
            mvi.set_midpoint()
            f1 = mvi.p1+D1L2(mvi)+fm2(mvi)-np.dot(Dh1T,mvi.lambda1)
            f3 = tau+D3L2(mvi)-np.dot(D2h1, mvi.lambda1)+Ec

            return np.hstack((f1, f2, f3))

        def fprime(x):
            mvi.q2 = np.append(x[:nd], mvi.kin_config(mvi.t2))
            mvi.lambda1 = x[nd:nd+m]
            Ec = x[-1] 

            mvi.set_midpoint()
            Df11 = D2D1L2(mvi) + D2fm2(mvi)
            Df31 = D2D3L2(mvi)

            set_system_state(mvi, 2)
            Df21 = Dh(mvi)
            Df32 = -D2h(mvi)

            Df1 = np.column_stack((Df11, -Dh1T, np.zeros((nd,1))))
            Df2 = np.hstack((Df21, np.zeros((m, m+1))))
            Df3 = np.hstack((Df31, Df32, (1.0,)))

            return np.vstack((Df1, Df2, Df3))


    (x, infodict, ier, mesg) = \
            scipy.optimize.fsolve(func, x0, full_output=True, fprime=fprime)

    if (ier != 1):
        if np.linalg.norm(func(x))<mvi.tolerance:
            pass
        else:
            print mesg  
            raise trep.ConvergenceError("Discrete plastic impact update failed\
                                         to converge.")  

    mvi.q2 = np.append(x[:nd], mvi.kin_config(mvi.t2))
    mvi.lambda1 = x[nd:(nd+m)]
    mvi.set_midpoint()    
    mvi.p2 = D2L2(mvi)
    mvi.tau2 = D4L2(mvi)

    return {'t': mvi.t2, 'q': mvi.q2, 'lambda1': mvi.lambda1, 
            'lambdac': x[(nd+m):-1], 'Ec': x[-1]}

 
def find_first_impact(mvi, impact_surfaces, method):
    """ 
    When given one or more impacts, finds the first impact to occur and the
    state at that point. Returns the surface(s) and state. 
    solve impact returns mvi with the same states with which it was passed. 
    """
        
    # Endpoint method: all impacts occur at state 1. 
    if method == 'endpoint':
        #TODO: make sure system state is correct for calculating D4L2.
        state = {'t': mvi.t1, 'q': mvi.q1, 'p': mvi.p1, 
                 'lambda1': mvi.lambda0, 'tau': D4L2(mvi)}
        return impact_surfaces, state

    # Interp method: linearly interpolate the time of release.
    elif method == 'interp':

        times = []
        for surface in impact_surfaces:
            mvi._surface = surface
            times.append(linear_interpolate_impact(mvi)[0])
        ti = np.mean(times)

        if abs(mvi.t1-ti)<TIME_TOL:
            state = {'t': mvi.t1, 'q': mvi.q1, 'p': mvi.p1,
                     'lambda1': mvi.lambda0, 'tau': D4L2(mvi)}
            return impact_surfaces, state

        initial_state = save_state(mvi)
        mvi.set_single_state(2, mvi.t1, mvi.q1, mvi.p1)
        mvi.step(ti, k2=mvi.kin_config(ti))
        mvi.set_midpoint()
        state = {'t': mvi.t2, 'q': mvi.q2, 'p': mvi.p2, 
                 'lambda1': mvi.lambda1, 'tau': D4L2(mvi)}
        restore_state(mvi, initial_state)

        return impact_surfaces, state


    # Root method: root-solve the time of impact. Solve separately and then
    # check for multiple impacts.
    elif method == 'root':

        # First solve for each impact separately, saving the states.
        state = save_state(mvi)
        impacts = {}
        for surface in impact_surfaces:
            s = solve_impact_time(mvi, surface)
            impacts.update({surface: s})
            restore_state(mvi, state)

        # Find the first surface to impact assuming distinct impacts.
        impact_surface = impacts.keys()[np.argmin([impacts[k]['t'] 
                                        for k in impacts.keys()])]
        impact_state = impacts[impact_surface]

        # Finally, check for multiple impacts. Multiple impacts are said to 
        # have occurred if the impact configuration lies on more that one 
        # surface.
        surfaces = [impact_surface]
        mvi.system.set_state(q=impact_state['q'], t=impact_state['t'])
        for surface in impact_surfaces:
            if surface != impact_surface:
                if abs(surface.phi())<surface.tolerance:
                    surfaces.append(surface)         

        return surfaces, impact_state   

def add_constraints(mvi, impact_surfaces):
    """Add constraints to the system."""
    state = save_state(mvi) 

    mvi.system.hold_structure_changes()
    for surface in impact_surfaces:
        surface.activate_constraint()
    mvi.system.resume_structure_changes()

    restore_state(mvi, state)
    mvi.lambda0 = np.append(state['lambda0'], [0.0]*len(impact_surfaces)) 
    mvi.lambda1 = np.append(state['lambda1'], [0.0]*len(impact_surfaces))

def continuous_plastic_impact(sys, surfaces, q1, dq1):
    """ 
    Fully continuous impact update. q1 = q(t-), dq1 = qdot(t-)  
    """
    p = len(surfaces)
    nd = sys.nQd

    sys.q = q1
    sys.dq = dq1
    
    p1 = np.array([sys.L_ddq(q1) for q1 in sys.dyn_configs])
    Dh = np.reshape(np.array([[c.h_dq(q) for q in sys.dyn_configs] 
                               for c in sys.constraints]), (sys.nc, sys.nQd)) 
    Dphi = np.array([[surf.phi_dq(q) for q in sys.dyn_configs]
                      for surf in surfaces])
    Dh2 = np.append(Dh,  Dphi).reshape(sys.nc+p,nd)
    Dh2T = np.transpose(Dh2)

    def func(x):
        sys.dqd = x[:nd]
        lam = x[nd:]
        p2 = np.array([sys.L_ddq(q1) for q1 in sys.dyn_configs]) 
        f1 = p2-p1-np.dot(Dh2T, lam)
        f2 = np.dot(Dh2, sys.dqd)
        return np.hstack((f1, f2))

    def fprime(x):
        sys.dqd = x[:nd]
        lam = x[nd:]
        Df11 = np.array([[sys.L_ddqddq(q1, q2) for q2 in sys.dyn_configs]
                          for q1 in sys.dyn_configs])
        Df1 = np.hstack((Df11, -Dh2T))
        Df2 = np.hstack((Dh2, np.zeros((sys.nc+p, sys.nc+p))))
        Df = np.vstack((Df1, Df2))
        return Df

    dq_guess = dq1[:nd]
    x0 = np.hstack((dq_guess, np.zeros(p+sys.nc)))

    (x, infodict, ier, mesg) = \
            scipy.optimize.fsolve(func, x0, full_output=True, fprime=fprime)

    if (ier != 1) and np.linalg.norm(func(x))>(1e-10):
        print mesg
        raise trep.ConvergenceError("Continuous plastic impact update failed\
                                     to converge.")

    sys.dqd = x[:nd]
    lam = x[nd:]

    return sys.dqd

