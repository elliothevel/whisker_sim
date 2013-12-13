import numpy as np

"""
 The slot derivatives refer to the arguments of the discrete Lagrangian:
    L2(q1, q2, t1, t2) = (t2-t1)*L((q1+q2)/2, (q2-q1)/(t2-t1), (t1+t2)/2)
 and the discrete forcing:
    f2_minus(q1, q2, t1, t2) = (t2-t1)*f((q1+q2)/2, (q2-q1)/(t2-t1), (t1+t2)/2).
"""

def D1L2(mvi):
    dt = mvi.t2-mvi.t1
    return np.array([0.5*mvi.system.L_dq(q)*dt-mvi.system.L_ddq(q) for q in
        mvi.system.dyn_configs])

def D2L2(mvi):
    dt = mvi.t2-mvi.t1
    return np.array([0.5*mvi.system.L_dq(q)*dt+mvi.system.L_ddq(q) for q in
        mvi.system.dyn_configs])

def D4L2(mvi):
    return mvi.system.L()-np.dot(mvi.system.dq[:mvi.nd], [mvi.system.L_ddq(q)
        for q in mvi.system.dyn_configs])

def D3L2(mvi):
    return -mvi.system.L()+np.dot(mvi.system.dq[:mvi.nd], [mvi.system.L_ddq(q)
        for q in mvi.system.dyn_configs])    

def D2D1L2(mvi):
    dt = mvi.t2-mvi.t1
    return np.array([[0.25*dt*mvi.system.L_dqdq(q1,q2)+0.5*mvi.system.L_ddqdq(q1,q2)
                     -0.5*mvi.system.L_ddqdq(q2,q1)-1.0/dt*mvi.system.L_ddqddq(q1,q2)
                    for q1 in mvi.system.dyn_configs] for q2 in mvi.system.dyn_configs])

def fm2(mvi):
    dt = mvi.t2-mvi.t1
    return np.array([sum([force.f(q)*dt for force in mvi.system.forces]) for q
        in mvi.system.dyn_configs])

def D2fm2(mvi):
    dt = mvi.t2-mvi.t1
    return np.array([[sum([0.5*force.f_dq(q1,q2)*dt+force.f_ddq(q1,q2) for
        force in mvi.system.forces]) for q2 in mvi.system.dyn_configs] for q1
        in mvi.system.dyn_configs])

def D1fm2(mvi):
    dt = mvi.t2-mvi.t1
    return np.array([[sum([0.5*force.f_dq(q1,q2)*dt-force.f_ddq(q1,q2) for force
        in mvi.system.forces]) for q2 in mvi.system.dyn_configs] for q1 in
        mvi.system.dyn_configs])    

def D4fm2(mvi):
    dt = mvi.t2-mvi.t1
    return np.array([[sum([force.f(q1) for force in mvi.system.forces])
        -np.dot(mvi.system.dq[:mvi.nd], [sum([force.f_dq(q1,q2) for force in
            mvi.system.forces]) for q2 in mvi.system.dyn_configs])]
                     for q1 in mvi.system.dyn_configs])
    
def D4D1L2(mvi):
    dt = mvi.t2-mvi.t1
    return np.array([[0.5*mvi.system.L_dq(q1)
        -0.5*np.dot(mvi.system.dq[:mvi.nd], [mvi.system.L_ddqdq(q2,q1) for q2
            in mvi.system.dyn_configs]) +1.0/dt*np.dot(mvi.system.dq[:mvi.nd],
                [mvi.system.L_ddqddq(q1,q2) for q2 in mvi.system.dyn_configs])]
            for q1 in mvi.system.dyn_configs])

def D2D3L2(mvi):
    dt = mvi.t2-mvi.t1
    return np.array([-0.5*mvi.system.L_dq(q1) + np.dot(mvi.system.dq[:mvi.nd],
        [0.5*mvi.system.L_ddqdq(q1,q2) +1.0/dt*mvi.system.L_ddqddq(q1,q2) for
            q2 in mvi.system.dyn_configs]) for q1 in mvi.system.dyn_configs])

def D4D3L2(mvi):
    dt = mvi.t2-mvi.t1
    return 1/dt*np.dot((2*np.array([mvi.system.L_ddq(q1) for q1 in
        mvi.system.dyn_configs])+ np.dot([[mvi.system.L_ddqddq(q1, q2) for q1
            in mvi.system.dyn_configs] for q2 in mvi.system.dyn_configs],
            mvi.system.dq[:mvi.nd])), mvi.system.dq[:mvi.nd])

def h(mvi):
    return np.array([c.h() for c in mvi.system.constraints])

def Dh(mvi):
    dhq = np.array([[c.h_dq(q) for q in mvi.system.dyn_configs] for c in
        mvi.system.constraints])
    return np.reshape(dhq, (mvi.nc, mvi.nd))

def phi(mvi, surfaces=[]):
    if not surfaces:
        return mvi._surface.phi()
    phi_vec = []
    for surface in surfaces:
        phi_vec.append(surface.phi())
    return phi_vec    

def Dphi(mvi, q_i=None, surfaces=[]):
    if not surfaces:
        if q_i is not None:
            return mvi._surface.h_dq(q_i)
        return np.array([mvi._surface.phi_dq(q) for q in mvi.system.dyn_configs])
    p = len(surfaces)
    dphi = np.empty((p,mvi.nd))
    for j in range(p):
        surf = surfaces[j]
        dphi[j] = [surfaces[j].phi_dq(q) for q in mvi.system.dyn_configs]
    return dphi 

def D2LAM2(mvi, con):
    dt = mvi.t1-mvi.t2
    return (mvi.system.lambda_dq(constraint=con) + 1/dt*mvi.system.lambda_ddq(constraint=con))

def D4LAM2(mvi, con):
    dt = mvi.t2-mvi.t1
    return -1/dt*np.dot([mvi.system.lambda_ddq(constraint=con, dq1=q1) for q1
        in mvi.system.dyn_configs], mvi.system.dq[:mvi.nd])

def D2h(mvi):
    Dh = np.zeros(mvi.nc)
    for i in range(mvi.nc):
        constraint = mvi.system.get_constraint(i)
        try:
            Dh[i] = constraint.h_dt()
        except:
            pass
    return Dh

def D1L2_lim(mvi):
    return np.array([-mvi.system.L_ddq(q) for q in mvi.system.dyn_configs])

def D2L2_lim(mvi):
    return np.array([mvi.system.L_ddq(q) for q in mvi.system.dyn_configs])

def DQ(mvi):
    if mvi.t1 != mvi.t2:
        return (mvi.q2-mvi.q1)/(mvi.t2-mvi.t1)
    else:
        return (mvi.q1-mvi.q0)/(mvi.t1-mvi.t0)

def TAU(mvi):
    if mvi.t1 != mvi.t2:
        mvi.set_midpoint()
    else:
        mvi.system.set_state(q=0.5*(mvi.q0+mvi.q1), dq=DQ(mvi), t=0.5*(mvi.t1+mvi.t2))
    return D4L2(mvi)    

def get_dq_and_tau(mvi):
    if mvi.t1 != mvi.t2: 
        mvi.set_midpoint()
        dq = (mvi.q2-mvi.q1)/(mvi.t2-mvi.t1)
    else: 
        dq = mvi._dq1
        mvi.system.set_state(q=mvi.q2, dq=mvi._dq1, t=mvi.t2)
    tau = D4L2(mvi)   
    return dq, tau

def set_system_state(mvi, i):
    """Sets the state of the integrator's system to 1 or 2."""
    if i == 1:
        mvi.system.set_state(q=mvi.q1, dq=mvi._dq1, t=mvi.t1)
    elif i == 2:
        mvi.system.set_state(q=mvi.q2, dq=DQ(mvi), t=mvi.t2)

def save_state(mvi):
    """Saves the integrator's states in a dictionary."""
    return {'s1': mvi._s1, 's2': mvi._s2,
            'lambda0': mvi.lambda0, 'lambda1': mvi.lambda1}

def restore_state(mvi, state):
    """Resets the integrator's state from a dictionary."""
    mvi._s1 = state['s1']
    mvi._s2 = state['s2']

