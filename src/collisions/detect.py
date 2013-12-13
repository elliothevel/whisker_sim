from util import *

def detect_impacts(mvi):
    """Detects which surfaces have been crossed."""
    set_system_state(mvi, 2)
    impact_list = []
    for surface in mvi.inactive_surfaces: 
        if surface.phi() < 0:
            impact_list.append(surface)     
    return impact_list

def detect_releases(mvi, q=None, dq=None):
    """
    Detects constraint releases by checking the discrete-time
    Lagrange multipliers.
    """
    if mvi._releases_off:
        return []

    release_list = []
    lam2 = mvi.lambda1
    for surface in mvi.active_surfaces:
        if lam2[surface.index]*surface.sign < 0.0:
            release_list.append(surface)
    return release_list  

def sign_change(a, b):
    if (a > 0) and (b < 0):
        return True
    if (a < 0) and (b > 0):
        return True
    return False

