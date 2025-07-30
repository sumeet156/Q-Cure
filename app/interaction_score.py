import numpy as np

def compute_fidelity(statevector1, statevector2):
    fidelity = np.abs(np.dot(np.conj(statevector1), statevector2)) ** 2
    return fidelity
