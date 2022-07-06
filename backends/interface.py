# quantum-computed-chemistry/backend/interface.py
# Provides functions to retrieve backends

#qcchem
from .vector_backendwrapper import StatevectorBackend
from .counts_backendwrapper import CountsBackend
from .base_backendwrapper import is_statevector
from .account import all_backends

def get_backendwrapper(name):
    """
    Get the wrapper for the backend with a given name:
    arguments:
        name (str): the name of the desired backend, allowed names are:
            IBMQ backends (e.g. 'ibm_washington', 'ibmq_montreal', 'ibmq_qasm_simulator'),
            'vec': shortcut to statevector simulator,
            'fake(XXX)' where XXX is the name of an IBMQ backend (without prefix), e.g. 'fake(montreal)':
                alternatively 'fake_montreal' or 'qasm_simulator(fake_montreal)'
            'simX' where 'X' is replaced by a number: a simulated backend with an depolarising and readout errors at X times device level
    returns:
        a CountsBackend or StatevectorBackend wrapper to the IBMQ device or simulator
    """
    if isinstance(name,str):
        if is_statevector(name): return StatevectorBackend(name)
        return CountsBackend(name)
    elif is_statevector(name.name()): return StatevectorBackend(name)
    return CountsBackend(name)
