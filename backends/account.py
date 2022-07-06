# quantum-computed-chemistry/backend/account.py
# Handles connecting to IBMQ backends

#core
from warnings import warn
from importlib import import_module

#qiskit
from qiskit import IBMQ, Aer
from qiskit.providers.aer import QasmSimulator, noise

#define ACCOUNT as a global variable so that login is required only once per session
global ACCOUNT
ACCOUNT=None

def load(**kwargs):
    """
    load IBMQ account:
    arguments:
        kwargs: can be used to specify a provider (e.g. hub, group, project), (default is; hub=ibm-q-melbourne)
    returns:
        the qiskit account
    """
    global ACCOUNT
    if not len(kwargs): kwargs={'hub':'ibm-q-melbourne'}
    IBMQ.load_account()
    ACCOUNT=IBMQ.get_provider(**kwargs)
    return ACCOUNT

def all_backends(name=None,**kwargs):
    """
    Retrieve a list of all available backends
    arguments:
        name (str): the name of a backend as a string, for allowed names see get_backend
        **kwargs: can be used to specify a backend by properties other than it's name (e.g. 'hub')
    returns:
        tuple of IBMQ backends matching the criteria
    """
    global ACCOUNT
    if not ACCOUNT: load()
    if name: kwargs['name']=name
    backends=ACCOUNT.backends(**kwargs)
    if not len(backends): raise ValueError(f'No available backends match the specified criteria: {kwargs}')
    return backends

def get_backend(name):
    """
    retrieve a single backend
    arguments:
        name (str): the name of the desired backend, allowed names are:
            IBMQ backends (e.g. 'ibm_washington', 'ibmq_montreal', 'ibmq_qasm_simulator'),
            'vec': shortcut to statevector simulator,
            'fake(XXX)' where XXX is the name of an IBMQ backend (without prefix), e.g. 'fake(montreal)':
                alternatively 'fake_montreal' or 'qasm_simulator(fake_montreal)'
            'simX' where 'X' is a number: a simulated backend with depolarising and readout errors at X times device level
    returns: the first backend matching the criteria
    """
    if name[:5].lower()=='fake_': return _get_aer_from_device(name[5:])
    elif name[:5].lower()=='fake(': return _get_aer_from_device(name[5:-1])
    elif name[:20].lower()=='qasm_simulator(fake_': return _get_aer_from_device(name[20:-1])
    elif name[:3].lower()=='sim': return _get_aer_from_model(name)
    elif name.lower() in ('statevector_simulator','vec'): return _get_aer_statevector('statevector_simulator')
    else:
        backends=all_backends(name)
        if len(backends)>1: warn(f'Multiple backends found with the name "{name}"')
        try: return backends[0]
        except IndexError: raise(f'No backends found with the name "{name}"')

def _get_aer_from_device(name):
    """
    return a simulated version of the device (should be called through get_backend)
    """
    return QasmSimulator.from_backend(getattr(import_module('qiskit.test.mock'),'Fake'+name[0].upper()+name[1:])(),method='automatic')
    
def _get_aer_from_model(name):
    """
    return a simulated device with depolarising and readout error (should be called through get_backend)
    """
    try: level=float(name[3:])
    except ValueError: (level,name)=(0,'sim0')
    p1=0.001*level
    p2=0.01*level
    p01=0.03*level
    p10=0.015*level
    limited=set()
    if level<0: raise ValueError('Cannot have a negative error level')
    if p1>1: limited.add('1-qubit gate'); p1=1
    if p2>1: limited.add('2-qubit gate'); p2=2
    if p10>0.5: limited.add('Readout 1->0'); p10=0.5
    if p01>0.5: limited.add('Readout 0->1'); p01=0.5
    if len(limited): warn('The following error types cannot be increased further: '+', '.join(limited))
    # Depolarizing quantum errors
    error_1 = noise.depolarizing_error(p1, 1)
    error_2 = noise.depolarizing_error(p2, 2)
    # Readout error
#     error_r = noise.ReadoutError([[1-p01,p01],[p10,1-p10]])
    error_r = noise.ReadoutError([[1-p10,p10],[p01,1-p01]])
    # Add errors to noise model
    noise_model = noise.NoiseModel()
    noise_model.add_all_qubit_quantum_error(error_1, ['u1', 'u2', 'u3'])
    noise_model.add_all_qubit_quantum_error(error_2, ['cx'])
    noise_model.add_all_qubit_readout_error(error_r)
    sim=QasmSimulator(noise_model=noise_model)
    sim.custom_name=name
    return sim
def _get_aer_statevector(name):
    """
    return a statevector simulator
    """
    return Aer.get_backend(name)

