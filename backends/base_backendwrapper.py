# quantum-computed-chemistry/backends.base_backendwrapper
# A base class for backend wrappers, defines initialisation and methods related to checking for statevector capabilities and whether the backend is simulated

#qiskit
from qiskit import assemble

#qcchem
from .account import get_backend

class BackendWrapper():
    """
    A base class for backend wrappers. Should not be used directly, defines properties common to CountsBackend and StatevectorBackend
    """
    def __init__(self,name):
        if isinstance(name,str): self.backend=get_backend(name)
        else: self.backend=name
    def __getattr__(self,attr):
        if attr=='name':
            try: return self.backend.custom_name
            except AttributeError: return self.backend.name()
        if attr=='n_qubits':
            try: return len(self.backend.properties().qubits)
            except AttributeError: return None
        if attr=='connectivity':
            try: gates=self.properties().gates
            except AttributeError: return None
            connectivity={}
            for g in gates:
                if len(g.qubits)==2:
                    q0,q1=g.qubits
                    if q0 in connectivity: connectivity[q0]+=(q1,)
                    else: connectivity[q0]=(q1,)
            return connectivity
        if attr=='backend': return self.__getattribute__(attr)
        return getattr(self.backend,attr)
    @property
    def is_aer(self):
        return ('qasm_simulator' in self.backend.name() and not self.name=='ibmq_qasm_simulator') or self.is_statevector
    @property
    def is_statevector(self):
        return is_statevector(self.name)
    @property
    def is_noiseless(self):
        return self.is_statevector or self.name=='sim0'
    def circuits2jobs(self,circuits,shots=8192):
        """
        Batch (if required) and submit a set of circuits to the backend
        arguments:
            circuits (iterable of qcchem.Circuit and/or qiskit.QuantumCircuit): the circuits to be submitted
            shots (int, default 8192): the number of measurements to be taken in each basis
        returns:
            tuple of qiskit jobs: the jobs being run on the device
        """
        if self.is_aer: max_circs=len(circuits)
        else:
            try: max_circs=self.backend.configuration().max_experiments
            except AttributeError: max_circs=300
        chunked_circuits=[]
        while len(circuits)>max_circs:
            chunked_circuits.append(circuits[:max_circs])
            circuits=circuits[max_circs:]
        chunked_circuits.append(circuits)
        qobjs=(assemble(list(circuit_chunk),self.backend,shots=shots) for circuit_chunk in chunked_circuits)
        # note that assemble won't accept a tuple of circuits it has to be a list
        return tuple(self.backend.run(qobj) for qobj in qobjs)
    
def is_statevector(name):
    return name in ('statevector_simulator','vec')
