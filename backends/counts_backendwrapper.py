# quantum-computed-chemistry/backends/counts_backendwrapper.py
# A class to wrap backends that return measurement counts (as opposed to a statevector)

#core
from time import sleep

#qiskit
from qiskit import transpile,ClassicalRegister

#qcchem
from .base_backendwrapper import BackendWrapper


class CountsBackend(BackendWrapper):
    """
    A class to wrap backends that return measurement counts (as opposed to a statevector)
        main methods:
            CountsBackend.run
            CountsBackend.circ2circuits
            CountsBackend.circuits2jobs
            CountsBackend.jobs2counts
    """
    
    def circ2circuits(self,circuit,bases):
        """
        Apply a set of measurements to a base circuit
        arguments:
            circuit (qcchem.Circuit or qiskit.QuantumCircuit): the state preparation circuit to be measured
            bases (iterable of qcchem.measurement.PauliBasis and/or qcchem.measurement.FermionicBasis): the bases in which the prepared state should be measured
        returns:
            tuple of circuits (matching the type of input circuit): the complete circuits to be run on the backend, if bases is ordered (e.g. tuple or list as opposed to set) the order of the output circuits will match the order of the input bases
        """
        if len(circuit.clbits)==0: circuit.add_classicalregister()
        circuit=transpile(circuit,self.backend)
        return tuple(transpile(basis.apply(circuit.copy()),self.backend) for basis in bases)
    def jobs2counts(self,jobs,bases=None):
        """
        Retrieve the results for a previously submitted batch of jobs
        arguments:
            jobs (iterable of qiskit jobs and/or str): the jobs for which results should be retrieved, if given as strings, they are assumed to be job IDs for jobs run on this device
            bases, optional (iterable of qcchem.measurement.PauliBasis and/or qcchem.measurement.FermionicBasis): the bases being measured in each circuit if provided, these are used to modify the counts to remove any non-locality in fermionic measurements
        returns:
            bases, if bases is provided (the same object containing the same bases): the results can be retrieved using bases[ix].results
            results, if bases is not provided (tuple of dictionaries): each dictionary has bitstrings as keys and frequencies as values, the dictionaries are in the order the circuits were initially submitted
        """
        results=()
        if isinstance(jobs[0],str): jobs=[self.retrieve_job(job) for job in jobs]
        for job in jobs:
            while not job.in_final_state(): sleep(5)
            results+=(job.result().get_counts(),)
        results=(r for res in results for r in (res if isinstance(res,list) else [res,]))
        if bases is None: return tuple(results)
        for basis,res in zip(bases,results):
            basis.results=res
        return bases
    def run(self,circuit,bases='def',shots=8192):
        """
        Convert a state preparation circuit to counts in a given set of measurement bases
        see also:
            CountsBackend.circ2circuits
            CountsBackend.circuits2jobs
            CountsBackend.jobs2counts
        arguments:
            circuit (qcchem.Circuit or qiskit.QuantumCircuit): the state preparation circuit to be measured
            bases (iterable of qcchem.measurement.PauliBasis and/or qcchem.measurement.FermionicBasis): the bases in which the prepared state should be measured
            shots (int, default 8192): the number of measurements to be taken in each basis
        returns:
            bases (same as input): the results can be retrieved using bases[ix].results
        """
        if isinstance(circuit,(list,tuple)): return type(circuit)(self.run(circ,bases,shots) for circ in circuit)
        if bases=='def': bases=[PauliBasis('Z'*circuit.num_qubits)]
        elif bases is None: bases=[PauliBasis('I'*circuit.num_qubits)]
        return self.jobs2counts(self.circuits2jobs(self.circ2circuits(circuit,bases),shots),bases)
    
from itertools import zip_longest

class CommutationError(Exception):
    pass

class PauliBasis():
    def __init__(self,elements):
        self.basis=elements
    def apply(self,circuit,barriers=False):
        if self.active is None: self.initialise()
        self._apply(circuit,barriers)
        if not len(circuit.clbits): circuit.add_classicalregister()
        active=sorted(self.active)
        circuit.measure(active,[len(circuit.clbits)-i-1 for i in active])
        return circuit
    def __lt__(self,other):
        try: other=other.basis
        except AttributeError: pass
        return self.basis<other
    def __gt__(self,other):
        try: other=other.basis
        except AttributeError: pass
        return self.basis>other
    def __eq__(self,other):
        try: other=other.basis
        except AttributeError: pass
        return self.basis==other
    def __hash__(self):
        return hash(self.basis)
    def __len__(self):
        return len(self.basis)
    def contains(self,other):
        try: diff=self.difference(other,check=False)
        except CommutationError: return False
        return not diff
    @property
    def basis(self):
        return self._basis
    @basis.setter
    def basis(self,basis):
        self._basis=basis.upper()
        self.active=None
    @property
    def results(self):
        try: return self._res
        except AttributeError: raise AttributeError('No results assigned yet')
    @results.setter
    def results(self,counts):
        self._res=counts
    def initialise(self):
        self.active=[i for i,P in enumerate(self.basis) if not P=='I']
    def _apply(self,circuit,barriers=False):
        if barriers: circuit.barrier()
        for j,P in enumerate(self.basis):
            if P in 'IZ': continue
            if P=='Y': circuit.sdg(j)
            if P in 'XY':
                circuit.h(j)
                continue
            raise KeyError(f'Unknown Pauli operator {P}')
        return circuit
    def __str__(self):
        return 'Basis('+str(self.basis)+')'
    def difference(self,other,check=True):
        try: other=other.basis
        except AttributeError: pass
        diff=0
        for p1,p2 in zip_longest(self.basis,other,fillvalue='I'):
            if p2 in (p1,'I'): continue
            if p1=='I':
                diff+=1
                continue
            raise CommutationError('The operators do not commute')
        return diff
