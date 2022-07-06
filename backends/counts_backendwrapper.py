# quantum-computed-chemistry/backends/counts_backendwrapper.py
# A class to wrap backends that return measurement counts (as opposed to a statevector)

#core
from time import sleep

#qiskit
from qiskit import transpile,ClassicalRegister
from qiskit.ignis.mitigation.measurement import complete_meas_cal,CompleteMeasFitter

#qcchem
from .base_backendwrapper import BackendWrapper
from ..measurement.base_measurement import BaseMeasurement
from ..measurement.pauli_measurement import PauliBasis


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
    
    def calibrate_pairwise_qrem(self,qubits=None):
        if qubits is None: qubits=self.n_qubits
        if qubits is None: raise ValueError('The number of qubits must be specified for a simulator')
        if isinstance(qubits,int): qubits=range(qubits)
        self.filters=[]
        for i,q1 in enumerate(qubits):
            self.filters.append([])
            for q2 in qubits[i+1:]:
                cal_circuits,states=complete_meas_cal([q1,q2])
                res=self.circuits2jobs(cal_circuits)[0].result()
                fitter=CompleteMeasFitter(res,states)
                qrem_filter=fitter.filter
                self.filters[-1].append(qrem_filter)
        self.filters.pop(-1)
    def apply_pairwise_qrem(self,counts,qubits):
        assert len(qubits)==2
        assert qubits[0]<qubits[1]
        reduced_counts={}
        if isinstance(counts,BaseMeasurement): counts=counts.results
        for key,val in counts.items():
            key=key[qubits[0]]+key[qubits[1]]
            reduced_counts[key]=reduced_counts.get(key,0)+val
        return self.filters[qubits[0]][qubits[1]].apply(reduced_counts)
        
#     def circ2counts(self,circuit,bases,shots=8192):
#         return self.jobs2counts(self.circuits2jobs(self.circ2circuits(circuit,bases),shots),bases)
