# quantum-computed-chemistry/backends/statevector_backendwrapper.py
# A class to wrap backends that return statevectors (as opposed to counts)

#qiskit
from qiskit.compiler import transpile

#qcchem
from .base_backendwrapper import BackendWrapper

class StatevectorBackend(BackendWrapper):
    """
    A class to wrap backends that return statevectors (as opposed to measurement counts)
        main methods:
            StatevectorBackend.run
            StatevectorBackend.circ2circuits
            StatevectorBackend.circuits2jobs
            StatevectorBackend.jobs2vector
    """
    def circ2circuits(self,circuit,bases=None):
        """
        Prepare a circuit to be submitted
        arguments:
            circuit (qcchem.Circuit or qiskit.QuantumCircuit): the state preparation circuit to be measured
            bases (NA): NB. this argument exists to provide consistency between CountsBackend and StatevectorBackend
        returns:
            (circuit,)
        """
        #qubits and meas_basis are not required but are included so the call signature matches CountsBackend.circ2circuits
        return (transpile(circuit,self),)
    def jobs2vector(self,jobs):
        """
        Retrieve the results for a previously submitted batch of jobs
        arguments:
            jobs (iterable of qiskit jobs and/or str): the jobs for which results should be retrieved, if given as strings, they are assumed to be job IDs for jobs run on this device
        returns:
            statevector (qiskit statevector)
        """
        assert len(jobs)==1
        if isinstance(jobs[0],str): jobs=[self.retrieve_job(job) for job in jobs]
        return jobs[0].result().get_statevector()
    def jobs2counts(self,jobs,basis=None):
        """
        alias to jobs2vector
            NB. this method exists to provide consistency between CountsBackend and StatevectorBackend
        """
        return self.jobs2vector(jobs)
    def run(self,circuit,bases=None,shots=None):
        """
        Convert a state preparation circuit to counts in a given set of measurement bases
        see also:
            StatevectorBackend.circ2circuits
            StatevectorBackend.circuits2jobs
            StatevectorBackend.jobs2vector
        arguments:
            circuit (qcchem.Circuit or qiskit.QuantumCircuit): the state preparation circuit to be measured
            bases (NA): NB. this argument exists to provide consistency between CountsBackend and StatevectorBackend
            shots (NA): NB. this argument exists to provide consistency between CountsBackend and StatevectorBackend
        returns:
            statevector (qiskit statevector)
        """
        #qubits, meas_basis and measurment_method are not required but are included so the call signature matches CountsBackend.run
        return self.jobs2vector(self.circuits2jobs(self.circ2circuits(circuit)))
