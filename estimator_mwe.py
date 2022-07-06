from qiskit import QuantumCircuit
from qiskit.circuit.parameter import Parameter
from backends.interface import get_backendwrapper
from backends.counts_backendwrapper import PauliBasis

def evaluate_pauli(pauli,results):
    E=0
    for bitstring,weight in results.items():
        sign=1
        for P,B in zip(pauli,bitstring):
            if B=='1' and not P=='I': sign=-sign
        E+=sign*weight
    return E/sum(results.values())

backend=get_backendwrapper('sim0')
theta=0

theta_=Parameter('t')
circ=QuantumCircuit(2,2)
circ.x(0)
circ.h(0)
circ.cx(0,1)
circ.ry(theta_,0)
circ.ry(theta_,1)
circ.cx(0,1)
circ.h(0)

Ham={'XX':1,'YY':1,'ZZ':1}
bases=tuple(PauliBasis(b) for b in Ham)

counts=backend.run(circ.assign_parameters({theta_:theta}),bases)

E=0
for pauli in Ham:
    for c in counts:
        if c.contains(pauli):
            E+=evaluate_pauli(pauli,c.results)
            continue
print(E)
