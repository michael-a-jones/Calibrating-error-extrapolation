from backends.interface import get_backendwrapper
from baceknds.counts_backendwrapper import PauliBasis

from qiskit import QuantumCircuit
from qiskit.circuit.parameter import Parameter
from qiskit import transpile

from warnings import warn
from copy import copy

from numpy import pi,average,NaN,exp
from statistics import stdev
from lmfit import Model
from matplotlib import pyplot as plt

#define the function to calculate the energy of Pauli strings
def evaluate_pauli(pauli,results):
    E=0
    for bitstring,weight in results.items():
        sign=1
        for P,B in zip(pauli,bitstring):
            if B=='1' and not P=='I': sign=-sign
        E+=sign*weight
    return E/sum(results.values())
#define the function to use in fitting
def exponential(x,A,b,d):
    return A*exp(b*x)+d

#get the backend
backend=get_backendwrapper('sim1')
#define the number of evaluations to average over
number=2
#define the points to evaluate at
points=[0,pi/2,3*pi/4,pi,3*pi/2]
#define the ideal values of the reference states (or None for the trial state)
reference=[1,1,None,1,1]
#define the circuit repetition depths to evaluate at
reps_range=range(1,51,3)

#define the hamiltonian
Ham={'XX':1,'YY':1,'ZZ':1}
#define the measurement bases
bases=tuple(PauliBasis(b) for b in Ham)
#define the high noise limit
noise_limit=0

#create a dictionary to hold the data
data={}
for point in points:
    data[point]=[]
    for circuit_reps in reps_range:
      
        #generate the circuit
        circ=QuantumCircuit(2,2)
        circ.x(0)
        for i in range(circuit_reps):
            circ.h(0)
            circ.cx(0,1)
            circ.ry(point/circuit_reps,0)
            circ.ry(point/circuit_reps,1)
            circ.barrier()
            circ.cx(0,1)
            circ.h(0)
            circ.barrier()
            
        es=[]
        for i in range(number):
            #run the circuit
            counts=backend.run(circ,bases)
            e=0
            #evaluate each term in the Hamiltonian
            for pauli in Ham:
                #check each measurement to see if it measured that term
                for c in counts:
                    if c.contains(pauli):
                        #evaluate the term
                        e+=evaluate_pauli(pauli,c.results)
                        continue
            es.append(e)
        #store the data
        data[point].append([average(es),stdev(es)])
       
#make the fits
fits={}
x=list(reps_range)
for point in points:
    y=[av for av,dev in data[point]]
    model=Model(exponential)
    p=model.make_params(A=y[0],b=-1,d=noise_limit)
    try: result=model.fit(y,p,x=x)
    except ValueError:
        fits[point]=[NaN,NaN,NaN]
        warn(f'Fitting failed for {point}')
    else: fits[point]=[p.value for s,p in result.params.items()]


#plot the data
for point in points:
    y,y_err=zip(*data[point])
    fig_trial=plt.figure()
    ax=fig_trial.add_subplot(111)
    ax.errorbar(x,y,y_err,fmt='bx',ecolor='k')
    X=[0]+list(circuit_reps_range)
    Y=[exponential(x,*fits[point]) for x in X]
    ax.plot(X,Y,'k-')
    ax.set_xlabel('Noise level')
    ax.set_ylabel('Energy')
    ax.set_title(f'Theta={point}')
