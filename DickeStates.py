## Explicit Method

# Imports
from qiskit import IBMQ
from qiskit.compiler import transpile
from qiskit import QuantumCircuit, ClassicalRegister, QuantumRegister, Aer, execute
from qiskit.tools.visualization import plot_histogram
import math
import time
import matplotlib.pyplot as plt
import numpy as np
from fractions import Fraction

# Variables
pi = np.pi

# Functions
def ccr(alpha,theta,c1,c2,u,qc,n): # control-control-rotation gate gate
    
    pi = np.pi
    sim = QuantumRegister(n,'sim')
    # alpha is x,y,z. u is qubit that u acts on. c_1,c_2 are the control locations
    if alpha == 'x':
        qc.cu3(theta/2,-pi/2,pi/2,sim[c2],sim[u])
        qc.cx(sim[c1],sim[c2])
        qc.cu3(-theta/2,-pi/2,pi/2,sim[c2],sim[u])
        qc.cx(sim[c1],sim[c2])
        qc.cu3(theta/2,-pi/2,pi/2,sim[c1],sim[u])
    elif alpha == 'y':
        qc.cu3(theta/2,0,0,sim[c2],sim[u])
        qc.cx(sim[c1],sim[c2])
        qc.cu3(-theta/2,0,0,sim[c2],sim[u])
        qc.cx(sim[c1],sim[c2])
        qc.cu3(theta/2,0,0,sim[c1],sim[u])
    else:
        qc.cu1(theta/2,sim[c2],sim[u])
        qc.cx(sim[c1],sim[c2])
        qc.cu1(-theta/2,sim[c2],sim[u])
        qc.cx(sim[c1],sim[c2])
        qc.cu1(theta/2,sim[c1],sim[u])
    
def scs(x,y,qc,n): # s is starting qubit, qc is quantum circuit
    sim = QuantumRegister(n,'sim')
    for i in range(1,y+1):
        if i == 1:
            qc.cx(sim[x-1-i],sim[x-1])
            qc.cu3(2*math.acos(math.sqrt(i/x)),0,0,sim[x-1],sim[x-1-i])
            qc.cx(sim[x-1-i],sim[x-1])
        else:
            qc.cx(sim[x-1-i],sim[x-1])
            ccr('x',2*math.acos(math.sqrt(i/x)),x-1,x-i,x-1-i,qc,n)
            qc.cx(sim[x-1-i],sim[x-1])

# Returns Dicke circuit
def Dicke_exp(n,k):
    
    # Create circuit
    sim = QuantumRegister(n,'sim')
    meas = ClassicalRegister(n,'meas')
    circ = QuantumCircuit(sim, meas)
    
    # Prepare initial bit string
    for i in range(n-1,n-k-1,-1):
        circ.x(sim[i])
    
    for l in range(n,k,-1):
        scs(l,k,circ,n)
    for l in range(k,1,-1):
        scs(l,l-1,circ,n)
        
    # Measure
    circ.measure(sim,meas)
    
    trans = transpile(circ, basis_gates=['u1', 'u2', 'u3', 'cx'], optimization_level=3)

    # Draw
    return circ