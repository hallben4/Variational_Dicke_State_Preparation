import math
import cmath
get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib as mp
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from numpy import linalg as la
from scipy import linalg as sla
from IPython.display import clear_output
import time
import pickle
from qiskit.extensions import XGate, UnitaryGate
from qiskit import *
from qiskit.compiler import transpile
from qiskit.visualization import plot_histogram

provider = IBMQ.load_account()
from qiskit.providers.aer.noise import NoiseModel

from qiskit.quantum_info.synthesis import two_qubit_cnot_decompose as two_qubit_decomp

X = np.matrix([[0,1],[1,0]],dtype=complex)
Y = np.matrix([[0,-1j],[1j,0]],dtype=complex)

swap = np.array([[1,0,0,0],
                 [0,0,1,0],
                 [0,1,0,0],
                 [0,0,0,1]])

## pswap
# Partial swap matrix
def pswap(theta):
    
    mat = np.array([[1,0,0,0],[0,np.cos(theta),-np.sin(theta),0],[0,np.sin(theta),np.cos(theta),0],[0,0,0,1]])
    
    return mat


## DickeCirc
# Creates variational quantum circuit to prepare the Dicke state D(n,k)
# n is the number of qubits
# k is the Hamming weight
# theta is list of variational parameters

def DickeCirc(n,k,layer,theta,init='pure'):
    
    k_bad = False
    if k > int(n/2):
        k = n - k
        k_bad = True
    
    # Initialize quantum circuit
    q = QuantumRegister(n,'q')
    c = ClassicalRegister(n,'c')
    circ = QuantumCircuit(q,c)
    
    # Initialize list for indeces of qubits that will have X applied to them
    index_list = []
    for i in range(k):
        index_list.append(2*i)
    
    # Non-variational Section: non-variationally prepares some 
    # superpostion of quantum states, each with
    # Hamming weight k
    
    # Option 1: Pure state (apply X-gates)   
    if init == 'pure':
        for i in index_list:
            circ.x(q[i])
        
    # Option 2: Mixed state
    # Create Bell |Psi^+> and copy
    if init == 'mixed':
        circ.h(q[0])
        circ.x(q[1])
        for i in range(int(n/2)-1):    
            circ.cx(q[0],q[2*i+2])
            circ.cx(q[1],q[2*i+3])
        circ.barrier()

    # Variational Section:

    # Initialize theta_index
    theta_index = 0

    for l in range(layer):
        
        # First layer (l=0)
        if l != 0:
            # Subsequent layers (l>0)
            new_index_list = []
            for i in index_list:
                new_index_list.append(i-1)
            new_index_list.append(i+1)

            index_list = new_index_list
        
        for elm in index_list:
            if elm >= 0 and elm < n-1:
                circ.append(two_qubit_decomp(pswap(theta[theta_index])),[q[int(elm)],q[int(elm+1)]])
                theta_index += 1 
             
    # k too large case 
    if k_bad:
        circ.x(q)
        
    return circ

