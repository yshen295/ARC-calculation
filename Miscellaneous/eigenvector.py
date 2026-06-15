#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  1 16:24:04 2025

@author: yeddie
"""
from arc import *
import numpy as np
import sys
from scipy.optimize import curve_fit
from scipy.sparse.linalg import eigsh
from arc.alkali_atom_functions import printStateString
%matplotlib inline
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import os


def arc_eigenstate(a1,n1,a2,n2):
    if a1 == "Cs":
        atom1 = Caesium()
    elif a1 == "Rb":
        atom1 = Rubidium87()
    # n1 = 61
    l1 = 1
    j1 = 1.5
    m1 = 1.5
    if a2 == "Cs":
        atom2 = Caesium()
    elif a2 == "Rb":
        atom2 = Rubidium87()
    # n2 = 61
    l2 = 1
    j2 = 1.5
    m2 = 1.5
    

    Bz = 0.0001# #Tesla
    
    rmin = 1.5
    rmax = 6
    
    theta = 0 # Polar Angle [0-pi]
    phi = 0 # Azimuthal Angle [0-2pi]
    dn = 3 # Range of n to consider (n0-dn:n0+dn)
    dl = 3 # Range of l values
    deltaMax = 25e9  # Max pair-state energy difference [Hz]
    
    nEig = 250  # Number of eigenstates to extract
    
    
    calc = PairStateInteractions(atom1, n1, l1, j1, n2, l2, j2, m1, m2, s=0.5, 
                                 s2=0.5, atom2=atom2)
    
    r = np.linspace(rmin, rmax, 500)

    calc.defineBasis(theta, phi, dn, dl, deltaMax, Bz = Bz, 
                     progressOutput=True, debugOutput=False)
    
    calc.diagonalise(r, nEig, sortEigenvectors=False, progressOutput=True)

    i = np.argmin(np.abs(calc.r - 2.5))  # nearest to R = 2.5 μm
    
    overlap = np.asarray(calc.highlight[i]) 

    k = 10
    
    idx = np.argpartition(overlap, -k)[-k:]
    
    top_eigvals = calc.y[i][idx]
    
    print("Top-k index (unsorted):", idx)
    print("Top-k values:", overlap[idx])
    
    mask = np.zeros_like(overlap, dtype=bool)
    mask[idx] = True

    egvector = np.array(calc.composition)
    
    dim = egvector.shape[1]
    
    print(dim)
    
    print(egvector.shape)
    
    with open("output.txt", "w", encoding="utf-8") as f:
        for j in range(dim):
            if j in idx:
                f.write(str(egvector[i][j]) + "\n")
    # R = 2.5  # Target distance  


    return

arc_eigenstate("Cs", 50, "Rb", 55)