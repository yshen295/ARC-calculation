from arc import *
import numpy as np
import sys
from scipy.sparse.linalg import eigsh
from arc.alkali_atom_functions import printStateString
import os

def hamiltonian_at_R(calc, R_um):

    R_m = R_um * 1e-6
    m = calc.matDiagonal.copy()
    rX = R_m**3
    for matRX in calc.matR:
        m += matRX / rX
        rX *= R_m
    return m

def arc_eigenstate(a1, n1, a2, n2):
    atom1 = Caesium() if a1 == "Cs" else Rubidium87()
    atom2 = Caesium() if a2 == "Cs" else Rubidium87()
    previousEigenvectors = []

    
    l1, j1, m1 = 1, 1.5, 1.5
    l2, j2, m2 = 1, 1.5, 1.5
    Bz = 0.0001  # Tesla
    rmin, rmax = 1.5, 6
    theta, phi = 0, 0
    dn = 7
    dl = 7
    deltaMax = 25e9
    nEig = 250

    calc = PairStateInteractions(atom1, n1, l1, j1, n2, l2, j2, m1, m2,
                                 s=0.5, s2=0.5, atom2=atom2)
    r = np.linspace(rmin, rmax, 500)
    calc.defineBasis(theta, phi, dn, dl, deltaMax, Bz=Bz,
                     progressOutput=True, debugOutput=False)

    m = hamiltonian_at_R(calc, 2.5)

    eigenstateDetuning = 0
    ev, egvector = eigsh(m, k=nEig, sigma=eigenstateDetuning * 1.0e-9, which="LM", tol=1e-6)

    output_dir = "saved-eigenstates"
    os.makedirs(output_dir, exist_ok=True)

    # construct basis label
    basis_labels = []
    for bs in calc.basisStates:
        state1 = printStateString(bs[0], bs[1], bs[2], bs[3])
        state2 = printStateString(bs[4], bs[5], bs[6], bs[7])
        basis_labels.append(f"|{state1}; {state2}>")


    output_path = os.path.join(output_dir, f"eigenpairs_{n1}_{n2}.csv")
    with open(output_path, "w", encoding="utf-8") as f:
        f.write("Eigenvalue," + ",".join(basis_labels) + "\n")
        for k in range(len(ev)):
            row = [str(ev[k])] + [str(egvector[i, k]) for i in range(egvector.shape[0])]
            f.write(",".join(row) + "\n")

    print(f"saved to {output_path}")

arc_eigenstate("Cs",eval(sys.argv[1]),"Rb",eval(sys.argv[2]))
