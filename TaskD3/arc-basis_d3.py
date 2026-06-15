#ARC Program functions
#This cell generates the Hamiltonian, solves it, and exports the data

from arc import *
import numpy as np
import sys
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import os  # Import the 'os' module to handle file paths

def state_to_string(n,l,j,m):
    if l == 0:
        lstring = "s"
    if l == 1:
        lstring = "p"
    if l == 2:
        lstring = "d"
    return str(int(n)) + "$" + lstring + "_{" + str(j) + "},$ "+ "m=" + str(m)

def arc_calculate(a1,n1, a2,n2):
    if a1 == "Cs":
        atom1 = Caesium()
    elif a1 == "Rb":
        atom1 = Rubidium87()
    # n1 = 61
    l1 = 1
    j1 = 1.5
    m1 = 0.5
    if a2 == "Cs":
        atom2 = Caesium()
    elif a2 == "Rb":
        atom2 = Rubidium87()
    # n2 = 61
    l2 = 1
    j2 = 1.5
    m2 = 0.5

    Bz = 0.0001# #Tesla

    rmin = 1.5
    rmax = 6

    theta = 0  # Polar Angle [0-pi]
    phi = 0 # Azimuthal Angle [0-2pi]
    dn = 4# Range of n to consider (n0-dn:n0+dn)
    dl = 4 # Range of l values
    deltaMax = 25e9  # Max pair-state energy difference [Hz]

    nEig = 250  # Number of eigenstates to extract


    calc = PairStateInteractions(atom1, n1, l1, j1, n2, l2, j2, m1, m2, s=0.5, s2=0.5, atom2=atom2)
    r = np.linspace(rmin, rmax, 100)
    # Generate pair-state interaction Hamiltonian
    calc.defineBasis(theta, phi, dn, dl, deltaMax, Bz = Bz, progressOutput=True, debugOutput=False)
    # Diagonalize the Hamiltonian
    calc.diagonalise(r, nEig, progressOutput=True)
    # Plota
    calc.plotLevelDiagram()
    calc.ax.set_xlim(rmin, rmax)
    calc.ax.set_ylim(-0.06, 0.06)

    # --- PROBLEM 1 FIX: Commented out interactive plot command ---
    # calc.showPlot() # This command will hang the script in a batch job.

    # --- PROBLEM 2 FIX: Use correct path separator and ensure directory exists ---
    output_dir = "saved-calcs_d3"
    os.makedirs(output_dir, exist_ok=True) # Ensure the directory is created
    # Use forward slash '/' for Linux paths
    calc.exportData(f"{output_dir}/{n1}{n2}", exportFormat='csv')

    separations_file = f"{output_dir}/{n1}{n2}_r.csv"
    datapoints_file = f"{output_dir}/{n1}{n2}_energyLevels.csv"
    highlights_file = f"{output_dir}/{n1}{n2}_highlight.csv"

    separations = np.genfromtxt(separations_file, delimiter=',')
    datapoints = np.genfromtxt(datapoints_file, delimiter=',') #* 1000 #Units are now in MHz
    highlights = np.genfromtxt(highlights_file, delimiter=',')
    # print(datapoints[0].shape)

    highest_overlap_points = []

    fig, ax = plt.subplots(figsize=(8, 6))
    norm = mcolors.LogNorm(vmin=0.001, vmax=1)
    colors = ["whitesmoke", "orange", "darkred"]
    cmap = mcolors.LinearSegmentedColormap.from_list("custom_colormap", colors, N=256)

    for i in range(len(separations)):#for each r
        xplot = np.ones(len(datapoints[i])) * separations[i]
        #Sort the data so that the most overlap data is on top.
        sorted_pairs = sorted(zip(highlights[i], datapoints[i]))  # Sort by highlights
        sorted_highlights, sorted_ydata = zip(*sorted_pairs)  # Unzip into separate lists
        # Convert back to lists
        sorted_ydata = list(sorted_ydata)
        sorted_highlights = list(sorted_highlights)
        # print(xplot)
        # print(sorted_ydata)
        highest_overlap_points.append(sorted_ydata[-1])  # Appends the strongest eigenvalue for each x point
        plt.scatter(xplot, sorted_ydata,c=sorted_highlights, cmap=cmap, norm = norm, s=8)

    plt.ylim((-0.25,0.4))
    plt.xlim(np.min(separations),np.max(separations))
    plt.axhline(y=0.0, color='black', linestyle='-')
    cbar = plt.colorbar()
    cbar.set_label('Overlap', fontsize = 16)
    plt.xlabel("$r$ Inter-Nuclear Distance ($\\mu$m)",fontsize=20)
    plt.ylabel("Pair State Energy (GHz)",fontsize=20)
    a1label = a1 +": " + state_to_string(n1,l1,j1,m1)
    a2label = a2 +": " + state_to_string(n2,l2,j2,m2)
    calc_data = str(theta * 180 / np.pi) + "$\\degree,$ $B_z = $ " + str(Bz * 1e4) + " G"
    plt.text(4,0.12,a1label)
    plt.text(4,0.1,a2label)
    plt.text(4,0.08,calc_data)
    
    # --- FIX: Save the plot to a file instead of showing it interactively ---
    plt.savefig(f"{output_dir}/{n1}{n2}.png")
    # plt.show() # This command will hang the script in a batch job.

arc_calculate("Cs",eval(sys.argv[1]),"Rb",eval(sys.argv[2]))
#for i in range(40,81):
    #for j in range(40,81):
        #arc_calculate("Cs", i, "Rb", j)
#arc_calculate("Cs", 40, "Rb", 40)