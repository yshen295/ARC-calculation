# Cs–Rb Rydberg Pair-State Interaction Calculations

This repository contains Python scripts for calculating and visualizing Rydberg pair-state interaction spectra for heteronuclear cesium–rubidium atom pairs using the [ARC: Alkali Rydberg Calculator](https://arc-alkali-rydberg-calculator.readthedocs.io/).

The scripts build the pair-state interaction Hamiltonian for selected Cs and Rb Rydberg states, diagonalize the Hamiltonian over a range of interatomic separations, export numerical data, and generate level-diagram plots colored by overlap with the target pair state.

## Project Overview

The calculation focuses on pair-state interactions of the form

```text
Cs(n1 l j m) + Rb(n2 l j m)
```

where `n1` and `n2` are supplied from the command line. Each script fixes the angular momentum state and geometric configuration, then scans the interaction energy as a function of interatomic distance.

For each calculation, the code:

1. Defines Cs and Rb atom objects using ARC.
2. Constructs a pair-state basis with nearby Rydberg states.
3. Diagonalizes the pair-state interaction Hamiltonian.
4. Exports the calculated energy levels and overlap data as CSV files.
5. Generates a plot of pair-state energy versus interatomic distance.

## Files

| File              | Target state  | Magnetic sublevel | Angle `theta` | Output folder     |
| ----------------- | ------------- | ----------------: | ------------: | ----------------- |
| `arc-basis_d1.py` | Cs/Rb `nS1/2` |         `m = 1/2` |         `π/2` | `saved-calcs_d1/` |
| `arc-basis_d2.py` | Cs/Rb `nP3/2` |         `m = 3/2` |         `π/2` | `saved-calcs_d2/` |
| `arc-basis_d3.py` | Cs/Rb `nP3/2` |         `m = 1/2` |           `0` | `saved-calcs_d3/` |
| `arc-basis_d4.py` | Cs/Rb `nP3/2` |         `m = 1/2` |         `π/2` | `saved-calcs_d4/` |
| `arc-basis_d5.py` | Cs/Rb `nP3/2` |         `m = 3/2` |           `0` | `saved-calcs_d5/` |
| `arc-basis_d6.py` | Cs/Rb `nS1/2` |         `m = 1/2` |           `0` | `saved-calcs_d6/` |

## Physical and Numerical Parameters

The scripts use the following shared parameters:

| Parameter  |          Value | Meaning                                                         |
| ---------- | -------------: | --------------------------------------------------------------- |
| `atom1`    |    `Caesium()` | First atom species                                              |
| `atom2`    | `Rubidium87()` | Second atom species                                             |
| `Bz`       |     `0.0001 T` | Magnetic field along the quantization axis, equal to 1 G        |
| `rmin`     |       `1.5 μm` | Minimum interatomic separation                                  |
| `rmax`     |         `6 μm` | Maximum interatomic separation                                  |
| `r`        |     100 points | Distance grid from `rmin` to `rmax`                             |
| `phi`      |            `0` | Azimuthal angle                                                 |
| `dn`       |            `4` | Range of nearby principal quantum numbers included in the basis |
| `dl`       |            `4` | Range of orbital angular momentum states included in the basis  |
| `deltaMax` |      `25e9 Hz` | Maximum pair-state energy defect included in the basis          |
| `nEig`     |          `250` | Number of eigenstates extracted at each separation              |

## Requirements

Install the required Python packages:

```bash
pip install numpy scipy matplotlib ARC-Alkali-Rydberg-Calculator
```

Depending on your environment, the ARC package may also be installed as:

```bash
pip install arc-alkali-rydberg-calculator
```

## Usage

Each script takes two command-line arguments:

```bash
python arc-basis_d1.py <n_Cs> <n_Rb>
```

For example:

```bash
python arc-basis_d1.py 50 55
```

This calculates the interaction spectrum for

```text
Cs(50S1/2, m=1/2) + Rb(55S1/2, m=1/2)
```

using the geometry and parameters specified in `arc-basis_d1.py`.

Another example:

```bash
python arc-basis_d2.py 50 55
```

This calculates

```text
Cs(50P3/2, m=3/2) + Rb(55P3/2, m=3/2)
```

with `theta = π/2`.

## Output

Each script creates a corresponding output directory, such as:

```text
saved-calcs_d1/
saved-calcs_d2/
...
saved-calcs_d6/
```

For each pair of input quantum numbers, the generated files include:

| Output file                 | Description                                                      |
| --------------------------- | ---------------------------------------------------------------- |
| `{n1}{n2}_r.csv`            | Interatomic separations in micrometers                           |
| `{n1}{n2}_energyLevels.csv` | Pair-state eigenenergies as a function of separation             |
| `{n1}{n2}_highlight.csv`    | Overlap of each eigenstate with the target pair state            |
| `{n1}{n2}.png`              | Level diagram showing energy versus distance, colored by overlap |

For example, running

```bash
python arc-basis_d6.py 50 55
```

will generate files such as:

```text
saved-calcs_d6/5055_r.csv
saved-calcs_d6/5055_energyLevels.csv
saved-calcs_d6/5055_highlight.csv
saved-calcs_d6/5055.png
```

## Plot Description

The generated plots show:

* x-axis: interatomic separation `r` in micrometers
* y-axis: pair-state energy in GHz
* color scale: overlap with the selected target pair state
* darker/redder points: stronger overlap with the chosen Cs–Rb Rydberg pair state

The overlap coloring is useful for identifying which eigenenergy branches are most relevant to the initial target pair state.

## Notes

The current scripts are designed for batch-style calculations. Interactive plotting commands such as `plt.show()` and `calc.showPlot()` are disabled or omitted so that the scripts can run on remote servers or computing clusters without hanging.

The scripts currently use:

```python
eval(sys.argv[1])
eval(sys.argv[2])
```

For safety, this can be replaced with:

```python
int(sys.argv[1])
int(sys.argv[2])
```

because the command-line inputs are expected to be integer principal quantum numbers.

## Scientific Context

These calculations are useful for studying Rydberg–Rydberg interaction potentials, Förster resonances, blockade shifts, and state mixing in heteronuclear Cs–Rb Rydberg atom pairs. The exported overlap and energy data can be used to identify relevant pair-state channels for neutral-atom quantum information and Rydberg blockade gate studies.

## Repository Structure

A typical repository layout is:

```text
.
├── arc-basis_d1.py
├── arc-basis_d2.py
├── arc-basis_d3.py
├── arc-basis_d4.py
├── arc-basis_d5.py
├── arc-basis_d6.py
├── README.md
└── saved-calcs_*/        # generated output folders, usually not tracked by git
```

It is recommended to add generated output folders to `.gitignore` if the CSV and PNG files are large:

```text
saved-calcs_d*/
*.png
*.csv
```

Alternatively, keep selected plots or processed data in the repository if they are part of the final analysis.
