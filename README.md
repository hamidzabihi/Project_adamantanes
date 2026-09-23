# Project Adamantanes

A computational study of the electronic structure and excited-state relaxation dynamics of pristine and methylated adamantanes.

---

## Overview

This project — *Nonadiabatic Molecular Dynamics in Pristine and Methylated Adamantanes* — combines static electronic-structure calculations with non-adiabatic molecular dynamics (NAMD) simulations to investigate the photophysical properties of adamantane and its methylated derivatives.

The repository contains the complete computational workflow for:

- Quantum-chemistry calculations of the electronic structure of adamantane and its methylated derivatives.
- Non-adiabatic molecular dynamics (NAMD) simulations of excited-state relaxation.
- Analysis of radiative and non-radiative recombination, quantum yields, and their dependence on molecular structure.

In addition to the dynamic trajectory method, a single-optimization approach was also employed.

---

## Repository Structure

Project_adamantanes/
├── LICENSE
├── README.md                          # This file
├── NAMD/                              # NAMD workflow (STEP1 → STEP4)
│   ├── STEP1/                         # Structure preparation
│   ├── STEP2/                         # Electronic-structure calculations
│   ├── STEP3/                         # Vibronic Hamiltonians and NACs
│   └── STEP4/                         # NAMD simulations and analysis
└── Static_Optimization_Structure/     # Optimized geometries and IP/UV-Vis data

---

## Workflow Overview

### Static Optimization Structure

Optimize the initial crystal structure using the **B3LYP** functional.

- **Input:** Optimized CIF files, ORCA and CP2K input files with the appropriate `energy_type` setting (4 different systems).
- **Output:** Ground state and TDDFT output files and corresponding Molden files. The TDDFT output files can be further analyzed with Multiwfn to generate UV–Vis spectra in `.txt` format. Different spectral line shapes and degrees of peak broadening can be obtained by adjusting the full width at half maximum (FWHM).

Related scripts and examples are available in:

Static_Optimization_Structure/UV_Vis_spectra

---

### Nonadiabatic Molecular Dynamics Workflow

#### STEP 1: Adiabatic Molecular Dynamics

- **Objective:** Obtain equilibrated adiabatic trajectories at 300 K.
- **Output:** Equilibrated ground-state molecular-dynamics trajectories for subsequent electronic-structure calculations.

**Additional Note: Restarting a CP2K MD Simulation:**

To continue a previously interrupted CP2K molecular-dynamics simulation, include the following restart settings in the input file:

&EXT_RESTART
  RESTART_FILE_NAME adamantane_MD-1_5000.restart
  RESTART_DEFAULT T
&END EXT_RESTART

---

#### STEP 2: TDDFT Calculations

- **Input:** TDDFT calculations based on a 300 K temperature profile.
- **Output:** Excited-state energies, oscillator strengths, and Kohn–Sham orbital energies for visualizing electronic entropy and UV-Vis spectra.

---

#### STEP 3: Vibronic Hamiltonian Construction

- **Objective:** Generate the time-dependent vibronic Hamiltonian using the Step 2 outputs and the selected active space.
- **Active space:** Taking adamantane as an example, orbitals 18–46, with orbital 38 as the HOMO and orbital 39 as the LUMO.
- **Output:** Time-dependent electronic energies, probability distributions, and non-adiabatic couplings for the selected electronic states.

---

#### STEP 4: Nonadiabatic Molecular Dynamics

- **Objective:** Perform non-adiabatic molecular-dynamics simulations using the Step 3 Hamiltonian data with four different surface-hopping methods: **FSSH, FSSH2, IDA, and MSDM**.
- **Output:** Determine the timescale using relaxation dynamics, and compute PLQY, non-radiative lifetimes, and stretch exponents.

---

## Workflow Summary Table

| Step | Purpose | Key Outputs |
|------|---------|-------------|
| Static Optimization | B3LYP geometry optimization and UV-Vis | CIF, Molden, TDDFT output, spectra |
| STEP1 | Prepare and align molecular structures, run MD at 300 K | aligned.tar.gz, MD trajectories |
| STEP2 | TDDFT calculations along MD trajectories | Kohn–Sham energies, UV-Vis spectra |
| STEP3 | Build vibronic Hamiltonians and non-adiabatic couplings | Hvib_ci_*.npz, St_ci_*.npz |
| STEP4 | Run NAMD simulations and analyze kinetics | mem_data.hdf, analysis_results.pkl |

---

## Structures Studied

- Adamantane
- 1-Methyl-adamantane
- 1,3-Dimethyl-adamantane
- 1,3,5-Trimethyl-adamantane

---

## Surface Hopping Methods

The NAMD simulations in STEP4 employ four different trajectory surface-hopping (TSH) methods:

| Method | Description |
|--------|-------------|
| FSSH | Fewest Switches Surface Hopping (Tully, 1990) |
| FSSH2 | Alternative FSSH formulation with robust switching |
| IDA | Instantaneous Decoherence at Attempted hops |
| MSDM | Simplified Decay of Mixing |

---

## Dependencies

- Python 3.10+
- numpy, scipy, h5py, matplotlib
- libra_py (Libra library)
- CP2K, ORCA, Multiwfn (for STEP2 / static calculations)

---

## Contact

- Qingxin Zhang (qingxin@buffalo.edu)
- Hamid Zabihi Hashjin (h.zabihi@eng.uk.ac.ir)
- Alexey V. Akimov (alexeyak@buffalo.edu)