# STEP3: Vibronic Hamiltonian Construction

This step builds the vibronic Hamiltonians and non-adiabatic couplings (NACs) needed for the NAMD simulations in STEP4. It processes the output of STEP2 (Kohn–Sham orbital energies, TDDFT excitations) and produces the files used by the surface-hopping dynamics.

---

## Directory Structure

STEP3/
├── README.md
└── adamantane/
    ├── adamantane_step3.ipynb
    ├── step33three_B3LYP/
    │   ├── Hvib_ci_*_re.npz     # Real part of vibronic Hamiltonian (CI basis)
    │   ├── Hvib_ci_*_im.npz     # Imaginary part (NACs)
    │   ├── Hvib_sd_*_re.npz     # Real part (SD basis)
    │   ├── Hvib_sd_*_im.npz     # Imaginary part (SD basis)
    │   ├── St_ci_*_re.npz       # State energies (CI basis)
    │   ├── St_sd_*_re.npz       # State energies (SD basis)
    │   └── T_*.txt              # Additional time-dependent data
    └── step33three_B3LYP.tar.gz

---

## Workflow

The notebook adamantane_step3.ipynb performs the following tasks:

1. Active-Space Selection
   - Reads the Kohn–Sham orbital energies from STEP2.
   - Selects a reduced active space around the HOMO/LUMO.

   For adamantane, the active space is:
   - Orbitals: 18–46
   - HOMO: orbital 38
   - LUMO: orbital 39

2. Data Generation (step3.run_step3_sd_nacs_libint)
   - Computes the vibronic Hamiltonian (Hvib_ci_*) and NACs using the selected active space.
   - Produces both real and imaginary parts of the Hamiltonian.

3. Visualization of NACs
   - Plots the average absolute value of the NAC matrix as a heatmap.

4. Influence Spectrum
   - Computes the influence spectrum from the S1–S0 energy gap.
   - Quantifies the contribution of low-frequency modes.

5. Time-Dependent Excitation Energies
   - Extracts and plots the time-dependent excitation energies relative to the ground state.

6. Kohn–Sham Orbital-Energy Visualization
   - Visualizes the time evolution of Kohn–Sham orbital energies for a short trajectory segment.

7. S1–S0 Energy-Gap Distribution
   - Computes the probability density distribution of the S1–S0 gap.

8. NAC Distribution
   - Computes the probability density distribution of |NAC(S0, S1)|.

9. Distributions for All Structures
   - Runs the same analysis for all four structures and compares the results in a single figure.

---

## Key Outputs

| Output | Description |
|--------|-------------|
| Hvib_ci_*_re.npz | Real part of the vibronic Hamiltonian (CI basis) |
| Hvib_ci_*_im.npz | Imaginary part = non-adiabatic couplings (CI basis) |
| Hvib_sd_*_re.npz | Real part of the vibronic Hamiltonian (SD basis) |
| Hvib_sd_*_im.npz | Imaginary part (SD basis) |
| St_ci_*_re.npz | State energies (CI basis) |
| St_sd_*_re.npz | State energies (SD basis) |
| T_*.txt | Additional time-dependent data |

These files are the direct input to the NAMD simulations in STEP4.

---

## How to Run

Open the notebook and run all cells sequentially:

jupyter notebook adamantane_step3.ipynb

The notebook is fully self-contained and generates all the required files for STEP4.

---

## Contact

- Qingxin Zhang (qingxin@buffalo.edu)
- Hamid Zabihi Hashjin (h.zabihi@eng.uk.ac.ir)
- Alexey V. Akimov (alexeyak@buffalo.edu)