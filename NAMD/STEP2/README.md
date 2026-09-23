# STEP2: TDDFT Calculations

This step performs time-dependent density functional theory (TDDFT) calculations along the adiabatic MD trajectories from STEP1. The goal is to generate the Kohn–Sham orbital energies and excitation energies that will be used in STEP3 to build vibronic Hamiltonians and non-adiabatic couplings.

---

## Directory Structure

STEP2/
├── README.md
├── adamantane/
│   ├── AMS.inp
│   ├── UV_vis_spectra_part_1_and_radiative_time.ipynb
│   ├── UV_vis_spectra_part_2.ipynb
│   ├── all_logfiles/
│   ├── all_pdosfiles/
│   ├── distribute_jobs.py
│   ├── run.py
│   └── submit_0.slm
├── 1_methyladamantane/
│   ├── AMS.inp
│   ├── distribute_jobs.py
│   ├── run.py
│   └── submit_1.slm
├── 1_3_dimethyladamantane/
│   ├── AMS.inp
│   ├── distribute_jobs.py
│   ├── run.py
│   └── submit_2.slm
└── 1_3_5_trimethyladamantane/
    ├── AMS.inp
    ├── distribute_jobs.py
    ├── run.py
    └── submit_3.slm

---

## Workflow

### 1. Input File (AMS.inp)

Defines the CP2K / AMS calculation settings, including:
- Level of theory (e.g., B3LYP hybrid functional).
- Basis set (e.g., 6-311++G**).
- Number of excited states (30 electronic states).
- MD trajectory parameters.

### 2. Job Distribution (distribute_jobs.py)

Splits the MD trajectory into independent chunks, each handled by a separate compute job. This enables efficient parallelization on HPC clusters.

### 3. Main Driver (run.py)

Runs the TDDFT calculation for a single trajectory chunk and produces:
- Kohn–Sham orbital energies (E_ks_*.npz).
- Excitation energies and oscillator strengths.
- Log files and PDOS files for post-processing.

### 4. SLURM Submission (submit_*.slm)

SLURM batch scripts used to launch the calculations on the HPC cluster.

---

## Key Outputs

| Output | Description |
|--------|-------------|
| all_logfiles/step_*.log | CP2K / AMS log files for each trajectory frame |
| all_pdosfiles/*.pdos | Projected density of states for each frame |
| E_ks_*.npz | Kohn–Sham orbital energies (used in STEP3) |
| UV-Vis spectra notebooks | Analysis of thermally-averaged UV-Vis spectra and radiative lifetimes |

---

## UV-Vis Spectra Notebooks

- UV_vis_spectra_part_1_and_radiative_time.ipynb: Computes the thermally-averaged UV-Vis spectra and estimates radiative lifetimes using Einstein's formula.
- UV_vis_spectra_part_2.ipynb: Additional analysis of spectral features and comparison with experiment.

---

## How to Run

For each structure:

cd STEP2/<structure>
python distribute_jobs.py    # Generate input chunks
sbatch submit_<n>.slm        # Launch jobs on the cluster

After completion, open the UV-Vis notebooks to analyze the results.

---

## Contact

- Qingxin Zhang (qingxin@buffalo.edu)
- Hamid Zabihi Hashjin (h.zabihi@eng.uk.ac.ir)
- Alexey V. Akimov (alexeyak@buffalo.edu)