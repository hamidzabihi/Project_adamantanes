# STEP1: Structure Preparation and Adiabatic Molecular Dynamics

This step prepares the initial molecular structures of adamantane and its methylated derivatives, aligns them, and runs equilibrated adiabatic molecular-dynamics (MD) simulations at 300 K. The resulting trajectories serve as input for the TDDFT calculations in STEP2.

---

## Directory Structure

STEP1/
├── adamantane/
│   └── aligned.tar.gz
├── 1_methyladamantane/
│   └── aligned.tar.gz
├── 1_3_dimethyladamantane/
│   └── aligned.tar.gz
├── 1_3_5_trimethyladamantane/
│   └── aligned.tar.gz
└── align.py

---

## Contents

| File / Folder | Description |
|---------------|-------------|
| adamantane/aligned.tar.gz | Aligned structures for pristine adamantane |
| 1_methyladamantane/aligned.tar.gz | Aligned structures for 1-methyl-adamantane |
| 1_3_dimethyladamantane/aligned.tar.gz | Aligned structures for 1,3-dimethyl-adamantane |
| 1_3_5_trimethyladamantane/aligned.tar.gz | Aligned structures for 1,3,5-trimethyl-adamantane |
| align.py | Python script used to align the molecular structures |

---

## Workflow

### 1. Structure Alignment

Run the alignment script:

python align.py

The script reads the raw structural files, aligns them to a common reference frame, and writes the aligned structures to the corresponding aligned.tar.gz archives.

### 2. Adiabatic Molecular Dynamics (MD)

Using the aligned structures, run adiabatic MD simulations at 300 K to obtain equilibrated ground-state trajectories.

**Objective:**
Obtain equilibrated adiabatic trajectories at 300 K.

**Output:**
Equilibrated ground-state molecular-dynamics trajectories for subsequent electronic-structure calculations.

---

## Restarting a CP2K MD Simulation

To continue a previously interrupted CP2K molecular-dynamics simulation, include the following restart settings in the input file:

&EXT_RESTART
  RESTART_FILE_NAME adamantane_MD-1_5000.restart
  RESTART_DEFAULT T
&END EXT_RESTART

---

## Purpose

Consistent alignment and equilibrated MD trajectories are essential for:
- Reliable comparison of electronic properties across the methylated series.
- Stable molecular-dynamics trajectories in STEP2.
- Consistent active-space selection in STEP3.

---

## Contact

- Qingxin Zhang (qingxin@buffalo.edu)
- Hamid Zabihi Hashjin (h.zabihi@eng.uk.ac.ir)
- Alexey V. Akimov (alexeyak@buffalo.edu)