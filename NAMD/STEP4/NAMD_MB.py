#!/usr/bin/env python3

import os
import time
import multiprocessing as mp

import numpy as np

from liblibra_core import *
import libra_py
from libra_py import units, data_conv
import libra_py.dynamics.tsh.compute as tsh_dynamics
import libra_py.workflows.nbra.decoherence_times as decoherence_times

from recipes import (
    fssh_nbra,
    fssh2_nbra,
    ida_nbra,
    msdm_nbra,
)


# ============================================================
# USER SETTINGS
# ============================================================

ROOT = "/home/hamid/A"

SYSTEMS = [
    ("adamantane",
     os.path.join(ROOT, "step3_0", "adamantane")),

    ("1methyl_adamantane",
     os.path.join(ROOT, "step3_1", "1methyl_adamantane")),

    ("2methyl_adamantane",
     os.path.join(ROOT, "step3_2", "2methyl_adamantane")),

    ("3methyl_adamantane",
     os.path.join(ROOT, "step3_3", "3methyl_adamantane")),
]


# Electronic dynamics
istep = 0
fstep = 3000
DT = 1.0
ISTATE = 1
BASIS_SIZE = 15


# Production statistics
NTRAJ_PER_BATCH = 1000
ICONDS = list(range(1, 3001, 100))   # 30 batches
NTHREADS = 4


# Select ONE method
METHOD = "FSSH"
# METHOD = "FSSH2"
# METHOD = "IDA"
# METHOD = "MSDM"


# ============================================================
# GLOBAL WORKER CONTEXT
# ============================================================

CURRENT_STRUCTURE = None

Ham_adi = None
Hvib_adi = None
NAC_adi = None
Basis_transform = None
St_adi = None

NSTEPS = None
NSTATES = None

MODEL_PARAMS = None
DYN_GENERAL = None
ELEC_PARAMS = None
NUCL_PARAMS = None

PRF = None


# ============================================================
# VALIDATION
# ============================================================

VALID_METHODS = {
    "FSSH",
    "FSSH2",
    "IDA",
    "MSDM",
}

if METHOD not in VALID_METHODS:
    raise ValueError(
        f"Unknown METHOD = {METHOD}. "
        f"Choose one of: {sorted(VALID_METHODS)}"
    )

if len(ICONDS) != 30:
    raise RuntimeError(
        f"Expected 30 independent batches, got {len(ICONDS)}"
    )


# ============================================================
# LIBRA FILE READER
# ============================================================

def read_files_new(
    data_path,
    istep,
    fstep,
    dt,
    _nstates=None
):

    nsteps = fstep - istep

    print(f"Number of steps = {nsteps}")

    # --------------------------------------------------------
    # Determine number of electronic states
    # --------------------------------------------------------

    st0 = np.loadtxt(
        os.path.join(
            data_path,
            f"st_adi_{istep}.txt"
        )
    )

    max_nstates = st0.shape[0]

    if _nstates is None:
        NSTATES_local = max_nstates
    else:
        NSTATES_local = min(
            _nstates,
            max_nstates
        )

    rng = Py2Cpp_int(
        list(range(NSTATES_local))
    )

    St_adi_local = []
    Ham_adi_local = []
    Hvib_adi_local = []
    NAC_adi_local = []
    Basis_transform_local = []

    # --------------------------------------------------------
    # Read all electronic-structure data
    # --------------------------------------------------------

    for i in range(istep, fstep):

        # ----------------------------------------------------
        # Overlap
        # ----------------------------------------------------

        st = CMATRIX(
            max_nstates,
            max_nstates
        )

        st_red = CMATRIX(
            NSTATES_local,
            NSTATES_local
        )

        st.Load_Matrix_From_File(
            os.path.join(
                data_path,
                f"st_adi_{i}.txt"
            )
        )

        pop_submatrix(
            st,
            st_red,
            rng,
            rng
        )

        St_adi_local.append(
            st_red
        )

        # ----------------------------------------------------
        # Hamiltonian
        # ----------------------------------------------------

        ham_adi = CMATRIX(
            max_nstates,
            max_nstates
        )

        ham_adi_red = CMATRIX(
            NSTATES_local,
            NSTATES_local
        )

        ham_adi.Load_Matrix_From_File(
            os.path.join(
                data_path,
                f"ham_adi_{i}.txt"
            )
        )

        pop_submatrix(
            ham_adi,
            ham_adi_red,
            rng,
            rng
        )

        Ham_adi_local.append(
            ham_adi_red
        )

        # ----------------------------------------------------
        # NAC
        # ----------------------------------------------------

        nac_adi_red = (
            st_red - st_red.H()
        ) / (
            2.0 * dt * units.fs2au
        )

        NAC_adi_local.append(
            nac_adi_red
        )

        # ----------------------------------------------------
        # Vibronic Hamiltonian
        # ----------------------------------------------------

        hvib_adi_red = (
            ham_adi_red
            - nac_adi_red * (0.0 + 1.0j)
        )

        Hvib_adi_local.append(
            hvib_adi_red
        )

        # ----------------------------------------------------
        # Basis transformation
        # ----------------------------------------------------

        bas_trans_red = CMATRIX(
            NSTATES_local,
            NSTATES_local
        )

        bas_trans_red.identity()

        Basis_transform_local.append(
            bas_trans_red
        )

    return (
        nsteps,
        NSTATES_local,
        St_adi_local,
        Ham_adi_local,
        Hvib_adi_local,
        NAC_adi_local,
        Basis_transform_local,
    )


# ============================================================
# EMPTY MODEL CLASS
# ============================================================

class abstr_class:
    pass


# ============================================================
# GLOBAL compute_model
# ============================================================

def compute_model(q, params, full_id):

    timestep = params["timestep"]

    obj = abstr_class()

    obj.ham_adi = Ham_adi[timestep]
    obj.hvib_adi = Hvib_adi[timestep]
    obj.nac_adi = NAC_adi[timestep]

    obj.basis_transform = (
        Basis_transform[timestep]
    )

    obj.time_overlap_adi = (
        St_adi[timestep]
    )

    return obj


# ============================================================
# GLOBAL WORKER FUNCTION
# ============================================================

def function1(icond):

    global CURRENT_STRUCTURE
    global MODEL_PARAMS
    global DYN_GENERAL
    global ELEC_PARAMS
    global NUCL_PARAMS
    global PRF

    # Small stagger between workers
    time.sleep(
        (icond % NTHREADS) * 0.2
    )

    rnd = Random()

    mdl = dict(
        MODEL_PARAMS
    )

    mdl.update({
        "icond": icond
    })

    dyn_gen = dict(
        DYN_GENERAL
    )

    prefix = (
        f"{PRF}0_icond_{icond}"
    )

    dyn_gen.update({
        "prefix": prefix,
        "prefix2": prefix,
    })

    print(
        f"[START] {CURRENT_STRUCTURE} | "
        f"{prefix} | "
        f"{NTRAJ_PER_BATCH} trajectories",
        flush=True
    )

    res = tsh_dynamics.generic_recipe(
        dyn_gen,
        compute_model,
        mdl,
        ELEC_PARAMS,
        NUCL_PARAMS,
        rnd
    )

    print(
        f"[DONE ] {CURRENT_STRUCTURE} | "
        f"{prefix}",
        flush=True
    )

    return res


# ============================================================
# RUN ONE STRUCTURE
# ============================================================

def run_structure(
    structure_name,
    data_path
):

    global CURRENT_STRUCTURE

    global Ham_adi
    global Hvib_adi
    global NAC_adi
    global Basis_transform
    global St_adi

    global NSTEPS
    global NSTATES

    global MODEL_PARAMS
    global DYN_GENERAL
    global ELEC_PARAMS
    global NUCL_PARAMS

    global PRF

    CURRENT_STRUCTURE = structure_name

    print("\n")
    print("#" * 70)
    print(
        f"# STARTING STRUCTURE: "
        f"{structure_name}"
    )
    print("#" * 70)

    print(
        f"Data path: {data_path}"
    )

    print(
        f"Method   : {METHOD}"
    )

    print("#" * 70)

    # ========================================================
    # Validate directory
    # ========================================================

    if not os.path.isdir(data_path):

        raise FileNotFoundError(
            f"Data directory not found:\n"
            f"{data_path}"
        )

    # ========================================================
    # Validate required files
    # ========================================================

    required_files = [
        f"st_adi_{istep}.txt",
        f"ham_adi_{istep}.txt",
        f"st_adi_{fstep - 1}.txt",
        f"ham_adi_{fstep - 1}.txt",
    ]

    missing = []

    for fname in required_files:

        path = os.path.join(
            data_path,
            fname
        )

        if not os.path.isfile(path):
            missing.append(fname)

    if missing:

        raise FileNotFoundError(
            f"\nMissing required files for "
            f"{structure_name}:\n"
            + "\n".join(missing)
        )

    # ========================================================
    # Output directory
    # ========================================================

    output_dir = os.path.join(
        data_path,
        METHOD
    )

    os.makedirs(
        output_dir,
        exist_ok=True
    )

    print(
        f"Output directory: {output_dir}"
    )

    # ========================================================
    # Save original working directory
    # ========================================================

    original_cwd = os.getcwd()

    try:

        # ====================================================
        # All relative Libra output goes here
        # ====================================================

        os.chdir(
            output_dir
        )

        print(
            f"Working directory: "
            f"{os.getcwd()}"
        )

        # ====================================================
        # Load data
        # ====================================================

        print(
            "\nLoading electronic-structure data..."
        )

        (
            NSTEPS,
            NSTATES,
            St_adi,
            Ham_adi,
            Hvib_adi,
            NAC_adi,
            Basis_transform,
        ) = read_files_new(
            data_path,
            istep,
            fstep,
            DT,
            BASIS_SIZE
        )

        print(
            f"NSTEPS  = {NSTEPS}"
        )

        print(
            f"NSTATES = {NSTATES}"
        )

        # ====================================================
        # Decoherence times
        # ====================================================

        print(
            "\nCalculating decoherence times..."
        )

        tau, rates = (
            decoherence_times.decoherence_times_ave(
                [Hvib_adi],
                [0],
                len(Hvib_adi),
                0
            )
        )

        # ====================================================
        # Energy gaps
        # ====================================================

        print(
            "Calculating average energy gaps..."
        )

        dE = (
            decoherence_times.energy_gaps_ave(
                [Ham_adi],
                [0],
                len(Ham_adi)
            )
        )

        # ====================================================
        # Save decoherence times
        # ====================================================

        deco_times = (
            data_conv.MATRIX2nparray(tau)
            * units.au2fs
        )

        np.fill_diagonal(
            deco_times,
            0
        )

        np.savetxt(
            "decoherence_times.txt",
            deco_times.real
        )

        # ====================================================
        # Average gaps
        # ====================================================

        gaps = MATRIX(
            NSTATES,
            NSTATES
        )

        for step in range(NSTEPS):

            gaps += dE[step]

        gaps /= NSTEPS

        rates.show_matrix(
            "decoherence_rates.txt"
        )

        gaps.show_matrix(
            "average_gaps.txt"
        )

        # ====================================================
        # Model parameters
        # ====================================================

        MODEL_PARAMS = {
            "timestep": 0,
            "icond": 0,
            "model0": 0,
            "nstates": NSTATES,
        }

        # ====================================================
        # General dynamics parameters
        # ====================================================

        DYN_GENERAL = {

            "nsteps": NSTEPS,

            "ntraj": NTRAJ_PER_BATCH,

            "nstates": NSTATES,

            "dt": (
                DT * units.fs2au
            ),

            "decoherence_rates": rates,

            "ave_gaps": gaps,

            "progress_frequency": 0.1,

            "which_adi_states":
                range(NSTATES),

            "which_dia_states":
                range(NSTATES),

            "mem_output_level": 2,

            "properties_to_save": [
                "timestep",
                "time",
                "se_pop_adi",
                "sh_pop_adi",
            ],

            "prefix": "NBRA",

            "prefix2": "NBRA",

            "isNBRA": 0,

            "nfiles": NSTEPS,

            "maxstep": 1000,
        }

        # ====================================================
        # Select method
        # ====================================================

        if METHOD == "FSSH":

            fssh_nbra.load(
                DYN_GENERAL
            )

            PRF = "FSSH"

        elif METHOD == "FSSH2":

            fssh2_nbra.load(
                DYN_GENERAL
            )

            PRF = "FSSH2"

        elif METHOD == "IDA":

            ida_nbra.load(
                DYN_GENERAL
            )

            PRF = "IDA"

        elif METHOD == "MSDM":

            msdm_nbra.load(
                DYN_GENERAL
            )

            PRF = "MSDM"

        else:

            raise ValueError(
                f"Unknown METHOD: {METHOD}"
            )

        # ====================================================
        # Nuclear parameters
        # ====================================================

        NUCL_PARAMS = {

            "ndof": 1,

            "init_type": 3,

            "q": [-10.0],

            "p": [0.0],

            "mass": [2000.0],

            "force_constant": [0.01],

            "verbosity": -1,
        }

        # ====================================================
        # Electronic parameters
        # ====================================================

        ELEC_PARAMS = {

            "ndia": NSTATES,

            "nadi": NSTATES,

            "verbosity": -1,

            "init_dm_type": 0,

            "init_type": 1,

            "rep": 1,

            "istate": ISTATE,
        }

        # ====================================================
        # Run trajectories
        # ====================================================

        print("\n")
        print("-" * 70)

        print(
            f"Running {len(ICONDS)} batches × "
            f"{NTRAJ_PER_BATCH} trajectories"
        )

        print(
            f"Total trajectories = "
            f"{len(ICONDS) * NTRAJ_PER_BATCH}"
        )

        print(
            f"Processes = {NTHREADS}"
        )

        print("-" * 70)

        start_time = time.time()

        # ----------------------------------------------------
        # IMPORTANT:
        # function1 is GLOBAL, so multiprocessing can pickle it.
        # ----------------------------------------------------

        with mp.Pool(
            processes=NTHREADS
        ) as pool:

            pool.map(
                function1,
                ICONDS
            )

        elapsed = (
            time.time()
            - start_time
        )

        # ====================================================
        # Completion marker
        # ====================================================

        with open(
            "RUN_COMPLETED.txt",
            "w"
        ) as f:

            f.write(
                f"Structure: "
                f"{structure_name}\n"
            )

            f.write(
                f"Method: "
                f"{METHOD}\n"
            )

            f.write(
                f"Active space: "
                f"{BASIS_SIZE}\n"
            )

            f.write(
                f"Trajectories per batch: "
                f"{NTRAJ_PER_BATCH}\n"
            )

            f.write(
                f"Number of batches: "
                f"{len(ICONDS)}\n"
            )

            f.write(
                f"Total trajectories: "
                f"{len(ICONDS) * NTRAJ_PER_BATCH}\n"
            )

            f.write(
                f"Processes: "
                f"{NTHREADS}\n"
            )

            f.write(
                f"Time steps: "
                f"{NSTEPS}\n"
            )

            f.write(
                f"Elapsed seconds: "
                f"{elapsed:.2f}\n"
            )

        print("\n")
        print("#" * 70)

        print(
            f"# COMPLETED: "
            f"{structure_name}"
        )

        print(
            f"# Elapsed: "
            f"{elapsed / 3600:.2f} hours"
        )

        print("#" * 70)

    finally:

        os.chdir(
            original_cwd
        )


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":

    total_start = time.time()

    print("\n")
    print("=" * 70)
    print("NAMD / NBRA PRODUCTION RUN")
    print("=" * 70)

    print(
        f"METHOD              : {METHOD}"
    )

    print(
        f"Active space        : {BASIS_SIZE}"
    )

    print(
        f"Trajectories/batch  : "
        f"{NTRAJ_PER_BATCH}"
    )

    print(
        f"Number of batches   : "
        f"{len(ICONDS)}"
    )

    print(
        f"Total trajectories  : "
        f"{len(ICONDS) * NTRAJ_PER_BATCH}"
    )

    print(
        f"Processes           : "
        f"{NTHREADS}"
    )

    print(
        f"Time steps          : "
        f"{fstep - istep}"
    )

    print("=" * 70)

    print("\nSystems:")

    for name, path in SYSTEMS:

        print(
            f"  {name:22s} -> {path}"
        )

    print("=" * 70)

    # ========================================================
    # Run structures sequentially
    # ========================================================

    for structure_name, data_path in SYSTEMS:

        run_structure(
            structure_name,
            data_path
        )

    # ========================================================
    # Finished
    # ========================================================

    total_elapsed = (
        time.time()
        - total_start
    )

    print("\n")
    print("=" * 70)
    print("ALL FOUR STRUCTURES COMPLETED")
    print("=" * 70)

    print(
        f"Total elapsed time: "
        f"{total_elapsed / 3600:.2f} hours"
    )

    print("=" * 70)