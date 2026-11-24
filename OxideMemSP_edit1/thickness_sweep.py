# thickness_sweep.py
"""
Critical-thickness sweep: repeats the full simulation for a range of unit cells of
config["thickness_sweep"]["sweep_mat"] grown on config["thickness_sweep"]["fixed_stack"].
"""
import copy
import numpy as np

from simulation import run_heterostructure_simulation
from plotting import plot_thickness_sweep


def run_thickness_sweep(base_config):
    sweep_cfg = base_config.get("thickness_sweep", {})
    min_uc, max_uc = sweep_cfg["uc_range"]
    fixed_stack = sweep_cfg["fixed_stack"]
    sweep_mat = sweep_cfg["sweep_mat"]

    # STO Lattice constant in cm (approx 3.905 Angstroms)
    a_sto_cm = 3.905e-8
    # In-plane unit cell area at the interface in cm^2
    uc_area_cm2 = a_sto_cm ** 2

    uc_array = np.arange(min_uc, max_uc + 1)
    e_per_uc_array = []

    for uc in uc_array:
        print(f"\n==================================================")
        print(f" SWEEP RUN: {sweep_mat} * {uc} Unit Cells")
        print(f"==================================================")

        # Localize config for this iteration
        cfg_run = copy.deepcopy(base_config)
        cfg_run["stack_definition"] = [fixed_stack, f"{sweep_mat}*{uc}"]

        res = run_heterostructure_simulation(cfg_run)

        # Integrate N_3D (m^-3) over z (m) to get 2D sheet density (m^-2)
        n2d_m2 = np.trapz(res["N_3D"], res["z"])

        # Convert m^-2 to cm^-2
        n2d_cm2 = n2d_m2 * 1e-4

        # Convert cm^-2 to electrons per unit cell (e-/u.c.)
        e_per_uc = n2d_cm2 * uc_area_cm2
        e_per_uc_array.append(e_per_uc)

    # Plotting the Critical Thickness Curve
    plot_thickness_sweep(uc_array, e_per_uc_array, fixed_stack, sweep_mat)

    return uc_array, np.array(e_per_uc_array)
