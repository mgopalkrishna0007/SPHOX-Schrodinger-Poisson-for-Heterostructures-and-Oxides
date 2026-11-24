# master_run.py
"""
Entry point of the Schrödinger-Poisson heterostructure solver.

Step by step:
    1. pick a config (and optionally override the mixing method)
    2. run_heterostructure_simulation   -> stack, mesh, properties, charges, SC loop   (simulation.py)
    3. run_thickness_sweep (if enabled) -> critical-thickness curve                   (thickness_sweep.py)
    4. plot_all_results                 -> plots 1-7                                  (plotting.py)

Module map:
    master_helper.py       find_Ef_general, Schroed1D_dolfin, Poisson_dolfin  (unchanged)
    stack_geometry.py      stack expansion, layer boundaries, interface-relative layer lookup
    mesh_generation.py     exact + uniform z mesh, nodal widths
    material_properties.py band edges, masses, (field-dependent) permittivity on the mesh
    fixed_charges.py       polar planes, doping profiles, split-charge dipoles
    mixers.py              Anderson / adaptive-alpha mixing schemes
    scf_solver.py          FEniCS init + self-consistent loop
    simulation.py          run_heterostructure_simulation (ties the above together)
    thickness_sweep.py     critical-thickness sweep
    plotting.py            all figures
"""
import warnings
import dolfin as df

# from config_BBOiLNB import config
# from config_SIrOiLNB import config
# from config_LAOiSTO import config
# from config_STO100iSTO111 import config
# from config_LLOiSBO import config
from config_STOiSAP import config
# from config_STOiGTO import config

from simulation import run_heterostructure_simulation
from thickness_sweep import run_thickness_sweep
from plotting import plot_all_results

warnings.filterwarnings("ignore")
df.set_log_level(df.LogLevel.WARNING)

# Optional: override the mixing method chosen in the config file ("anderson" or "alpha")
MIXING_METHOD_OVERRIDE = None


if __name__ == "__main__":
    # --- 1. Config / mixing method ---
    if MIXING_METHOD_OVERRIDE is not None:
        config["mixing_method"] = MIXING_METHOD_OVERRIDE

    # --- 2. Single self-consistent run ---
    res = run_heterostructure_simulation(config)

    # --- 3. Thickness sweep (only if enabled in the config) ---
    if config.get("thickness_sweep", {}).get("enabled", False):
        run_thickness_sweep(config)

    # --- 4. Plots ---
    plot_all_results(res, config)
