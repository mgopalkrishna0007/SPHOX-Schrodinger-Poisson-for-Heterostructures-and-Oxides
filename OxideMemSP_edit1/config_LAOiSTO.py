# config.py
import numpy as np

# =============================================================================
# Dielectric Physics Functions
# =============================================================================
def sto_field_dependent_epsilon(E_V_per_m):
    """
    Calculates the field-dependent dielectric constant of STO.
    Equation derived from van der Berg et al.
    E_V_per_m: Local Electric field in V/m.
    """
    # Convert E from V/m to MV/cm
    E_MV_cm = np.abs(E_V_per_m) * 1e-8
    
    # Fitting parameters from the paper
    a = 358.50
    b = 1.84
    c = 83.84
    f = 0.15
    # C = 719.307
    C = 736.30


    if E_MV_cm <= 0.27:
        return 300.0
    else:
        return -(1.0 / E_MV_cm) * ((a / b) * np.exp(-b * E_MV_cm) + (c / f) * np.exp(-f * E_MV_cm) - C)


# =============================================================================
# Real Lattice Parameters in (001) Direction
# =============================================================================
# STO cubic lattice parameter a = 3.905 Å
a_STO = 0.3905e-9  
t_STO_sublayer = a_STO / 2.0  # Monolayer thickness (~0.19525 nm) for TiO2 or SrO

# LAO pseudocubic out-of-plane lattice parameter c = 3.789 Å
c_LAO = 0.3789e-9  
t_LAO_sublayer = c_LAO / 2.0  # Monolayer thickness (~0.18945 nm) for LaO or AlO2

# Unit cell area in-plane (LAO is epitaxially strained to the STO substrate)
area_uc = a_STO**2  # ~1.5249e-19 m^2

#  Charge Dens @ Interface Node : -5.4011e+25 m^-3 (z = 19.51 nm) for 300K
#  Charge Dens @ Interface Node : -5.3835e+25 m^-3 (z = 19.51 nm) for 1K
# =============================================================================
# Master Simulation Configuration
# =============================================================================
config = {
    # --- Solver Settings ---
    "Nloops": 500,
    "tol": 1e-7,
    "n": 100,                  
    "T": 1,                  
    "fully_ionized_flag": True,
    "E_bind_D": 0.005,
    "E_bind_A": 0.025,
    "points_per_layer": 10,     
    
    # --- Mixing Settings ---
    # Choose the self-consistent mixing scheme: "anderson" or "alpha"
    "mixing_method": "alpha",

    # Anderson mixing (used when mixing_method == "anderson")
    "mix_m": 6,             # number of previous iterates kept in the Anderson history
    "mix_beta": 0.1,        # linear mixing fraction inside the Anderson step
    "mix_max_step": 0.1,    # cap on the largest potential change per loop (V)
    "mix_beta_eps": 0.1,    # fixed damping of the field-dependent dielectric constant
    "eps_start": 20,        # eps is only updated after this many loops (lets V settle first)
    "eps_tol": 1e-3,        # relative eps change required for convergence

    # Field-dependent eps update (only matters for layers with "is_eps_dynamic": True)
    #   "input" : eps recomputed from the input potential's field every loop after eps_start
    #             (recommended: anderson then converges reliably, e.g. STOiGTO 87 loops vs never)
    #   "damped": eps from the previous output field, damped (anderson: mix_beta_eps after eps_start;
    #             alpha: with alpha) -- original behaviour
    "eps_update": "input",

    # Alpha mixing (adaptive linear mixing, used when mixing_method == "alpha")
    "alpha": 0.05,           # initial mixing fraction (V and eps)
    "alpha_min": 1e-4,
    "alpha_max": 0.3,
    "grow_factor": 1.05,     # alpha *= grow_factor when the error decreases
    "shrink_factor": 0.5,    # alpha *= shrink_factor when the error increases

    # --- Fixed Charge Profiles ---
    # Shape of every fixed charge inside its layer: "delta" | "step" | "gaussian"
    #   delta    : on the atomic plane (layer centre); dipoles: on centre -/+ delta_plus/minus
    #   step     : uniform over the whole layer;        dipoles: + over left half, - over right half
    #   gaussian : centred on the layer centre (dipoles: on the +/- points), truncated to the layer,
    #              sigma = sigma_frac * t_layer
    # Per-charge overrides: "profile" (+ optional "sigma_frac" or "sigma" in m) inside a
    # doping_profiles / dipole_configs entry, or "charge_profile" inside a layer dict (polar planes).

    # Dipole layers keep their constant "eps_r" even when the layer is field dependent (eps(E)
    # cannot carry the D of a dipole of more than ~0.6 e/uc). Add "eps_dynamic_inside": True to a
    # dipole_configs entry to let eps(E) act inside that dipole layer anyway.
    "charge_profile_defaults": {
        "polar":  {"shape": "step", "sigma_frac": 0.25},
        "doping": {"shape": "step", "sigma_frac": 0.25},
        "dipole": {"shape": "step", "sigma_frac": 0.25},
    },

    # --- Boundary Conditions ---
    "bc_left_type": "dirichlet", 
    "V_left_applied": 0.0,
    "bc_right_type": "neumann",
    "V_right_applied": 0.0,

    # --- Material Definitions ---
    "MATERIALS": {
        "STO": {
            "uc_area": area_uc,
            "layers": [
                {
                    "name": "SrO",
                    "Ec": 0.0, 
                    "Ev": -3.2, 
                    "meff_e": 0.4, 
                    "meff_h": 1.2, 
                    "eps_r": 300.0, 
                    "is_eps_dynamic": True, 
                    "eps_dynamic_func": sto_field_dependent_epsilon, # <--- Function injected here
                    "t_layer": t_STO_sublayer,
                    "layer_charge": 0.0   
                },
                {
                    "name": "TiO2", 
                    "Ec": 0.0, 
                    "Ev": -3.2, 
                    "meff_e": 0.4, 
                    "meff_h": 1.2, 
                    "eps_r": 300.0, 
                    "is_eps_dynamic": True, 
                    "eps_dynamic_func": sto_field_dependent_epsilon, # <--- Function injected here
                    "t_layer": t_STO_sublayer,
                    "layer_charge": 0.0  
                }
            ]
        },
        "LAO": {
            "uc_area": area_uc,
            "layers": [
                {
                    "name": "LaO", 
                    "Ec": 2.25, 
                    "Ev": -3.35, 
                    "meff_e": 0.4, 
                    "meff_h": 1.2, 
                    "eps_r": 24.0, 
                    "is_eps_dynamic": False, 
                    "t_layer": t_LAO_sublayer,
                    "layer_charge": 1.0   
                },
                {
                    "name": "AlO2", 
                    "Ec": 2.25, 
                    "Ev": -3.35, 
                    "meff_e": 0.4, 
                    "meff_h": 1.2, 
                    "eps_r": 24.0, 
                    "is_eps_dynamic": False, 
                    "t_layer": t_LAO_sublayer,
                    "layer_charge": -1.0  
                }
            ]
        }
    },

    # --- Stack Definition ---
    "stack_definition": ["STO*100", "LAO*20"],

    "thickness_sweep": {
        "enabled": False,
        "fixed_stack": "STO*50",   
        "sweep_mat": "LAO",         
        "uc_range": [1, 10]         
    },

    "dipole_configs": [],
    "doping_profiles": []
}