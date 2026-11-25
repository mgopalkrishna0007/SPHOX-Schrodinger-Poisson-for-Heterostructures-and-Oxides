# config.py
import math

# =============================================================================
# Real Lattice Parameters in (100) and (111) Directions
# =============================================================================
# STO cubic lattice parameter a = 3.905 Å
a_STO = 0.3905e-9  

# --- STO (100) Parameters ---
t_STO100_sublayer = a_STO / 2.0  # Distance between alternating SrO and TiO2 planes
area_uc_100 = a_STO**2           # In-plane unit cell area for (100) facet

# --- STO (111) Parameters ---
# In the [111] direction, planes are alternating Ti and SrO3.
# The distance between identical planes is a/sqrt(3).
t_STO111_sublayer = a_STO / (2.0 * math.sqrt(3)) 
# To conserve unit cell volume (a^3 = area * a/sqrt(3)), the in-plane area is sqrt(3)*a^2
area_uc_111 = math.sqrt(3) * (a_STO**2)


# =============================================================================
# Master Simulation Configuration
# =============================================================================
config = {
    # --- Solver Settings ---
    "Nloops": 500,
    "tol": 1e-10,
    "n": 100,                  
    "T": 300,                  
    "fully_ionized_flag": True,
    "E_bind_D": 0.005,
    "E_bind_A": 0.025,
    "points_per_layer": 5,     
    "print_device_check": True,   # mesh / layer / charge (e/uc) consistency report before the SC loop
    "print_polarization": True,   # inserted-dipole lines + polarization (screening) balance every loop
    
    # --- Mixing Settings ---
    # Choose the self-consistent mixing scheme: "anderson" or "alpha"
    "mixing_method": "anderson",

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
    "alpha_max": 0.03,
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
        "polar":  {"shape": "delta", "sigma_frac": 0.25},
        "doping": {"shape": "delta", "sigma_frac": 0.25},
        "dipole": {"shape": "delta", "sigma_frac": 0.25},
    },

    # --- Boundary Conditions ---
    "bc_left_type": "dirichlet", # Pinned to 0V at the start of STO (100)
    "V_left_applied": 0.0,
    "bc_right_type": "neumann",  # Neumann at the end of STO (111) film
    "V_right_applied": 0.0,

    # --- Material Definitions ---
    "MATERIALS": {
        "STO_100": {
            "uc_area": area_uc_100,
            "layers": [
                {
                    "name": "TiO2", 
                    "Ec": 0.0, 
                    "Ev": -3.2, 
                    "meff_e": 0.4, 
                    "meff_h": 1.2, 
                    "eps_r": 330.0, 
                    "is_eps_dynamic": False, 
                    "t_layer": t_STO100_sublayer,
                    "layer_charge": 0.0   # Neutral (Layer 0 charge)
                },
                {
                    "name": "SrO",  
                    "Ec": 0.0, 
                    "Ev": -3.2, 
                    "meff_e": 0.4, 
                    "meff_h": 1.2, 
                    "eps_r": 330.0, 
                    "is_eps_dynamic": False, 
                    "t_layer": t_STO100_sublayer,
                    "layer_charge": 0.0   # Neutral (Layer 0 charge)
                }
            ]
        },
        "STO_111": {
            "uc_area": area_uc_111,
            "layers": [
                {
                    "name": "Ti", 
                    "Ec": 0.0, 
                    "Ev": -3.2, 
                    "meff_e": 0.4, 
                    "meff_h": 1.2, 
                    "eps_r": 330.0, 
                    "is_eps_dynamic": False, 
                    "t_layer": t_STO111_sublayer,
                    "layer_charge": 4.0   # Ti4+
                },
                {
                    "name": "SrO3", 
                    "Ec": 0.0, 
                    "Ev": -3.2, 
                    "meff_e": 0.4, 
                    "meff_h": 1.2, 
                    "eps_r": 330.0, 
                    "is_eps_dynamic": False, 
                    "t_layer": t_STO111_sublayer,
                    "layer_charge": -4.0  # Sr2+ and 3*(O2-) = -4
                }
            ]
        }
    },


    # --- Stack Definition ---
    "stack_definition": ["STO_100*70", "STO_111*20"],


    # --- Interface-Relative Doping Profiles ---
    "dipole_configs": [],
    "doping_profiles": []
}