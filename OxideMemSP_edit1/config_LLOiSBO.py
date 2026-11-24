# config.py
"""
Configuration for LLO / SBO Heterostructure Simulation.

Interface Profile:
- LLO Substrate (Left): Dirichlet Boundary Condition (V = 0.0 V)
- SBO Film (Right): Neumann Boundary Condition
"""

# =============================================================================
# Lattice Parameters and Areas
# =============================================================================
# LLO pseudocubic lattice constant
a_LLO = 4.18e-10          
t_layer_half = a_LLO / 2.0 
area_uc = a_LLO**2

# =============================================================================
# Master Simulation Configuration
# =============================================================================
config = {
    # --- Solver Settings ---
    "Nloops": 1000,
    "tol": 1e-10,
    "n": 100,
    "T": 300,
    "fully_ionized_flag": True,
    "E_bind_D": 0.005,
    "E_bind_A": 0.025,
    "points_per_layer": 5,
    
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
    "bc_left_type": 'dirichlet', 
    "V_left_applied": 0.0,
    "bc_right_type": 'dirichlet',
    "V_right_applied": 0.0,

    # --- Material Definitions ---
    "MATERIALS": {
        "LLO": {  
            "uc_area": area_uc,
            "layers": [
                {
                    "name": "LaO",
                    "Ec": 2.0, "Ev": -2.5,  # Estimated from Figure 1(b) DOS 
                    "meff_e": 0.4, "meff_h": 1.2, 
                    "eps_r": 24.0, "is_eps_dynamic": False, # Placeholder: not in paper
                    "t_layer": t_layer_half,
                    "layer_charge": 1.0  
                },
                {
                    "name": "LuO2",
                    "Ec": 2.0, "Ev": -2.5, 
                    "meff_e": 0.4, "meff_h": 1.2, 
                    "eps_r": 24.0, "is_eps_dynamic": False, # Placeholder: not in paper
                    "t_layer": t_layer_half,
                    "layer_charge": -1.0  
                }
            ]
        },
        "SBO": {
            "uc_area": area_uc,
            "layers": [
                {
                    "name": "SrO",
                    "Ec": 0.65, "Ev": -0.65, # 1.3 eV gap centered at 0
                    "meff_e": 0.4, "meff_h": 1.2, 
                    "eps_r": 30.0, "is_eps_dynamic": False, # Placeholder: not in paper
                    "t_layer": t_layer_half,
                    "layer_charge": 0.0  
                },
                {
                    "name": "BiO2",
                    "Ec": 0.65, "Ev": -0.65, 
                    "meff_e": 0.4, "meff_h": 1.2, 
                    "eps_r": 30.0, "is_eps_dynamic": False, # Placeholder: not in paper
                    "t_layer": t_layer_half,
                    "layer_charge": 0.0
                }
            ]
        }
    },

    # --- Stack Definition ---
    # "stack_definition": ["SBO*50", "LLO*4"],
    "stack_definition": ["SBO*20","LLO*3", "SBO*20"],
    "doping_profiles": [],
    "dipole_configs": []
}