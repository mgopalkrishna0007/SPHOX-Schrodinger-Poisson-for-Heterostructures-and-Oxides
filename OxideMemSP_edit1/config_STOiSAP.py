# config.py
"""
Configuration for Sapphire / STO Heterostructure Simulation.Hi

Interface Profile:
- Sapphire Substrate (Left): Dirichlet Boundary Condition (V = 0.0 V)
- STO Surface (Right): Neumann Boundary Condition
"""

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
    C = 736.30
    
    if E_MV_cm <= 0.27:
        return 300.0
    else:
        return -(1.0 / E_MV_cm) * ((a / b) * np.exp(-b * E_MV_cm) + (c / f) * np.exp(-f * E_MV_cm) - C)


# =============================================================================
# Lattice Parameters and Areas (Based on provided values)
# =============================================================================
# STO (Mat1) parameters
a_STO = 0.39e-9          
t_STO_half = a_STO / 2.0 
area_STO = a_STO**2

# Sapphire / Al2O3 (Mat2) parameters
area_Sapphire = (np.sqrt(3)/2)*1.96e-10**2
t_Sapphire = 1.299e-9


# =============================================================================
# Master Simulation Configuration
# =============================================================================
config = {
    # --- Solver Settings ---
    "Nloops": 2000,
    "tol": 1e-6,
    "n": 200,
    "T": 300,
    "fully_ionized_flag": True,
    "E_bind_D": 0.005,
    "E_bind_A": 0.025,
    "points_per_layer": 10,
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
    #             (recommended: anderson then converges reliably)
    #   "damped": eps from the previous output field, damped (anderson: mix_beta_eps after eps_start;
    #             alpha: with alpha) -- original behaviour
    "eps_update": "input",
    # "mix_m": 3e-4,
    # "mix_beta": 0.05,
    # "mix_max_step": 0.05,
    # "mix_beta_eps": 0.1,
    # "eps_start": 20,

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

    # --- Plots to draw (True = draw, False = skip) ---
    "plots": {
        "band_diagram_wavefunctions": False,          # Band Diagram and Wavefunctions
        "carrier_distribution_heatmap": True,         # Carrier Distribution Heatmap
        "charge_distribution_shaded": True,           # Charge Distribution with Layer Shading
        "conduction_band_wavefunctions": True,        # Conduction Band and Electron Wavefunctions
        "band_diagram_carrier_density": True,         # Band Diagram and Carrier Density
        "charge_distribution": True,                  # Charge Distribution
        "electric_field_dielectric_constant": True,   # Electric Field and Dielectric Constant
    },

    # --- Plot Settings ---
    "plot_settings": {
        # Shading color of each material in the colored plots (any matplotlib color);
        # materials not listed get one automatic color per atomic layer
        "material_colors": {
            "STO": "#f4a261",
            "Sapphire": "#4a90d9",
            "Vac": "#d9d9d9",
            "SapphireTop": "#2f6fb5",
        },
        "background_alpha": 0.35,      # opacity of the layer shading (0 = invisible, 1 = solid)
        # Band plots: y-axis spans the material with the largest band gap (vacuum ignored).
        "band_ylim_exclude": [],       # materials to ignore for that, e.g. ["Vac"]
        "max_real_gap_eV": 20.0,       # gaps above this are treated as vacuum/barrier and ignored
        "band_ylim": None,             # or (ymin, ymax) in eV to fix the window by hand
        # Density plots (log scale, m^-3): lower limit
        "density_floor": None,         # fixed lower limit in m^-3, or None -> top * density_decades
        "density_decades": 1e-10,      # used when density_floor is None
        "charge_linear_ylim": None,    # charge-distribution left (linear) axis: (ymin, ymax) in m^-3, or None
        "charge_xlim": None,           # charge-distribution plots: (zmin, zmax) in nm to zoom in, or None
    },

    # --- Boundary Conditions ---
    "bc_left_type": 'dirichlet', 
    "V_left_applied": 0.0,
    "bc_right_type": 'neumann',
    "V_right_applied": 0.0,

    # --- Material Definitions ---
    "MATERIALS": {
        "STO": {  # Previously Mat1
            "uc_area": area_STO,
            "layers": [
                {
                    "name": "SrO",   # Previously L2
                    "Ec": -3.86, "Ev": -7.1, 
                    "meff_e": 0.7, "meff_h": 1.0, 
                    "eps_r": 300.0, 
                    "is_eps_dynamic": True, # while checking for bare band alignments, set to False to use a constant eps_r=300.0
                    "eps_dynamic_func": sto_field_dependent_epsilon, 
                    "t_layer": t_STO_half,
                    "layer_charge": 0.0  
                }
                ,           
                {
                    "name": "TiO2",  # Previously L1
                    "Ec": -3.86, "Ev": -7.1, 
                    "meff_e": 0.7, "meff_h": 1.0, 
                    "eps_r": 300.0, 
                    "is_eps_dynamic": True, # while checking for bare band alignments, set to False to use a constant eps_r=300.0
                    "eps_dynamic_func": sto_field_dependent_epsilon, 
                    "t_layer": t_STO_half,
                    "layer_charge": 0.0  
                }

            ]
        },
        "Sapphire": {  # Previously Mat2
            "uc_area": area_Sapphire,
            "layers": [
                {
                    "name": "Al2O3", # Previously L3
                    "Ec": -1.0, "Ev": -9.9, 
                    "meff_e": 0.41, "meff_h": 0.36, 
                    "eps_r": 9.3, "is_eps_dynamic": False, 
                    "t_layer": t_Sapphire,
                    "layer_charge": 0.0  
                }
            ]
        },
        "Vac": {
            "uc_area": area_STO,
            "layers": [
                {
                    "name": "Vacuum", 
                    "Ec": 100.0, "Ev": -100.0, 
                    "meff_e": 1.0, "meff_h": 1.0, 
                    "eps_r": 1.0, "is_eps_dynamic": False, 
                    "t_layer": 2e-9,
                    "layer_charge": 0.0
                }
            ]
        },

        "SapphireTop": {"uc_area": area_Sapphire,   # fix to 1.96e-19
            "layers": [
                {"name":"Al2O3",  "Ec":-1.0,"Ev":-9.9,"meff_e":0.41,"meff_h":0.36,
                "eps_r":9.3,"is_eps_dynamic":False,"t_layer":1.299e-9-0.2165e-9,"layer_charge":0.0},
                {"name":"Al_IF",  "Ec":-2.0,"Ev":-9.9,"meff_e":0.41,"meff_h":0.36,
                "eps_r":9.3,"is_eps_dynamic":False,"t_layer":0.2165e-9,"layer_charge":0.0}
            ]}

    },

    # --- Stack Definition ---
    # 10 layers of Sapphire substrate + 10 layers of STO film
    # "stack_definition": ["Sapphire*30", "STO*77"],

    "stack_definition": ["Sapphire*10","STO*77", "Vac*1"],


    "thickness_sweep": {
        "enabled": False,
        "fixed_stack": "Sapphire*30",   
        "sweep_mat": "STO",         
        "uc_range": [10, 77]         
    },
        # --- Interface-Relative Doping Profiles ---
    "doping_profiles": [
    {
        "type": "donor",
        "profile": "step",            # "delta" | "step" | "gaussian"
        "value": +0.3,                  # e/uc = 0.2 V_O²⁺ per uc; scan 0.25–0.55
        "interface": ["Sapphire", "STO"],
        "direction": "right",
        "layer_offset": [1, 1]         # TiO2#1, the first Ti plane
    },
        
    {
        "type": "donor",
        "profile": "step",            # "delta" | "step" | "gaussian"
        "value": -0.0,                  # e/uc = 0.2 V_O²⁺ per uc; scan 0.25–0.55
        "interface": ["Sapphire", "STO"],
        "direction": "right",
        "layer_offset": [2, 2]         # TiO2#1, the first Ti plane
    },
    
    {
        "type": "donor",
        "profile": "step",            # "delta" | "step" | "gaussian"
        "value": +0.0,                  # e/uc = 0.2 V_O²⁺ per uc; scan 0.25–0.55
        "interface": ["Sapphire", "STO"],
        "direction": "right",
        "layer_offset": [154, 154]         # TiO2#1, the first Ti plane
    }
    ],


    "dipole_configs": [
        
        {
        "interface": ["Sapphire", "STO"],
        "direction": "right",
        "profile": "step",            # "delta" | "step" | "gaussian"
        "layer_offset": [4, 4],        # centred on the Ti plane of UC1
        "value": 0.0,                  # e/uc, see below
        "delta_plus_0": 0.1e-9,
        "delta_minus_0": 0.1e-9,     # total separation 0.19 nm
        "decay_mode": "none",
        "decay_rate": 2
    }

    ]
}
