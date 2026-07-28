# config.py
"""
Configuration for Sapphire / STO Heterostructure Simulation.

Interface Profile:
- Sapphire Substrate (Left): Dirichlet Boundary Condition (V = 0.0 V)
- STO Surface (Right): Neumann Boundary Condition
"""

# =============================================================================
# Lattice Parameters and Areas (Based on provided values)
# =============================================================================
# STO (Mat1) parameters
a_STO = 0.39e-9          
t_STO_half = a_STO / 2.0 
area_STO = a_STO**2

# Sapphire / Al2O3 (Mat2) parameters
area_Sapphire = 1.96e-10**2
t_Sapphire = 1.299e-9


# =============================================================================
# Master Simulation Configuration
# =============================================================================
config = {
    # --- Solver Settings ---
    "Nloops": 2000,
    "tol": 1e-10,
    "n": 100,
    "T": 300,
    "fully_ionized_flag": True,
    "E_bind_D": 0.005,
    "E_bind_A": 0.025,
    "points_per_layer": 5,
    
    "alpha": 0.05,
    "alpha_min": 1e-4,
    "alpha_max": 0.3,
    "grow_factor": 1.05,
    "shrink_factor": 0.5,

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
                    "name": "TiO2",  # Previously L1
                    "Ec": -3.86, "Ev": -7.1, 
                    "meff_e": 0.4, "meff_h": 1.2, 
                    "eps_r": 330.0, "is_eps_dynamic": True, 
                    "eps_A": 3.030e-3, "eps_B": 0.0, 
                    "t_layer": t_STO_half,
                    "layer_charge": 0.0  # Assign sheet charge (e per uc) here if needed
                },
                {
                    "name": "SrO",   # Previously L2
                    "Ec": -3.86, "Ev": -7.1, 
                    "meff_e": 0.4, "meff_h": 1.2, 
                    "eps_r": 330.0, "is_eps_dynamic": True, 
                    "eps_A": 3.030e-3, "eps_B": 0.0, 
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
        }
    },

    # --- Stack Definition ---
    # 10 layers of Sapphire substrate + 10 layers of STO film
    "stack_definition": ["Sapphire*20", "STO*77"],

    # --- Interface-Relative Doping Profiles ---
    "doping_profiles": [
        {
            "type": "donor", 
            "profile": "square",
            "value": 1.0 / (area_STO * t_STO_half),
            "interface": ["Sapphire", "STO"],
            "direction": "right",          
            "target_layer": "SrO",          # Targeting the second layer of STO
            "uc_range": [1, 1]             
        }
    ],

    # --- Interface-Relative Dipole Splitting Arrays ---
    "dipole_configs": [
        {
            "interface": ["Sapphire", "STO"],
            "direction": "right",
            "target_layer": "TiO2",         # Targeting the first layer of STO
            "uc_range": [2, 2],            
            "split_m_charge": 4,
            "delta_plus_0": 0.02 * 1e-9,
            "delta_minus_0": 0.0,
            "decay_mode": "none", 
            "decay_rate": 1.0,
            "fixed_delta": "minus"
        }
    ]
}