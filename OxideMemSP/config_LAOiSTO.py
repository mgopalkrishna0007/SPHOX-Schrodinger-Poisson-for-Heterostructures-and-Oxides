# config.py

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
    "points_per_layer": 10,     
    
    "alpha": 0.05,             
    "alpha_min": 1e-4,
    "alpha_max": 0.3,
    "grow_factor": 1.05,
    "shrink_factor": 0.5,

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
                    "name": "SrO",   # <--- Put SrO First (A-site)
                    "Ec": 0.0, 
                    "Ev": -3.2, 
                    "meff_e": 0.4, 
                    "meff_h": 1.2, 
                    "eps_r": 330.0, 
                    "is_eps_dynamic": False, 
                    "t_layer": t_STO_sublayer,
                    "layer_charge": 0.0   
                },
                {
                    "name": "TiO2",  # <--- Put TiO2 Second (B-site). This will now touch LaO.
                    "Ec": 0.0, 
                    "Ev": -3.2, 
                    "meff_e": 0.4, 
                    "meff_h": 1.2, 
                    "eps_r": 330.0, 
                    "is_eps_dynamic": False, 
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
                    "layer_charge": 1.0   # Alternating +1e charge
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
                    "layer_charge": -1.0  # Alternating -1e charge
                }
            ]
        }
    },

    # --- Stack Definition ---
    "stack_definition": ["STO*50", "LAO*10"],

    "thickness_sweep": {
        "enabled": True,
        "fixed_stack": "STO*50",   # The substrate/fixed layer
        "sweep_mat": "LAO",         # The material layer to vary
        "uc_range": [1, 10]         # [min_unit_cells, max_unit_cells]
    },

    # Replaced by `layer_charge` above, leaving empty
    "dipole_configs": [],
    "doping_profiles": []

}