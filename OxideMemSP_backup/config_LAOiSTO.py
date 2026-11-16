# # config.py
# import numpy as np

# # =============================================================================
# # Dielectric Physics Functions
# # =============================================================================
# def sto_field_dependent_epsilon(E_V_per_m):
#     """
#     Calculates the field-dependent dielectric constant of STO.
#     Equation derived from van der Berg et al.
#     E_V_per_m: Local Electric field in V/m.
#     """
#     # Convert E from V/m to MV/cm
#     E_MV_cm = np.abs(E_V_per_m) * 1e-8
    
#     # Fitting parameters from the paper
#     a = 358.50
#     b = 1.84
#     c = 83.84
#     f = 0.15
#     # C = 719.307
#     C = 736.30


#     if E_MV_cm <= 0.27:
#         return 300.0
#     else:
#         return -(1.0 / E_MV_cm) * ((a / b) * np.exp(-b * E_MV_cm) + (c / f) * np.exp(-f * E_MV_cm) - C)


# # =============================================================================
# # Real Lattice Parameters in (001) Direction
# # =============================================================================
# # STO cubic lattice parameter a = 3.905 Å
# a_STO = 0.3905e-9  
# t_STO_sublayer = a_STO / 2.0  # Monolayer thickness (~0.19525 nm) for TiO2 or SrO

# # LAO pseudocubic out-of-plane lattice parameter c = 3.789 Å
# c_LAO = 0.3789e-9  
# t_LAO_sublayer = c_LAO / 2.0  # Monolayer thickness (~0.18945 nm) for LaO or AlO2

# # Unit cell area in-plane (LAO is epitaxially strained to the STO substrate)
# area_uc = a_STO**2  # ~1.5249e-19 m^2

# #  Charge Dens @ Interface Node : -5.4011e+25 m^-3 (z = 19.51 nm) for 300K
# #  Charge Dens @ Interface Node : -5.3835e+25 m^-3 (z = 19.51 nm) for 1K
# # =============================================================================
# # Master Simulation Configuration
# # =============================================================================
# config = {
#     # --- Solver Settings ---
#     "Nloops": 500,
#     "tol": 1e-10,
#     "n": 100,                  
#     "T": 1,                  
#     "fully_ionized_flag": True,
#     "E_bind_D": 0.005,
#     "E_bind_A": 0.025,
#     "points_per_layer": 10,     
    
#     "alpha": 0.05,             
#     "alpha_min": 1e-4,
#     "alpha_max": 0.3,
#     "grow_factor": 1.05,
#     "shrink_factor": 0.5,

#     # --- Boundary Conditions ---
#     "bc_left_type": "dirichlet", 
#     "V_left_applied": 0.0,
#     "bc_right_type": "neumann",
#     "V_right_applied": 0.0,

#     # --- Material Definitions ---
#     "MATERIALS": {
#         "STO": {
#             "uc_area": area_uc,
#             "layers": [
#                 {
#                     "name": "SrO",
#                     "Ec": 0.0, 
#                     "Ev": -3.2, 
#                     "meff_e": 0.4, 
#                     "meff_h": 1.2, 
#                     "eps_r": 300.0, 
#                     "is_eps_dynamic": True, 
#                     "eps_dynamic_func": sto_field_dependent_epsilon, # <--- Function injected here
#                     "t_layer": t_STO_sublayer,
#                     "layer_charge": 0.0   
#                 },
#                 {
#                     "name": "TiO2", 
#                     "Ec": 0.0, 
#                     "Ev": -3.2, 
#                     "meff_e": 0.4, 
#                     "meff_h": 1.2, 
#                     "eps_r": 300.0, 
#                     "is_eps_dynamic": True, 
#                     "eps_dynamic_func": sto_field_dependent_epsilon, # <--- Function injected here
#                     "t_layer": t_STO_sublayer,
#                     "layer_charge": 0.0  
#                 }
#             ]
#         },
#         "LAO": {
#             "uc_area": area_uc,
#             "layers": [
#                 {
#                     "name": "LaO", 
#                     "Ec": 2.25, 
#                     "Ev": -3.35, 
#                     "meff_e": 0.4, 
#                     "meff_h": 1.2, 
#                     "eps_r": 24.0, 
#                     "is_eps_dynamic": False, 
#                     "t_layer": t_LAO_sublayer,
#                     "layer_charge": 1.0   
#                 },
#                 {
#                     "name": "AlO2", 
#                     "Ec": 2.25, 
#                     "Ev": -3.35, 
#                     "meff_e": 0.4, 
#                     "meff_h": 1.2, 
#                     "eps_r": 24.0, 
#                     "is_eps_dynamic": False, 
#                     "t_layer": t_LAO_sublayer,
#                     "layer_charge": -1.0  
#                 }
#             ]
#         }
#     },

#     # --- Stack Definition ---
#     "stack_definition": ["STO*100", "LAO*100"],

#     "thickness_sweep": {
#         "enabled": False,
#         "fixed_stack": "STO*50",   
#         "sweep_mat": "LAO",         
#         "uc_range": [1, 10]         
#     },

#     "dipole_configs": [],
#     "doping_profiles": []
# }

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
# --- Solver Settings ---
    "Nloops": 500,
    "tol": 1e-7,
    "eps_tol": 1e-3,
    "n": 100,                  
    "T": 1,                  
    "fully_ionized_flag": True,
    "E_bind_D": 0.005,
    "E_bind_A": 0.025,
    "points_per_layer": 10,     
    
    "mix_m": 6, 
    "mix_beta": 0.1, 
    "mix_max_step": 0.1,   
    "mix_beta_eps": 0.1, 
    "eps_start": 20,

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
    "stack_definition": ["STO*100", "LAO*100"],

    "thickness_sweep": {
        "enabled": False,
        "fixed_stack": "STO*50",   
        "sweep_mat": "LAO",         
        "uc_range": [1, 10]         
    },

    "dipole_configs": [],
    "doping_profiles": []
}