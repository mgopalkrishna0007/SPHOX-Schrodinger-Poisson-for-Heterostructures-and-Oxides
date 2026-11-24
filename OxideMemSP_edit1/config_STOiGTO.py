# config_STOiGTO.py
# Benchmark of Peelaers et al., Appl. Phys. Lett. 107, 183505 (2015): STO/GTO 2DEG
import numpy as np

# =============================================================================
# Benchmark switches (edit these to walk through the paper's numbers)
# =============================================================================
FIELD_DEPENDENT_EPS = True   # True: Eq. (2) for STO | False: constant eps = 300 (paper's "traditional" case)

# 2DEG density set through the donor sheet in the first GTO monolayer.
# e/uc = n2D[cm^-2] * a_STO^2 ;  a_STO^2 = 1.5249e-15 cm^2
#   3.3e14 cm^-2 -> ~0.50  e/uc   (paper's main case; 0.5 e/uc = 3.28e14)
#   1.0e14 cm^-2 -> ~0.152 e/uc   (Table I, row 2)
#   5.0e13 cm^-2 -> ~0.076 e/uc   (Table I, row 1)
DONOR_E_PER_UC = 0.5


# =============================================================================
# Dielectric Physics Functions
# =============================================================================
def sto_field_dependent_epsilon(E_V_per_m):
    """
    Field-dependent dielectric constant of STO, Eq. (2) of Peelaers et al.
    (integrated D/E form, required because Poisson is solved as div(eps(E) grad V)).
    E_V_per_m: Local Electric field in V/m.
    """
    # Convert E from V/m to MV/cm
    E_MV_cm = np.abs(E_V_per_m) * 1e-8

    # Fitting parameters from the paper (b, f in cm/MV)
    a = 358.50
    b = 1.84
    c = 83.84
    f = 0.15
    # C = 719.307 as printed in the paper makes eps jump 300 -> 237 at 0.27 MV/cm (the SC loop then
    # oscillates). 736.30 makes eps(E) continuous there (same value as the LAOiSTO / STOiSAP configs).
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

# GTO pseudocubic lattice parameter taken equal to STO (coherent growth, ~0.2 nm per monolayer,
# consistent with the paper's 0.2 nm dopant region)
a_GTO = a_STO
t_GTO_sublayer = a_GTO / 2.0

# Unit cell area in-plane
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
    #   "damped": eps from the previous output field, damped (anderson: mix_beta_eps after eps_start;
    #             alpha: with alpha) -- original behaviour
    #   "input" : eps recomputed from the input potential's field every loop after eps_start
    "eps_update": "damped",

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
        "polar":  {"shape": "delta", "sigma_frac": 0.25},
        "doping": {"shape": "delta", "sigma_frac": 0.25},
        "dipole": {"shape": "delta", "sigma_frac": 0.25},
    },

    # --- Boundary Conditions ---
    "bc_left_type": "dirichlet",
    "V_left_applied": 0.0,
    "bc_right_type": "neumann",
    "V_right_applied": 0.0,

    # --- Material Definitions ---
    # Band numbers from the paper: Eg(STO) = 3.27 eV, Eg(GTO) = 2.05 eV, CB offset = 1.47 eV
    # (GTO conduction band above STO's). Reference: Ec(STO) = 0.
    "MATERIALS": {
        "STO": {
            "uc_area": area_uc,
            "layers": [
                {
                    "name": "SrO",
                    "Ec": 0.0,
                    "Ev": -3.27,
                    "meff_e": 1.0,   # paper: 1.0 me (fit to hybrid-DFT DOS)
                    "meff_h": 1.2,   # not given in paper; holes play no role here
                    "eps_r": 300.0,
                    "is_eps_dynamic": FIELD_DEPENDENT_EPS,
                    "eps_dynamic_func": sto_field_dependent_epsilon,
                    "t_layer": t_STO_sublayer,
                    "layer_charge": 0.0
                },
                {
                    "name": "TiO2",
                    "Ec": 0.0,
                    "Ev": -3.27,
                    "meff_e": 1.0,
                    "meff_h": 1.2,
                    "eps_r": 300.0,
                    "is_eps_dynamic": FIELD_DEPENDENT_EPS,
                    "eps_dynamic_func": sto_field_dependent_epsilon,
                    "t_layer": t_STO_sublayer,
                    "layer_charge": 0.0
                }
            ]
        },
        "GTO": {
            "uc_area": area_uc,
            "layers": [
                {
                    "name": "GdO",
                    "Ec": 1.47,           # STO Ec + CB offset
                    "Ev": 1.47 - 2.05,    # Ec - Eg(GTO) = -0.58 eV
                    "meff_e": 1.0,        # not given for GTO; electrons barely penetrate
                    "meff_h": 1.2,
                    "eps_r": 30.0,        # paper: static eps of GTO
                    "is_eps_dynamic": False,
                    "t_layer": t_GTO_sublayer,
                    "layer_charge": 0.0   # 2DEG is introduced by the donor sheet below, not by polar planes
                },
                {
                    "name": "TiO2_GTO",
                    "Ec": 1.47,
                    "Ev": 1.47 - 2.05,
                    "meff_e": 1.0,
                    "meff_h": 1.2,
                    "eps_r": 30.0,
                    "is_eps_dynamic": False,
                    "t_layer": t_GTO_sublayer,
                    "layer_charge": 0.0
                }
            ]
        }
    },

    # --- Stack Definition ---
    # Paper: 500 nm STO / 500 nm GTO. 500 nm / 0.3905 nm = 1280 uc.
    # GTO thickness is irrelevant physically (net-neutral stack -> field-free GTO), so it is
    # shortened to keep the mesh small. Use "GTO*1280" for an exact geometric match.
    "stack_definition": ["STO*50", "GTO*50"],

    "thickness_sweep": {
        "enabled": False,
        "fixed_stack": "STO*1280",
        "sweep_mat": "GTO",
        "uc_range": [1, 10]
    },

    "dipole_configs": [],

    # --- Doping: donor sheet in the first GTO monolayer (~0.2 nm region of the paper) ---
    "doping_profiles": [
        {
            "type": "donor",
            "profile": "delta",            # "delta" | "step" | "gaussian"
            "value": DONOR_E_PER_UC,          # electrons per interface unit cell
            "interface": ("STO", "GTO"),
            "direction": "right",
            "layer_offset": [1, 1]            # offset 0 = last STO layer, 1 = first GTO layer
        }
    ]
}