# config_W_YHO_LSMO_STO.py
"""
SPHOx configuration: W / Y:HfO2(111) / La0.67Sr0.33MnO3(001) [/ SrTiO3(001) substrate], T = 300 K.

ELECTRODE MODEL (this version)
  LSMO and W are metallic electrodes, so their bands are not solved. Each electrode is a fixed
  potential (Dirichlet boundary condition) placed behind a thin, charge-free "screening gap":

      LSMO electrode | LSMO_scr (lambda_LSMO, eps_LSMO) | Y:HfO2 (solved) | W_scr (lambda_W, eps_W) | W electrode
            ^ Dirichlet U_left                                                         Dirichlet U_right ^

  For a Thomas-Fermi metal with screening length lambda and background permittivity eps_e, the screening
  charge sits ~lambda inside the metal; outside the metal this is electrostatically identical to an
  ideal electrode surface moved back by lambda, with a gap of permittivity eps_e in between
  (Zhuravlev et al., PRL 94, 246802 (2005)). With ideal electrodes (lambda = 0) a uniform polarization
  is screened completely and the bands of Y:HfO2 would be flat for every POL; the finite lambda is
  what makes the barriers depend on the polarization direction.
  The STO substrate lies behind the metallic LSMO and does not affect the Y:HfO2 electrostatics,
  so it is not part of the solved stack.

  Fermi level: fixed by the electrodes (zero bias): Ef = -Phi_LSMO. Electrode potentials (electron
  potential energy, eV) follow from the work functions:  U_left = 0,  U_right = Phi_W - Phi_LSMO - V_bias.
  The resulting barriers are then Phi_B = Phi_M - chi_YHO (Schottky-Mott) plus whatever the
  polarization and its screening add.

FERROELECTRIC MODEL
  Every Y:HfO2 (111) unit cell carries the same fixed dipole, so the film is uniformly polarized
  (P = 36 uC/cm^2 [A], Yun et al., Nat. Mater. 21, 903 (2022)). Each cell: +q over one half, -q over the
  other half ("step" dipole, separation t_YHO/2), with q * (t_YHO/2) = P * area_YHO * t_YHO / e, i.e.
  q = 2 * 0.254 e = 0.508 e per cell. Inside the film the dipoles' charges cancel on average, leaving
  the bound sheets +-0.254 e per cell at the two faces (the macroscopic P.n), plus a small intra-cell ripple.
  POL = +1: P points toward LSMO (+ end of every dipole on the LSMO side); POL = -1: toward W; POL = 0: off.

Every number carries a confidence tier; the evidence is in SPHOx_parameter_dossier.pdf.
  [A] measured on this material, primary source checked     [B] close proxy or measured range
  [C] theory / compilation / second-hand                     [D] modelling estimate on a stated assumption
  [E] not found / placeholder -> replace before quoting results.   A trailing + marks a value computed here.
Units: energies in eV with the vacuum level at 0 eV; lengths in m; areas in m^2; masses in m0.

Suggested sweeps (dossier, section 7): YHO chi (-2.0 / -2.4 eV), YHO eps_r (25, 30, 35), POL (-1, 0, +1),
Phi_LSMO (4.2, 4.5, 4.8 eV), Phi_W (4.55, 5.25 eV); here also lambda_LSMO / lambda_W.
"""
import numpy as np

# =============================================================================
# Lattice / geometry
# =============================================================================
# SrTiO3 substrate (sets the in-plane lattice of the coherent LSMO): a = 3.905 A [A] (Vailionis 2011)
a_STO = 3.905e-10
area_STO = a_STO**2

# Y:HfO2, 5 mol% YO1.5, (111) on LSMO(001): a = 5.20, b = c = 5.07 A [A] (Yun 2022)
t_YHO = 2.952e-10                    # [A+] d111 from 1/d^2 = 1/a^2 + 1/b^2 + 1/c^2
area_YHO = 1.132e-19                 # [A+] area per HfO2 in a (111) plane = (a*b*c/4) / d111
N_YHO = 34                           # 34 cells = 10.0 nm, matching the 9-11 nm films of Yun 2022 [A]

# =============================================================================
# Electrodes
# =============================================================================
# LSMO: work function 4.5 eV, midpoint of 4.2-4.8 eV [B] (de Jong et al. 2003)
Phi_LSMO = 4.5
# Thomas-Fermi screening length 0.30 nm from n = 5.65e21 cm^-3 [D], m* = 2.5 [B], eps_r = 30 [C];
# Hikita PRB 79, 073101 (2009) quotes ~0.31 nm [C+]
lambda_LSMO = 0.30e-9
eps_LSMO = 30.0

# W: work function 4.55 eV polycrystalline [A] (Michaelson 1977); (110) texture 5.05-5.47 eV [B]
Phi_W = 4.55
# Thomas-Fermi length of a free-electron gas with 2 e/atom (n = 1.26e23 cm^-3, E_F = 9.17 eV), eps = 1:
# lambda = sqrt(2 eps0 E_F / (3 e^2 n)) = 0.052 nm [E] (electron count assumed)
lambda_W = 0.052e-9
eps_W = 1.0

# Applied bias on W relative to LSMO (V). Only 0 is an equilibrium state; the Fermi level below is
# the LSMO one, so keep V_bias = 0 unless you know what a non-equilibrium occupation means here.
V_bias = 0.0

# =============================================================================
# Ferroelectric polarization -> one dipole per Y:HfO2 cell
# =============================================================================
# Pr = 36 uC/cm^2 at 300 K [A] (Yun 2022) -> bound charge 0.36 C/m^2 * area_YHO / e = 0.254 e per cell [A+]
POL = +1
sigma_P = 0.254
q_dipole = 2.0 * sigma_P             # "step" dipole: +q / -q over the two half-cells, separation t_YHO / 2

# =============================================================================
# Master simulation configuration
# =============================================================================
config = {
    # --- Solver settings ---
    "Nloops": 2000,
    "tol": 1e-10,
    "n": 100,
    "T": 300,
    "fully_ionized_flag": True,
    "E_bind_D": 0.005,
    "E_bind_A": 0.025,
    "points_per_layer": 5,
    "print_device_check": True,   # mesh / layer / charge (e/uc) consistency report before the SC loop
    "print_polarization": True,   # inserted-dipole lines + polarization (screening) balance every loop

    # Fermi level fixed by the electrodes (eV, vacuum = 0). Remove / set None to use global charge
    # neutrality instead (only meaningful when the electrodes are part of the solved stack).
    "fermi_level": -Phi_LSMO,

    # --- Mixing Settings ---
    "mixing_method": "anderson",
    "mix_m": 6,
    "mix_beta": 0.1,
    "mix_max_step": 0.1,
    "mix_beta_eps": 0.1,
    "eps_start": 20,
    "eps_tol": 1e-3,
    "eps_update": "input",
    "alpha": 0.05,
    "alpha_min": 1e-4,
    "alpha_max": 0.3,
    "grow_factor": 1.05,
    "shrink_factor": 0.5,

    # --- Fixed Charge Profiles ---
    "charge_profile_defaults": {
        "polar":  {"shape": "delta", "sigma_frac": 0.25},
        "doping": {"shape": "delta", "sigma_frac": 0.25},
        "dipole": {"shape": "step", "sigma_frac": 0.25},
    },

    # --- Plots to draw (True = draw, False = skip) ---
    "plots": {
        "band_diagram_wavefunctions": False,          # Band Diagram and Wavefunctions
        "carrier_distribution_heatmap": True,         # Carrier Distribution Heatmap
        "charge_distribution_shaded": True,           # Charge Distribution with Layer Shading
        "conduction_band_wavefunctions": False,       # Conduction Band and Electron Wavefunctions (fixed window)
        "band_diagram_carrier_density": True,         # Band Diagram and Carrier Density
        "charge_distribution": True,                  # Charge Distribution
        "electric_field_dielectric_constant": True,   # Electric Field and Dielectric Constant
    },

    # --- Plot Settings ---
    "plot_settings": {
        "material_colors": {
            "LSMO_scr": "#2a9d8f",
            "YHO": "#9b5de5",
            "W_scr": "#8d99ae",
        },
        "background_alpha": 0.35,
        "linewidth": 2.0,
        "band_ylim_exclude": [],
        "max_real_gap_eV": 20.0,
        "band_ylim": None,
        "density_floor": None,
        "density_decades": 1e-10,
        "charge_linear_ylim": None,
        "charge_xlim": None,
    },

    # --- Boundary conditions: the two electrodes (electron potential energy, eV) ---
    "bc_left_type": 'dirichlet',
    "V_left_applied": 0.0,                              # LSMO electrode (reference)
    "bc_right_type": 'dirichlet',
    "V_right_applied": Phi_W - Phi_LSMO - V_bias,       # W electrode: contact potential difference - bias

    # --- Material definitions ---
    "MATERIALS": {
        # Screening gaps: no carriers (band edges far away -> hard walls for the wavefunctions),
        # no charge, only permittivity. They stand for the screening region inside each metal.
        "LSMO_scr": {
            "uc_area": area_STO,
            "layers": [
                {
                    "name": "LSMO_screen",
                    "Ec": 50.0, "Ev": -50.0,
                    "meff_e": 1.0, "meff_h": 1.0,
                    "eps_r": eps_LSMO, "is_eps_dynamic": False,
                    "t_layer": lambda_LSMO,
                    "layer_charge": 0.0
                }
            ]
        },
        "YHO": {
            "uc_area": area_YHO,
            "layers": [
                {
                    "name": "HfO2",
                    # Ec = -chi, chi = 2.0 +- 0.25 eV [B]: Monaghan et al., Solid-State Electron. 53, 438 (2009)
                    #   (undoped HfO2; CNL-model table value 2.4 eV [C] Robertson & Falabretti 2006) -> sweep 2.0-2.4
                    # Ev = Ec - 6.0 eV gap [C]: literature spread 5.6-6.3 eV
                    "Ec": -2.00, "Ev": -8.00,
                    # meff_e 0.11 +- 0.03 tunnelling mass [B]: Monaghan 2009; meff_h 0.58 [C]: secondary quote only
                    "meff_e": 0.11, "meff_h": 0.58,
                    # eps_r [E]: no usable source (Pavoni 2022 value is the electronic eps1(0), not static) -> sweep 25-35
                    "eps_r": 30.0, "is_eps_dynamic": False,
                    "t_layer": t_YHO,
                    "layer_charge": 0.0           # neutral HfO2 planes; the polarization is in dipole_configs
                }
            ]
        },
        "W_scr": {
            "uc_area": area_YHO,
            "layers": [
                {
                    "name": "W_screen",
                    "Ec": 50.0, "Ev": -50.0,
                    "meff_e": 1.0, "meff_h": 1.0,
                    "eps_r": eps_W, "is_eps_dynamic": False,
                    "t_layer": lambda_W,
                    "layer_charge": 0.0
                }
            ]
        },
    },

    # --- Stack: LSMO electrode (BC) | screening gap | Y:HfO2 | screening gap | W electrode (BC) ---
    "stack_definition": ["LSMO_scr*1", f"YHO*{N_YHO}", "W_scr*1"],

    "doping_profiles": [],

    # --- Ferroelectric: one identical dipole in every Y:HfO2 cell ---
    "dipole_configs": [
        {
            "interface": ["LSMO_scr", "YHO"],
            "direction": "right",
            "layer_offset": [1, N_YHO],          # offset 1 = first YHO cell ... N_YHO = last
            "profile": "step",                   # + over one half-cell, - over the other (separation t_YHO/2)
            "value": q_dipole * POL,             # POL = +1: + half on the LSMO side (P toward LSMO); -1 flips it
            "delta_plus_0": 0.0,                 # not used by the step profile
            "delta_minus_0": 0.0,
            "decay_mode": "none",
            "decay_rate": 1.0,
        }
    ],
}
