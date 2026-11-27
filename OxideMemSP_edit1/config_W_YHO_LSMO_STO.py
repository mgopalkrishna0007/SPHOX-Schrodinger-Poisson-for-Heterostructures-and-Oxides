# config_W_YHO_LSMO_STO.py
"""
SPHOx configuration: W / Y:HfO2(111) / La0.67Sr0.33MnO3(001) / SrTiO3(001), T = 300 K.

Built from SPHOx_parameter_dossier.pdf (Appendix A, compiled 27 Nov 2026). The physics values are
the dossier's; only the solver / mixing / plotting blocks were added to match the current code.

Interface profile:
- SrTiO3 substrate (left): Dirichlet boundary condition (V = 0.0 V)
- W top electrode (right): Neumann boundary condition (no field inside the metal)

Every number carries a confidence tier; the evidence is in SPHOx_parameter_dossier.pdf.
  [A] measured on this material, primary source checked
  [B] measured on a close proxy, or a measured range
  [C] theory, compilation, or a value only quoted second-hand
  [D] modelling estimate resting on a stated assumption
  [E] not found, unverified or placeholder -> replace before quoting results
  A trailing + (e.g. [A+]) marks a value computed here from cited inputs; it inherits
  the tier of its weakest input.
Units: energies in eV with the vacuum level at 0 eV; lengths in m; areas in m^2; masses in m0.

Suggested sweeps (dossier, section 7):
  YHO chi (Ec)        -2.0, -2.4 eV          measured vs theory; sets every YHO barrier
  YHO eps_r           25, 30, 35              no usable source
  POL                 -1, 0, +1 (and +-0.5)   up/down asymmetry; intrinsic vs ionic polarization
  LSMO work function  4.2, 4.5, 4.8 eV        surface-treatment dependence
  W work function     4.55, 5.25 eV           polycrystalline vs (110) texture
  STO eps_B           0 first, then fitted    only an effective value exists at 300 K
"""
import numpy as np

# =============================================================================
# Lattice parameters and areas
# =============================================================================
# SrTiO3: a = 3.905 A [A] (Vailionis et al., PRB 83, 064101 (2011), Table I)
a_STO = 3.905e-10
t_STO_half = a_STO / 2.0             # [A+] AO / BO2 plane spacing, 1.953 A
area_STO = a_STO**2                  # [A+] 1.525e-19 m^2

# La0.67Sr0.33MnO3: bulk R-3c a_r = 5.471 A, alpha = 60.43 deg [A] (Radaelli 1997, via Vailionis 2011)
#   volume-equivalent pseudocubic a_pc = 3.881 A [A+]
#   coherent on STO: in-plane a = a_STO (+0.62 % tensile)                     [D] assumes no relaxation
#   out-of-plane c_pc = 3.849 A from elasticity with Poisson ratio 0.4       [D] nu assumed
area_LSMO = a_STO**2
t_LSMO_half = 3.849e-10 / 2.0

# Y:HfO2, 5 mol% YO1.5, (111) on LSMO(001): a = 5.20, b = c = 5.07 A [A] (Yun et al., Nat. Mater. 21, 903 (2022))
t_YHO = 2.952e-10                    # [A+] d111 from 1/d^2 = 1/a^2 + 1/b^2 + 1/c^2
area_YHO = 1.132e-19                 # [A+] area per HfO2 in a (111) plane = (a*b*c/4) / d111

# W: alpha-bcc a = 3.165 A [B]; (110) planes, one atom per in-plane cell
t_W = 2.238e-10                      # [B+] d110 = a / sqrt(2)
area_W = 7.083e-20                   # [B+] a^2 * sqrt(2) / 2

# =============================================================================
# Derived charges (see dossier, section 6)
# =============================================================================
# Ferroelectric bound charge: Pr = 36 uC/cm^2 at 300 K [A] (Yun 2022) -> sigma_b = P.n
#   0.36 C/m^2 * area_YHO / e = 0.254 e per (111) cell [A+]
# POL = +1: P points down, toward LSMO (+ sheet at the LSMO face)
# POL = -1: P points up, toward W;  POL = 0: paraelectric reference run
POL = +1
sigma_P = 0.254 * POL

# LSMO electron gas with a compensating positive background ("jellium"):
#   n = x / a_pc^3 = 0.33 / (3.881 A)^3 = 5.65e21 cm^-3 [D]
#   (Hikita 2009 uses 5.1e21 cm^-3 = 0.30 / a_pc^3 for x = 0.3, no source cited)
#   background per half-cell layer = n * t_LSMO_half * area_LSMO = 0.166 e [D]
jel_LSMO = 0.166
# W: free-electron gas with 2 e per atom = 1.26e23 cm^-3 -> 2.0 e per (110) layer [E] assumption
jel_W = 2.0

# =============================================================================
# Dielectric functions
# =============================================================================
# STO: eps_r(E) = 1/(A + B|E|) (form used by Copie et al. 2009)
#   A = 1/330 = 3.030e-3 [A+]; B at 300 K not found -> 0 (linear, field independent) [E]
eps_A_STO = 3.030e-3
eps_B_STO = 0.0                      # in m/V (E in V/m); 0 -> eps_r = 330 at every field


def sto_eps_yamamoto(E_V_per_m):
    """eps_r(E) = 1 / (A + B |E|), E in V/m. With B = 0 this is the constant 330."""
    return 1.0 / (eps_A_STO + eps_B_STO * np.abs(E_V_per_m))


# =============================================================================
# Master simulation configuration
# =============================================================================
config = {
    # --- Solver settings (as in the dossier) ---
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
    #   "damped": eps from the previous output field, damped -- original behaviour
    "eps_update": "input",

    # Alpha mixing (adaptive linear mixing, used when mixing_method == "alpha") -- dossier values.
    # Note (dossier, section 7): alpha = 0.05 starts above alpha_max = 0.03; the code only clamps
    # alpha after the first loop, so the first mixing step uses 0.05.
    "alpha": 0.05,
    "alpha_min": 1e-4,
    "alpha_max": 0.03,
    "grow_factor": 1.05,
    "shrink_factor": 0.5,

    # --- Fixed Charge Profiles ---
    # Shape of every fixed charge inside its layer: "delta" | "step" | "gaussian"
    # (polar planes = the layer_charge values below: formal charges, jellium backgrounds, P sheets)
    "charge_profile_defaults": {
        "polar":  {"shape": "delta", "sigma_frac": 0.25},
        "doping": {"shape": "delta", "sigma_frac": 0.25},
        "dipole": {"shape": "delta", "sigma_frac": 0.25},
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
        "material_colors": {
            "STO": "#f4a261",
            "LSMO": "#2a9d8f",
            "YHO_bot": "#6a4c93",
            "YHO": "#9b5de5",
            "YHO_top": "#6a4c93",
            "W": "#8d99ae",
        },
        "background_alpha": 0.35,      # opacity of the layer shading (0 = invisible, 1 = solid)
        "linewidth": 2.0,              # one linewidth for every line in every plot
        "band_ylim_exclude": [],       # materials ignored for the band-plot energy window
        "max_real_gap_eV": 20.0,       # gaps above this are treated as vacuum/barrier and ignored
        "band_ylim": None,             # or (ymin, ymax) in eV to fix the window by hand
        "density_floor": None,         # fixed lower limit in m^-3, or None -> top * density_decades
        "density_decades": 1e-10,      # used when density_floor is None
        "charge_linear_ylim": None,    # charge-distribution left (linear) axis: (ymin, ymax) in m^-3, or None
        "charge_xlim": None,           # charge-distribution plots: (zmin, zmax) in nm to zoom in, or None
    },

    # --- Boundary conditions ---
    "bc_left_type": 'dirichlet',
    "V_left_applied": 0.0,
    "bc_right_type": 'neumann',
    "V_right_applied": 0.0,

    # --- Material definitions ---
    "MATERIALS": {
        "STO": {
            "uc_area": area_STO,
            "layers": [
                {
                    "name": "SrO",
                    # Ec = -chi, chi = 3.9 eV [C]: Hikita APL 90, 143507 (2007), citing Robertson & Chen APL 74, 1168 (1999)
                    # Ev = Ec - 3.25 eV indirect gap [A gap -> C+]: van Benthem et al., JAP 90, 6156 (2001)
                    "Ec": -3.90, "Ev": -7.15,
                    # meff_e 1.3 [B]: Hikita 2007 citing Sroubek PRB 2, 3170 (1970);
                    #   ARPES 1% Nb:STO: 0.63 (light) / 8.0 (heavy) [A] Wakabayashi et al., JAP 138, 115707 (2025)
                    # meff_h 1.2 [E]: no source found (template value)
                    "meff_e": 1.3, "meff_h": 1.2,
                    # eps_r 330 along [001] at 300 K [A]: Neville, Hoeneisen & Mead, JAP 43, 2124 (1972);
                    #   316.3 +- 2.2, isotropic [A]: Goryachev et al., arXiv:1508.07550
                    # eps_r(E) = 1/(A + B|E|): A = 1/330 [A]; B at 300 K not found -> 0 (linear) [E]
                    "eps_r": 330.0, "is_eps_dynamic": True,
                    "eps_dynamic_func": sto_eps_yamamoto,
                    "eps_A": eps_A_STO, "eps_B": eps_B_STO,     # dossier fields (used by sto_eps_yamamoto)
                    "t_layer": t_STO_half,
                    "layer_charge": 0.0           # [A] formal ionic charge of SrO
                },
                {
                    "name": "TiO2",
                    "Ec": -3.90, "Ev": -7.15,
                    "meff_e": 1.3, "meff_h": 1.2,
                    "eps_r": 330.0, "is_eps_dynamic": True,
                    "eps_dynamic_func": sto_eps_yamamoto,
                    "eps_A": eps_A_STO, "eps_B": eps_B_STO,
                    "t_layer": t_STO_half,
                    "layer_charge": 0.0           # [A] formal ionic charge of TiO2
                }
            ]
        },
        "LSMO": {
            "uc_area": area_LSMO,
            "layers": [
                {
                    "name": "LaSrO",
                    # Work function 4.5 eV, midpoint of 4.2-4.8 eV [B]: de Jong et al. (2003), DOI 10.1063/1.1625081
                    # Ec = -WF - E_F, E_F = hbar^2 (3 pi^2 n)^(2/3) / 2m* = 0.46 eV (m* = 2.5, n = 5.65e21 cm^-3) [D]
                    # Ev only keeps hole states empty (1.46 eV below E_F)                                  [E]
                    "Ec": -4.96, "Ev": -5.96,
                    # meff_e 2.5 [B]: Hikita PRB 79, 073101 (2009), citing Okuda PRL 81, 3203 (1998); meff_h unused [E]
                    "meff_e": 2.5, "meff_h": 1.0,
                    # eps_r 30 [C]: Hikita 2009 citing Shannon JAP 73, 348 (1993) (polarizability sum)
                    "eps_r": 30.0, "is_eps_dynamic": False,
                    "t_layer": t_LSMO_half,
                    # formal (La0.67Sr0.33)O = +0.67 [A] plus electron-gas background [D]
                    "layer_charge": +0.67 + jel_LSMO
                },
                {
                    "name": "MnO2",
                    "Ec": -4.96, "Ev": -5.96,
                    "meff_e": 2.5, "meff_h": 1.0,
                    "eps_r": 30.0, "is_eps_dynamic": False,
                    "t_layer": t_LSMO_half,
                    # formal MnO2 = -0.67 [A] plus electron-gas background [D]
                    "layer_charge": -0.67 + jel_LSMO
                }
            ]
        },
        # Y:HfO2 split into bottom / bulk / top so the ferroelectric bound charge sits on the two faces.
        "YHO_bot": {
            "uc_area": area_YHO,
            "layers": [
                {
                    "name": "HfO2_bot",
                    # Ec = -chi, chi = 2.0 +- 0.25 eV [B]: Monaghan et al., Solid-State Electron. 53, 438 (2009)
                    #   (undoped HfO2; CNL-model table value 2.4 eV [C] Robertson & Falabretti 2006) -> sweep 2.0-2.4
                    # Ev = Ec - 6.0 eV gap [C]: literature spread 5.6-6.3 eV
                    "Ec": -2.00, "Ev": -8.00,
                    # meff_e 0.11 +- 0.03 tunnelling mass [B]: Monaghan 2009; meff_h 0.58 [C]: secondary quote only
                    "meff_e": 0.11, "meff_h": 0.58,
                    # eps_r [E]: no usable source (Pavoni 2022 value is the electronic eps1(0), not static) -> sweep 25-35
                    "eps_r": 30.0, "is_eps_dynamic": False,
                    "t_layer": t_YHO,
                    "layer_charge": +sigma_P      # [A+] bound charge at the LSMO face
                }
            ]
        },
        "YHO": {
            "uc_area": area_YHO,
            "layers": [
                {
                    "name": "HfO2",
                    "Ec": -2.00, "Ev": -8.00,
                    "meff_e": 0.11, "meff_h": 0.58,
                    "eps_r": 30.0, "is_eps_dynamic": False,
                    "t_layer": t_YHO,
                    "layer_charge": 0.0           # neutral HfO2 planes
                }
            ]
        },
        "YHO_top": {
            "uc_area": area_YHO,
            "layers": [
                {
                    "name": "HfO2_top",
                    "Ec": -2.00, "Ev": -8.00,
                    "meff_e": 0.11, "meff_h": 0.58,
                    "eps_r": 30.0, "is_eps_dynamic": False,
                    "t_layer": t_YHO,
                    "layer_charge": -sigma_P      # [A+] bound charge at the W face
                }
            ]
        },
        "W": {
            "uc_area": area_W,
            "layers": [
                {
                    "name": "W",
                    # Work function 4.55 eV, polycrystalline [A]: Michaelson JAP 48, 4729 (1977)
                    #   (110)-textured film: 5.05-5.47 eV [B]; beta-W in thin sputtered films is possible
                    # Ec = -WF - 9.17 eV free-electron Fermi energy for 2 e/atom [E]; Ev keeps holes empty [E]
                    "Ec": -13.72, "Ev": -14.72,
                    "meff_e": 1.0, "meff_h": 1.0,   # [E] free-electron placeholder
                    "eps_r": 1.0, "is_eps_dynamic": False,   # [E] background only
                    "t_layer": t_W,
                    "layer_charge": +jel_W          # [E] positive background for 2 e/atom
                }
            ]
        }
    },

    # --- Stack definition (substrate first) ---
    # YHO total 34 layers = 10.0 nm, matching the 9-11 nm films of Yun 2022 [A]; other counts are your choice.
    # Option without W layers: drop "W*8" and set bc_right_type to 'dirichlet' (see dossier, section 7).
    "stack_definition": ["STO*30", "LSMO*20", "YHO_bot*1", "YHO*32", "YHO_top*1", "W*8"],

    "doping_profiles": [],
    "dipole_configs": []
}
