# fixed_charges.py
"""
Fixed (non-mobile) charges placed on the z mesh.

- add_polar_plane_charges : intrinsic layer-based charges ("layer_charge", polar planes)
- add_doping_profiles     : interface-relative donor / acceptor injection ("doping_profiles")
- add_dipole_charges      : split-charge dipole polarization ("dipole_configs")

Every charge can be given one of three spatial profiles inside its layer (see charge_profiles.py):
    "delta"    : on the atomic plane (layer centre) / on the dipole points   -- original behaviour
    "step"     : uniform over the whole layer (dipole: + over the left half, - over the right half)
    "gaussian" : Gaussian centred on the layer centre (dipole: on the + and - points), truncated to the layer
"""
import numpy as np

from stack_geometry import (
    E_CHARGE, get_layer_indices_by_offset_range, layer_center, clamp_dipole_offsets
)
from charge_profiles import owned_edges, profile_weights, resolve_profile, deposit


def _layer_weights(z, edges, shape, sigma, layer_boundaries, idx, z0):
    """Weights of a charge confined to layer idx, centred (delta/gaussian) on z0."""
    return profile_weights(z, edges, shape, layer_boundaries[idx], layer_boundaries[idx + 1],
                           z0=z0, sigma=sigma)


def add_polar_plane_charges(config, flat_layers, layer_boundaries, z, dz_node_arr, N_D_array, N_A_array):
    """--- 4. Intrinsic Layer-Based Charges (Polar Planes) --- (modifies arrays in place)"""
    edges = owned_edges(z)
    for idx, L in enumerate(flat_layers):
        if "layer_charge" in L and L["layer_charge"] != 0.0:
            val_2d = L["layer_charge"] / L["uc_area"]
            z_center = layer_center(flat_layers, layer_boundaries, idx)

            shape, sigma = resolve_profile(config, "polar", L["t_layer"], layer=L)
            w = _layer_weights(z, edges, shape, sigma, layer_boundaries, idx, z_center)

            if L["layer_charge"] > 0:
                deposit(N_D_array, val_2d, w, dz_node_arr)
            else:
                deposit(N_A_array, abs(val_2d), w, dz_node_arr)


def add_doping_profiles(config, flat_layers, layer_boundaries, z, dz_node_arr, N_D_array, N_A_array):
    """--- 5. Interface-Relative Doping Injection (Array Support) --- (modifies arrays in place)"""
    edges = owned_edges(z)
    if "doping_profiles" in config:
        for dop in config["doping_profiles"]:
            target_arr = N_D_array if dop["type"] == "donor" else N_A_array
            val_e_per_uc = dop["value"]

            indices = get_layer_indices_by_offset_range(
                flat_layers, dop["interface"], dop["direction"], dop["layer_offset"]
            )

            for step, idx in enumerate(indices):
                L = flat_layers[idx]
                z_center = layer_center(flat_layers, layer_boundaries, idx)

                shape, sigma = resolve_profile(config, "doping", L["t_layer"], entry=dop)
                w = _layer_weights(z, edges, shape, sigma, layer_boundaries, idx, z_center)

                deposit(target_arr, val_e_per_uc / L["uc_area"], w, dz_node_arr)


def add_dipole_charges(config, flat_layers, layer_boundaries, z, dz_node_arr, info=None):
    """--- 6. Split Charge / Dipole Polarization (Array Support) ---
    Returns (N_D_dipole, N_A_dipole, P_layers).
    info : optional list; one dict per deposited dipole is appended (entry, offset, layer, q, separation,
           dipole moment, node weights) for the diagnostics printed by the solver.

    delta    : +q on z_center - delta_plus, -q on z_center + delta_minus
    step     : +q uniform over [layer start, z_center], -q uniform over [z_center, layer end]
    gaussian : +q Gaussian centred on z_center - delta_plus, -q on z_center + delta_minus,
               both truncated to the layer
    P_layers uses the centroid of the deposited + and - charge, so it follows the chosen shape.
    """
    Nz = len(z)
    edges = owned_edges(z)
    N_D_dipole, N_A_dipole = np.zeros(Nz), np.zeros(Nz)
    P_layers = []

    if "dipole_configs" in config:
        for entry, dip in enumerate(config["dipole_configs"]):
            base_val_e_per_uc = dip["value"]

            indices = get_layer_indices_by_offset_range(
                flat_layers, dip["interface"], dip["direction"], dip["layer_offset"]
            )

            for step, idx in enumerate(indices):
                L = flat_layers[idx]
                z_center = layer_center(flat_layers, layer_boundaries, idx)
                z_lo, z_hi = layer_boundaries[idx], layer_boundaries[idx + 1]

                if dip.get("decay_mode") == 'exponential':
                    decay = np.exp(-step / dip.get("decay_rate", 1.0))
                else:
                    decay = 1.0

                current_val_e_per_uc = base_val_e_per_uc * decay

                d_plus, d_minus = clamp_dipole_offsets(dip, L["t_layer"])

                z_pos = z_center - d_plus
                z_neg = z_center + d_minus

                shape, sigma = resolve_profile(config, "dipole", L["t_layer"], entry=dip)
                if shape == "step":
                    w_pos = profile_weights(z, edges, "step", z_lo, z_center)
                    w_neg = profile_weights(z, edges, "step", z_center, z_hi)
                else:
                    w_pos = profile_weights(z, edges, shape, z_lo, z_hi, z0=z_pos, sigma=sigma)
                    w_neg = profile_weights(z, edges, shape, z_lo, z_hi, z0=z_neg, sigma=sigma)

                # Deposit using local nodal widths
                sheet = current_val_e_per_uc / L["uc_area"]
                deposit(N_D_dipole, sheet, w_pos, dz_node_arr)
                deposit(N_A_dipole, sheet, w_neg, dz_node_arr)

                zc_pos = np.sum(w_pos * z)
                zc_neg = np.sum(w_neg * z)
                dipole_moment = current_val_e_per_uc * E_CHARGE * (zc_pos - zc_neg)
                P_layers.append(dipole_moment / (L["uc_area"] * L["t_layer"]))

                if info is not None:
                    info.append({
                        "entry": entry, "offset": dip["layer_offset"][0] + step, "layer": idx,
                        "name": L["name"], "mat_name": L["mat_name"], "profile": shape,
                        "q_e_per_uc": current_val_e_per_uc, "uc_area": L["uc_area"], "t_layer": L["t_layer"],
                        "z_pos": zc_pos, "z_neg": zc_neg, "separation": zc_neg - zc_pos,
                        "p_Cm": dipole_moment / L["uc_area"],          # dipole moment per area (C/m)
                        "w_pos": w_pos, "w_neg": w_neg,
                    })

    return N_D_dipole, N_A_dipole, P_layers
