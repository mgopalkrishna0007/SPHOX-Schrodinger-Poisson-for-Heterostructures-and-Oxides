# fixed_charges.py
"""
Fixed (non-mobile) charges placed on the z mesh.

- add_polar_plane_charges : intrinsic layer-based charges ("layer_charge", polar planes)
- add_doping_profiles     : interface-relative donor / acceptor injection ("doping_profiles")
- add_dipole_charges      : split-charge dipole polarization ("dipole_configs")
"""
import numpy as np

from stack_geometry import (
    E_CHARGE, get_layer_indices_by_offset_range, layer_center, clamp_dipole_offsets
)


def add_polar_plane_charges(flat_layers, layer_boundaries, z, dz_node_arr, N_D_array, N_A_array):
    """--- 4. Intrinsic Layer-Based Charges (Polar Planes) --- (modifies arrays in place)"""
    for idx, L in enumerate(flat_layers):
        if "layer_charge" in L and L["layer_charge"] != 0.0:
            val_2d = L["layer_charge"] / L["uc_area"]
            z_center = layer_center(flat_layers, layer_boundaries, idx)
            idx_center = np.argmin(np.abs(z - z_center))

            val_3d_delta = val_2d / dz_node_arr[idx_center]

            if L["layer_charge"] > 0:
                N_D_array[idx_center] += val_3d_delta
            else:
                N_A_array[idx_center] += abs(val_3d_delta)


def add_doping_profiles(config, flat_layers, layer_boundaries, z, dz_node_arr, N_D_array, N_A_array):
    """--- 5. Interface-Relative Doping Injection (Array Support) --- (modifies arrays in place)"""
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
                idx_center = np.argmin(np.abs(z - z_center))

                # Apply the exact 3D mesh density strictly for this node
                target_arr[idx_center] += (val_e_per_uc / L["uc_area"]) / dz_node_arr[idx_center]


def add_dipole_charges(config, flat_layers, layer_boundaries, z, dz_node_arr):
    """--- 6. Split Charge / Dipole Polarization (Array Support) ---
    Returns (N_D_dipole, N_A_dipole, P_layers).
    """
    Nz = len(z)
    N_D_dipole, N_A_dipole = np.zeros(Nz), np.zeros(Nz)
    P_layers = []

    if "dipole_configs" in config:
        for dip in config["dipole_configs"]:
            base_val_e_per_uc = dip["value"]

            indices = get_layer_indices_by_offset_range(
                flat_layers, dip["interface"], dip["direction"], dip["layer_offset"]
            )

            for step, idx in enumerate(indices):
                L = flat_layers[idx]
                z_center = layer_center(flat_layers, layer_boundaries, idx)

                if dip.get("decay_mode") == 'exponential':
                    decay = np.exp(-step / dip.get("decay_rate", 1.0))
                else:
                    decay = 1.0

                current_val_e_per_uc = base_val_e_per_uc * decay

                d_plus, d_minus = clamp_dipole_offsets(dip, L["t_layer"])

                z_pos = z_center - d_plus
                z_neg = z_center + d_minus

                mesh_idx_pos = np.argmin(np.abs(z - z_pos))
                mesh_idx_neg = np.argmin(np.abs(z - z_neg))

                # Deposit Deltas using local nodal widths
                N_D_dipole[mesh_idx_pos] += (current_val_e_per_uc / L["uc_area"]) / dz_node_arr[mesh_idx_pos]
                N_A_dipole[mesh_idx_neg] += (current_val_e_per_uc / L["uc_area"]) / dz_node_arr[mesh_idx_neg]

                dipole_moment = current_val_e_per_uc * E_CHARGE * (z_pos - z_neg)
                P_layers.append(dipole_moment / (L["uc_area"] * L["t_layer"]))

    return N_D_dipole, N_A_dipole, P_layers
