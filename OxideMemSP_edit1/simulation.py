# simulation.py
"""
run_heterostructure_simulation(config): the full single-run pipeline.

    build_device(config)            -> stack, mesh, material properties, fixed charges
    run_self_consistent_loop(...)   -> self-consistent Poisson-Schrödinger solution
"""
import numpy as np

from stack_geometry import expand_stack, compute_layer_boundaries
from mesh_generation import collect_exact_points, generate_mesh, compute_nodal_widths
from material_properties import (
    build_layer_index_map, find_interface_index, assign_material_properties,
    build_eps_frozen_mask, check_eps_saturation
)
from fixed_charges import add_polar_plane_charges, add_doping_profiles, add_dipole_charges
from scf_solver import run_self_consistent_loop


def build_device(config):
    """Steps 1-6: everything assigned on the z domain before the SC loop starts."""
    # --- 1. Expand Stack Definition + Layer Boundaries ---
    flat_layers = expand_stack(config)
    layer_boundaries, z_max = compute_layer_boundaries(flat_layers)

    # --- 2. Exact Mesh Generation ---
    exact_z = collect_exact_points(config, flat_layers, layer_boundaries, z_max)
    z = generate_mesh(config, flat_layers, exact_z, z_max)
    dz_node_arr = compute_nodal_widths(z)

    # --- 3. Property Array Initialization / Fill material properties ---
    layer_idx_map = build_layer_index_map(z, layer_boundaries, flat_layers)
    interface_idx = find_interface_index(flat_layers, layer_idx_map)
    eps_frozen = build_eps_frozen_mask(config, flat_layers, layer_idx_map)
    check_eps_saturation(config, flat_layers, layer_idx_map, eps_frozen)
    props = assign_material_properties(flat_layers, layer_idx_map, eps_frozen)

    # --- 4. Intrinsic Layer-Based Charges (Polar Planes) ---
    N_D_array, N_A_array = np.zeros(len(z)), np.zeros(len(z))
    add_polar_plane_charges(config, flat_layers, layer_boundaries, z, dz_node_arr, N_D_array, N_A_array)

    # --- 5. Interface-Relative Doping Injection ---
    add_doping_profiles(config, flat_layers, layer_boundaries, z, dz_node_arr, N_D_array, N_A_array)

    # --- 6. Split Charge / Dipole Polarization ---
    N_D_dipole, N_A_dipole, P_layers = add_dipole_charges(config, flat_layers, layer_boundaries, z, dz_node_arr)

    return {
        "flat_layers": flat_layers,
        "layer_boundaries": layer_boundaries,
        "z": z,
        "dz_node_arr": dz_node_arr,
        "layer_idx_map": layer_idx_map,
        "interface_idx": interface_idx,
        "eps_frozen": eps_frozen,
        **props,
        "N_D_array": N_D_array, "N_A_array": N_A_array,
        "N_D_dipole": N_D_dipole, "N_A_dipole": N_A_dipole,
        "P_layers": P_layers,
    }


def run_heterostructure_simulation(config):
    device = build_device(config)
    # --- 7./8. FEniCS initialization + Self-Consistent Loop ---
    return run_self_consistent_loop(config, device)
