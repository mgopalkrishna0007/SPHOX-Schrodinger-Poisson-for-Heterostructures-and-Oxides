# mesh_generation.py
"""
1D mesh generation in the z domain.

- collect_exact_points : layer boundaries, layer centres and dipole split-charge positions
- generate_mesh        : uniform mesh + tolerance-based merging with the exact points
- compute_nodal_widths : effective nodal width of every node (for 2D -> 3D density mapping)
"""
import numpy as np

from stack_geometry import (
    get_layer_indices_by_offset_range, layer_center, clamp_dipole_offsets
)


def collect_exact_points(config, flat_layers, layer_boundaries, z_max):
    """--- 2. Exact Mesh Generation ---
    Collect exact geometric points to prevent snapping offsets.
    """
    exact_z = []
    for i, L in enumerate(flat_layers):
        exact_z.append(layer_boundaries[i])                       # Exact boundary
        exact_z.append(layer_boundaries[i] + (L["t_layer"]/2.0))  # Exact center
    exact_z.append(z_max)

    # Pre-calculate explicit dipole offset points
    if "dipole_configs" in config:
        for dip in config["dipole_configs"]:
            indices = get_layer_indices_by_offset_range(
                flat_layers, dip["interface"], dip["direction"], dip["layer_offset"]
            )
            for idx in indices:
                L = flat_layers[idx]
                z_center = layer_center(flat_layers, layer_boundaries, idx)
                d_plus, d_minus = clamp_dipole_offsets(dip, L["t_layer"])

                exact_z.append(z_center - d_plus)
                exact_z.append(z_center + d_minus)

    return exact_z


def generate_mesh(config, flat_layers, exact_z, z_max, tol=1e-15):
    """Uniform base mesh merged with the exact geometric points.

    tol : uniform points closer than this (meters) to an exact point are dropped.
    Returns the sorted node array z.
    """
    # Generate base uniform mesh
    min_t = min(L["t_layer"] for L in flat_layers)
    dz_fine = min_t / config["points_per_layer"]
    Nz_uniform = int(round(z_max / dz_fine)) + 1
    z_uniform = np.linspace(0, z_max, Nz_uniform)

    # --- TOLERANCE-BASED MESH MERGING ---
    # Start with the exact geometric points
    z_exact = np.unique(np.array(exact_z))
    z_exact.sort()

    # Add uniform points only if they are not already sufficiently close to an exact point
    z_list = list(z_exact)

    for zu in z_uniform:
        if np.min(np.abs(z_exact - zu)) > tol:
            z_list.append(zu)

    z = np.sort(np.array(z_list))
    Nz = len(z)

    # Check minimum element size
    dz = np.diff(z)
    print(f"Mesh generated with {Nz} nodes")
    print(f"Minimum dz = {np.min(dz):.6e} m")
    print(f"Maximum dz = {np.max(dz):.6e} m")
    print(f"Target dz  = {dz_fine:.6e} m")

    return z


def compute_nodal_widths(z):
    """Local effective nodal volumes for 1D density mapping."""
    Nz = len(z)
    dz_node_arr = np.zeros(Nz)
    dz_node_arr[0] = (z[1] - z[0]) / 2.0
    dz_node_arr[-1] = (z[-1] - z[-2]) / 2.0
    dz_node_arr[1:-1] = (z[2:] - z[:-2]) / 2.0
    return dz_node_arr
