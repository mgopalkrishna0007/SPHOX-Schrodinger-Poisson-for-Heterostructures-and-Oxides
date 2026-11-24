# material_properties.py
"""
Assignment of material properties to the mesh nodes.

- build_layer_index_map     : which flat layer every mesh node belongs to
- find_interface_index      : first node where the material changes (tracked in the SC loop)
- assign_material_properties: band edges, effective masses and initial permittivity arrays
- compute_permittivity      : (field-dependent) dielectric constant for a given E-field profile
"""
import numpy as np


def build_layer_index_map(z, layer_boundaries, flat_layers):
    layer_idx_map = np.searchsorted(layer_boundaries, z, side='right') - 1

    # Now every mesh point knows which physical layer it belongs to.
    layer_idx_map[layer_idx_map >= len(flat_layers)] = len(flat_layers) - 1
    return layer_idx_map


def find_interface_index(flat_layers, layer_idx_map):
    """Find the primary interface node for tracking."""
    interface_idx = 0
    Nz = len(layer_idx_map)
    for i in range(Nz - 1):
        if flat_layers[layer_idx_map[i]]["mat_name"] != flat_layers[layer_idx_map[i + 1]]["mat_name"]:
            # "Does the material change between these two neighbouring nodes?"
            interface_idx = i
            break
    return interface_idx


def compute_permittivity(flat_layers, layer_idx_map, E_Vm):
    """Relative permittivity on every node.

    Layers flagged with "is_eps_dynamic" and providing "eps_dynamic_func" use the
    field-dependent function from the config (the solver doesn't know the physics,
    it just calls the function); all others use their constant "eps_r".
    E_Vm : absolute electric field (V/m) on every node.
    """
    Nz = len(layer_idx_map)
    Epsi = np.zeros(Nz)
    for k in range(Nz):
        L = flat_layers[layer_idx_map[k]]
        if L.get("is_eps_dynamic", False) and "eps_dynamic_func" in L:
            Epsi[k] = L["eps_dynamic_func"](E_Vm[k])
        else:
            Epsi[k] = L["eps_r"]
    return Epsi


def assign_material_properties(flat_layers, layer_idx_map):
    """--- 3. Property Array Initialization + Fill material properties ---
    Returns dict with V0_c, V0_v, meff_e_arr, meff_h_arr, Epsi_arr.
    The dielectric constant is initialized with an electric field of 0.0 V/m.
    """
    Nz = len(layer_idx_map)
    V0_c, V0_v = np.zeros(Nz), np.zeros(Nz)
    meff_e_arr, meff_h_arr = np.zeros(Nz), np.zeros(Nz)

    for k in range(Nz):
        L = flat_layers[layer_idx_map[k]]
        V0_c[k] = L["Ec"]
        V0_v[k] = L["Ev"]
        meff_e_arr[k] = L["meff_e"]
        meff_h_arr[k] = L["meff_h"]

    Epsi_arr = compute_permittivity(flat_layers, layer_idx_map, np.zeros(Nz))

    return {
        "V0_c": V0_c, "V0_v": V0_v,
        "meff_e_arr": meff_e_arr, "meff_h_arr": meff_h_arr,
        "Epsi_arr": Epsi_arr,
    }
