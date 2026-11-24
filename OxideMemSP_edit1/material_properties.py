# material_properties.py
"""
Assignment of material properties to the mesh nodes.

- build_layer_index_map     : which flat layer every mesh node belongs to
- find_interface_index      : first node where the material changes (tracked in the SC loop)
- assign_material_properties: band edges, effective masses and initial permittivity arrays
- compute_permittivity      : (field-dependent) dielectric constant for a given E-field profile
- build_eps_frozen_mask     : nodes that keep a constant eps (dipole layers by default)
- check_eps_saturation      : warns if a fixed charge needs more D than eps(E) can provide
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


def compute_permittivity(flat_layers, layer_idx_map, E_Vm, eps_frozen=None):
    """Relative permittivity on every node.

    Layers flagged with "is_eps_dynamic" and providing "eps_dynamic_func" use the
    field-dependent function from the config (the solver doesn't know the physics,
    it just calls the function); all others use their constant "eps_r".
    E_Vm       : absolute electric field (V/m) on every node.
    eps_frozen : optional boolean mask; nodes set True always use the constant "eps_r"
                 (used inside dipole layers, see build_eps_frozen_mask).
    """
    Nz = len(layer_idx_map)
    Epsi = np.zeros(Nz)
    for k in range(Nz):
        L = flat_layers[layer_idx_map[k]]
        frozen = eps_frozen is not None and eps_frozen[k]
        if not frozen and L.get("is_eps_dynamic", False) and "eps_dynamic_func" in L:
            Epsi[k] = L["eps_dynamic_func"](E_Vm[k])
        else:
            Epsi[k] = L["eps_r"]
    return Epsi


def build_eps_frozen_mask(config, flat_layers, layer_idx_map):
    """Nodes whose permittivity must NOT follow the field-dependent eps(E).

    Inside a dipole layer the displacement field between the +q and -q sheets is
    D = e * value / uc_area, which for a few e/uc is far above what a saturating eps(E)
    can carry (see check_eps_saturation): eps(E) -> 0 and the SC loop diverges.
    By default every dipole layer therefore keeps its constant "eps_r".
    Set "eps_dynamic_inside": True in a dipole_configs entry to let eps(E) act there.
    """
    from stack_geometry import get_layer_indices_by_offset_range

    frozen_layers = set()
    for dip in config.get("dipole_configs", []):
        if not dip.get("eps_dynamic_inside", False):
            frozen_layers.update(get_layer_indices_by_offset_range(
                flat_layers, dip["interface"], dip["direction"], dip["layer_offset"]))
    return np.isin(layer_idx_map, list(frozen_layers))


def max_displacement(eps_func, E_max_MV_cm=1e4):
    """Largest |D| = eps(E) * eps0 * E (C/m^2) a field-dependent eps can produce (scanned up to E_max)."""
    eps0 = 8.854187817620E-12
    E = np.logspace(-4, np.log10(E_max_MV_cm), 3000) * 1e8      # V/m
    return max(eps_func(x) * eps0 * x for x in E)


def check_eps_saturation(config, flat_layers, layer_idx_map, eps_frozen):
    """Warn about fixed sheet charges in field-dependent layers that need a larger D than eps(E) allows.

    A sheet of sigma = value/uc_area (e/m^2) needs |D| up to e*sigma next to it. If that exceeds the
    maximum D of the layer's eps(E) model, Poisson + eps(E) has no solution and the loop diverges.
    Returns the list of warning strings (also printed).
    """
    from stack_geometry import get_layer_indices_by_offset_range
    e = 1.602176487E-19
    in_frozen = set(np.unique(layer_idx_map[eps_frozen])) if eps_frozen is not None else set()
    dmax_cache, warnings_out = {}, []

    def dmax(L):
        f = L["eps_dynamic_func"]
        if f not in dmax_cache:
            dmax_cache[f] = max_displacement(f)
        return dmax_cache[f]

    def check(idx, value, what):
        L = flat_layers[idx]
        if idx in in_frozen or not (L.get("is_eps_dynamic", False) and "eps_dynamic_func" in L):
            return
        D_need, D_max = abs(value) / L["uc_area"] * e, dmax(L)
        if D_need > D_max:
            warnings_out.append(
                f"[eps(E) WARNING] {what} in layer {idx} ({L['mat_name']}:{L['name']}) needs |D| = {D_need:.3f} C/m^2 "
                f"> max D = {D_max:.3f} C/m^2 of its eps(E) -> no self-consistent solution (eps -> 0, loop diverges). "
                f"Freeze eps there or reduce the charge.")

    for idx, L in enumerate(flat_layers):
        if L.get("layer_charge", 0.0) != 0.0:
            check(idx, L["layer_charge"], f"polar charge {L['layer_charge']:+g} e/uc")
    for dop in config.get("doping_profiles", []):
        for idx in get_layer_indices_by_offset_range(flat_layers, dop["interface"], dop["direction"], dop["layer_offset"]):
            check(idx, dop["value"], f"doping {dop['value']:+g} e/uc")
    for dip in config.get("dipole_configs", []):
        for idx in get_layer_indices_by_offset_range(flat_layers, dip["interface"], dip["direction"], dip["layer_offset"]):
            check(idx, dip["value"], f"dipole {dip['value']:g} e/uc")

    for w in warnings_out:
        print(w)
    return warnings_out


def assign_material_properties(flat_layers, layer_idx_map, eps_frozen=None):
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

    Epsi_arr = compute_permittivity(flat_layers, layer_idx_map, np.zeros(Nz), eps_frozen)

    return {
        "V0_c": V0_c, "V0_v": V0_v,
        "meff_e_arr": meff_e_arr, "meff_h_arr": meff_h_arr,
        "Epsi_arr": Epsi_arr,
    }
