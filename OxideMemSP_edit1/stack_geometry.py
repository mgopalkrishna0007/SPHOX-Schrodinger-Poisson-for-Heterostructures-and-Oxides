# stack_geometry.py
"""
Stack expansion and layer geometry along the growth direction (z).

- expand_stack                       : "STO*100" style stack definition -> flat list of atomic layers
- compute_layer_boundaries           : z positions of every layer boundary
- get_layer_indices_by_offset_range  : interface-relative layer selection (used by doping & dipoles)
- layer_center                       : z position of the centre of a given flat layer
- clamp_dipole_offsets               : delta_plus / delta_minus limited to fit inside a layer
"""

# Elementary charge e = 1.602176487×10^−19 C
E_CHARGE = 1.602176487E-19


def get_layer_indices_by_offset_range(flat_layers, interface_pair, direction, offset_range):
    mat_left, mat_right = interface_pair
    idx = -1
    for i in range(len(flat_layers) - 1):
        if flat_layers[i]["mat_name"] == mat_left and flat_layers[i+1]["mat_name"] == mat_right:
            idx = i
            break
    if idx == -1: return []

    start_off, end_off = offset_range
    indices = []
    if direction == 'right':
        max_off = len(flat_layers) - 1 - idx
        if end_off == -1 or end_off > max_off: end_off = max_off
        for off in range(start_off, end_off + 1): indices.append(idx + off)
    elif direction == 'left':
        max_off = idx + 1
        if end_off == -1 or end_off > max_off: end_off = max_off
        for off in range(start_off, end_off + 1): indices.append(idx - off + 1)
    return indices


def expand_stack(config):
    """--- 1. Expand Stack Definition ---
    Returns flat_layers: every atomic layer of the stack, e.g. ababababcdcdcdcdcd
    """
    flat_layers = []
    MATERIALS = config["MATERIALS"]

    for item in config["stack_definition"]:
        mat_name, count = item.split('*')
        count = int(count)
        mat = MATERIALS[mat_name]
        for _ in range(count):
            for layer in mat["layers"]:  # adds the below extra information to all the layers
                layer_copy = layer.copy()
                layer_copy["mat_name"] = mat_name
                layer_copy["uc_area"] = mat["uc_area"]
                flat_layers.append(layer_copy)

    return flat_layers


def compute_layer_boundaries(flat_layers):
    """--- Pre-calculate Layer Boundaries ---
    Returns (layer_boundaries, z_max).
    """
    layer_boundaries = [0.0]
    curr_z = 0.0
    for L in flat_layers:
        curr_z += L["t_layer"]
        layer_boundaries.append(curr_z)
    z_max = layer_boundaries[-1]
    return layer_boundaries, z_max


def layer_center(flat_layers, layer_boundaries, idx):
    """Exact z position of the centre of flat layer `idx`."""
    return layer_boundaries[idx] + (flat_layers[idx]["t_layer"] / 2.0)


def clamp_dipole_offsets(dip, t_layer):
    """Split-charge offsets of a dipole, scaled down if they do not fit in the layer."""
    d_plus = dip.get("delta_plus_0", 0.0)
    d_minus = dip.get("delta_minus_0", 0.0)

    if (d_plus + d_minus) > t_layer:
        scale = (t_layer * 0.99) / (d_plus + d_minus)
        d_plus *= scale
        d_minus *= scale

    return d_plus, d_minus
