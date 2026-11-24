# charge_profiles.py
"""
Spatial profile of a fixed sheet charge inside its layer.

Every fixed charge (polar plane, doping, dipole half) is a sheet charge sigma_2D = value / uc_area
(e/m^2) that has to be written onto the mesh nodes as a volume density (e/m^3). Three shapes:

    "delta"    : the whole charge on ONE node (the atomic plane / dipole point)
                     N[k] += sigma_2D / dz_node[k]
    "step"     : uniform density sigma_2D / width over an interval [z_lo, z_hi] (e.g. the whole layer)
    "gaussian" : Gaussian centred on z0 with width sigma, truncated to [z_lo, z_hi] and renormalised
                 so that no charge leaks out of the layer

How step/gaussian are put on the nodes (control-volume / "owned width" picture):
    node k owns [edge_left_k, edge_right_k] (halfway to its neighbours, length dz_node[k]).
    weight_k = (integral of the normalised profile over the part of node k's owned stretch that
                lies inside [z_lo, z_hi])
    N[k]    += sigma_2D * weight_k / dz_node[k]
The weights sum to 1, so sum_k N[k]*dz_node[k] (= np.trapz(N, z)) is exactly sigma_2D, and only
nodes whose owned stretch overlaps [z_lo, z_hi] receive charge. A layer boundary node owns half
in each layer, so a step profile gives it half the plateau value (no double counting when two
neighbouring layers are both stepped).
"""
import numpy as np
from scipy.special import erf

PROFILE_SHAPES = ("delta", "step", "gaussian")
DEFAULT_SIGMA_FRAC = 0.25   # gaussian sigma as a fraction of t_layer if nothing else is given


def owned_edges(z):
    """Left / right edge of the stretch of z owned by every node (they tile [z0, z_end])."""
    mid = 0.5 * (z[1:] + z[:-1])
    left = np.concatenate([[z[0]], mid])
    right = np.concatenate([mid, [z[-1]]])
    return left, right


def profile_weights(z, edges, shape, z_lo, z_hi, z0=None, sigma=None):
    """Fraction of the sheet charge that goes to every node (sums to 1).

    shape  : "delta" | "step" | "gaussian"
    z_lo, z_hi : interval the charge is confined to (step: spread over it; gaussian: truncated to it)
    z0     : delta position / gaussian centre
    sigma  : gaussian standard deviation (m)
    """
    w = np.zeros(len(z))

    if shape == "delta":
        w[np.argmin(np.abs(z - z0))] = 1.0
        return w

    left, right = edges
    a = np.clip(left, z_lo, z_hi)          # part of each owned stretch inside [z_lo, z_hi]
    b = np.clip(right, z_lo, z_hi)

    if shape == "step":
        w = (b - a) / (z_hi - z_lo)
    elif shape == "gaussian":
        s = np.sqrt(2.0) * sigma
        cdf = lambda x: 0.5 * (1.0 + erf((x - z0) / s))
        w = cdf(b) - cdf(a)
        total = cdf(z_hi) - cdf(z_lo)
        w = w / total
    else:
        raise ValueError(f"Unknown charge profile '{shape}'. Choose one of: {', '.join(PROFILE_SHAPES)}")

    w[w < 0.0] = 0.0
    return w / np.sum(w)


def resolve_profile(config, kind, t_layer, entry=None, layer=None):
    """Pick shape and gaussian sigma for one charge.

    Priority (highest first):
        entry["profile"]          (a doping_profiles / dipole_configs entry)
        layer["charge_profile"]   (a layer dict in MATERIALS, for polar planes)
        config["charge_profile_defaults"][kind]["shape"]
        "delta"
    sigma: entry/layer "sigma" (m) > entry/layer "sigma_frac" * t_layer > defaults[kind]["sigma_frac"] * t_layer
    kind : "polar" | "doping" | "dipole"
    """
    defaults = config.get("charge_profile_defaults", {}).get(kind, {})
    src = entry if entry is not None else (layer if layer is not None else {})
    shape_key = "profile" if entry is not None else "charge_profile"

    shape = str(src.get(shape_key, defaults.get("shape", "delta"))).lower()
    if shape not in PROFILE_SHAPES:
        raise ValueError(f"Unknown charge profile '{shape}' ({kind}). Choose one of: {', '.join(PROFILE_SHAPES)}")

    if "sigma" in src:
        sigma = src["sigma"]
    else:
        sigma = src.get("sigma_frac", defaults.get("sigma_frac", DEFAULT_SIGMA_FRAC)) * t_layer

    return shape, sigma


def deposit(target, sheet_density, weights, dz_node_arr):
    """Add a sheet charge (e/m^2) to a nodal volume-density array (e/m^3) in place."""
    nz = np.nonzero(weights)[0]
    target[nz] += sheet_density * weights[nz] / dz_node_arr[nz]
