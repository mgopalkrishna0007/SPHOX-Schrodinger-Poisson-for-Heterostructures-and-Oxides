# diagnostics.py
"""
Printed diagnostics for the heterostructure solver.

1. verify_device(...)        -- once, before the SC loop: is the stack / mesh / layer map / every fixed
                                charge (in e/uc) assigned correctly?
2. polarization_report(...)  -- every loop: inserted dipoles and the polarization (dipole-moment) balance
                                that shows how they are screened.

Polarization balance (all quantities are dipole moments per unit area, C/m):

    Poisson solved by the code:   d/dz( eps0*eps*E ) = e*rho ,   E = dU/dz   (U = electron potential energy / e)

    Multiply by (z - z_ref) and integrate over the device (z_ref = the Dirichlet side, where the potential
    is fixed; on the other, Neumann side D = 0). Integrating by parts gives an exact identity:

        p_inserted + p_fixed + p_mobile  +  P_lattice  +  eps0*dU_total  =  0

        p_X        = e * integral( rho_X * (z - z_ref) dz )     explicit charges (inserted dipoles,
                                                                 fixed ions/doping, electrons+holes)
        P_lattice  = integral( eps0*(eps-1)*E dz )              polarization of the dielectric itself
                                                                 (the part of STO's eps = 300 response)
        eps0*dU    = eps0 * (U_end - U_start)                    what is left over: shows up as a net
                                                                 potential step across the whole stack

    "Screening" = everything that is not the inserted dipole (charges that moved or were placed
    elsewhere + lattice polarization). Net = inserted + screening = -eps0*dU, i.e. only the
    remainder that actually shifts the bands. The identity is evaluated with the same
    piecewise-linear integrals FEniCS assembles, so the printed residual should be ~1e-12 or smaller.
"""
import numpy as np

from stack_geometry import E_CHARGE, get_layer_indices_by_offset_range
from charge_profiles import owned_edges, profile_weights, resolve_profile

EPS0 = 8.854187817620E-12
ANG2 = 1e-20   # m^2 per square angstrom


# =============================================================================
# Helpers
# =============================================================================
def fe_moment(rho, z, z_ref):
    """Exact integral of rho_h(z) * (z - z_ref) for piecewise-linear rho_h (as assembled by FEniCS)."""
    g = z - z_ref
    h = np.diff(z)
    fa, fb, ga, gb = rho[:-1], rho[1:], g[:-1], g[1:]
    return np.sum(h / 6.0 * (2 * fa * ga + fa * gb + fb * ga + 2 * fb * gb))


def reference_side(config, z):
    """Dirichlet side is the natural origin of the dipole moments (potential fixed there)."""
    if str(config.get("bc_left_type", "neumann")).lower() == "dirichlet":
        return z[0], "left"
    if str(config.get("bc_right_type", "neumann")).lower() == "dirichlet":
        return z[-1], "right"
    return z[0], "left"


def reference_area(device, dipole_info):
    """In-plane unit-cell area used to express moments in e*nm/uc (host of the first dipole,
    otherwise the material right after the first interface)."""
    if dipole_info:
        return dipole_info[0]["uc_area"], dipole_info[0]["mat_name"]
    fl, lmap = device["flat_layers"], device["layer_idx_map"]
    L = fl[lmap[min(device["interface_idx"] + 1, len(lmap) - 1)]]
    return L["uc_area"], L["mat_name"]


def uc_thickness(config, mat_name):
    return sum(L["t_layer"] for L in config["MATERIALS"][mat_name]["layers"])


def e_nm_per_uc(p_Cm, area):
    """Dipole moment per area (C/m) -> e*nm per in-plane unit cell."""
    return p_Cm / E_CHARGE * area * 1e9


# =============================================================================
# 1. Device check (printed once before the SC loop)
# =============================================================================
def verify_device(config, device, N_D_polar, N_A_polar, dipole_info):
    z = device["z"]; fl = device["flat_layers"]; lmap = device["layer_idx_map"]
    B = np.asarray(device["layer_boundaries"]); dzn = device["dz_node_arr"]
    edges = owned_edges(z)
    M = config["MATERIALS"]
    ok_all = True
    lines = []

    def flag(ok):
        nonlocal ok_all
        ok_all &= bool(ok)
        return "OK" if ok else "<-- CHECK"

    # --- Stack ---
    exp_layers, exp_T, parts = 0, 0.0, []
    for item in config["stack_definition"]:
        name, count = item.split('*'); count = int(count)
        n = count * len(M[name]["layers"]); T = count * sum(L["t_layer"] for L in M[name]["layers"])
        exp_layers += n; exp_T += T
        parts.append(f"{item}: {n} layers, {T*1e9:.3f} nm")
    lines.append(f" [Stack ] {len(fl)} atomic layers built / {exp_layers} expected, thickness "
                 f"{B[-1]*1e9:.4f} / {exp_T*1e9:.4f} nm  {flag(len(fl) == exp_layers and abs(B[-1]-exp_T) < 1e-15)}")
    lines.append("          " + " | ".join(parts))

    # --- Mesh / points per layer ---
    min_t = min(L["t_layer"] for L in fl)
    dz_target = min_t / config["points_per_layer"]
    dz = np.diff(z)
    elems = np.array([np.count_nonzero((z >= B[i] - 1e-18) & (z <= B[i + 1] + 1e-18)) - 1 for i in range(len(fl))])
    imin = int(np.argmin(elems))
    slivers = int(np.count_nonzero(dz < 0.05 * dz_target))
    lines.append(f" [Mesh  ] {len(z)} nodes | dz min {dz.min()*1e9:.2e} nm, max {dz.max()*1e9:.4f} nm, target "
                 f"{dz_target*1e9:.4f} nm (thinnest layer {min_t*1e9:.4f} nm / points_per_layer {config['points_per_layer']})")
    lines.append(f"          elements per atomic layer: min {elems.min()} ({fl[imin]['mat_name']}:{fl[imin]['name']}), "
                 f"max {elems.max()}  {flag(elems.min() >= config['points_per_layer'] - 1)}"
                 + (f" | {slivers} sliver element(s) < 5% of target dz (delta spikes there look taller; harmless)" if slivers else ""))

    # --- Exact points ---
    centres = B[:-1] + np.array([L["t_layer"] for L in fl]) / 2
    gap = lambda pts: max((np.min(np.abs(z - p)) for p in pts), default=0.0)
    dip_pts = [d["z_pos"] for d in dipole_info if d["profile"] == "delta"] + \
              [d["z_neg"] for d in dipole_info if d["profile"] == "delta"]
    g_b, g_c, g_d = gap(B), gap(centres), gap(dip_pts)
    lines.append(f" [Points] largest distance to the nearest node: boundaries {g_b:.1e} m, layer centres {g_c:.1e} m"
                 + (f", delta-dipole charges {g_d:.1e} m" if dip_pts else "") + f"  {flag(max(g_b, g_c, g_d) == 0.0)}")

    # --- Layer map / properties ---
    inside = np.all((B[lmap] <= z) & ((z < B[lmap + 1]) | (lmap == len(fl) - 1)))
    roundtrip = all(lmap[np.argmin(np.abs(z - centres[i]))] == i for i in range(len(fl)))
    props = all(device["V0_c"][k] == fl[lmap[k]]["Ec"] and device["meff_e_arr"][k] == fl[lmap[k]]["meff_e"]
                for k in range(len(z)))
    lines.append(f" [Layers] node inside its layer: {inside} | centre node -> same layer: {roundtrip} | "
                 f"Ec, m* on nodes = layer values: {props}  {flag(inside and roundtrip and props)}")

    # --- Fixed charges, one by one, in e/uc ---
    lines.append(" [Charge] e/uc = electrons per in-plane unit cell of the HOST layer (sheet density = value / uc_area)")
    rebuilt_D, rebuilt_A = np.zeros(len(z)), np.zeros(len(z))
    expected_sheet = 0.0   # e/m^2, signed

    def contained(w, i):
        nz = np.nonzero(w)[0]
        return bool(np.all((z[nz] >= B[i] - 1e-18) & (z[nz] <= B[i + 1] + 1e-18)))

    # polar planes (grouped per material:layer name)
    groups = {}
    for i, L in enumerate(fl):
        q = L.get("layer_charge", 0.0)
        if q == 0.0:
            continue
        shape, sig = resolve_profile(config, "polar", L["t_layer"], layer=L)
        w = profile_weights(z, edges, shape, B[i], B[i + 1], z0=centres[i], sigma=sig)
        (rebuilt_D if q > 0 else rebuilt_A)[:] += abs(q) / L["uc_area"] * w / dzn
        expected_sheet += q / L["uc_area"]
        g = groups.setdefault((L["mat_name"], L["name"], q, shape),
                              {"n": 0, "dep": [], "in": True, "area": L["uc_area"]})
        g["n"] += 1; g["dep"].append(q * np.sum(w)); g["in"] &= contained(w, i)
    if not groups:
        lines.append("   polar planes : none")
    for (mat, name, q, shape), g in groups.items():
        dep = np.array(g["dep"])
        lines.append(f"   polar planes : {mat}:{name} {q:+.3f} e/uc x {g['n']} layers ({shape}), area {g['area']/ANG2:.2f} A^2 "
                     f"-> {q/g['area']*1e-4:+.3e} cm^-2 each | deposited {dep.min():+.6f}..{dep.max():+.6f} e/uc, "
                     f"inside layer {flag(g['in'] and np.allclose(dep, q, rtol=1e-12, atol=1e-15))}")

    # doping
    for j, dop in enumerate(config.get("doping_profiles", [])):
        idxs = get_layer_indices_by_offset_range(fl, dop["interface"], dop["direction"], dop["layer_offset"])
        if not idxs:
            lines.append(f"   doping #{j+1}    : interface {tuple(dop['interface'])} not found in the stack -> NOT placed  {flag(dop['value'] == 0)}")
            continue
        for step, i in enumerate(idxs):
            L = fl[i]
            shape, sig = resolve_profile(config, "doping", L["t_layer"], entry=dop)
            w = profile_weights(z, edges, shape, B[i], B[i + 1], z0=centres[i], sigma=sig)
            target = rebuilt_D if dop["type"] == "donor" else rebuilt_A
            target += dop["value"] / L["uc_area"] * w / dzn
            sgn = 1 if dop["type"] == "donor" else -1
            expected_sheet += sgn * dop["value"] / L["uc_area"]
            dep = dop["value"] * np.sum(w)
            lines.append(f"   doping #{j+1}    : {dop['type']:8s} {dop['value']:+.3f} e/uc, offset {dop['layer_offset'][0]+step} "
                         f"({dop['direction']}) -> layer {i} {L['mat_name']}:{L['name']} ({shape}), area {L['uc_area']/ANG2:.2f} A^2 "
                         f"-> {sgn*dop['value']/L['uc_area']*1e-4:+.3e} cm^-2 | deposited {dep:+.6f} e/uc on "
                         f"{np.count_nonzero(w)} node(s), inside layer {flag(contained(w, i) and abs(dep - dop['value']) < 1e-12)}")

    # dipoles
    rebuilt_dD, rebuilt_dA = np.zeros(len(z)), np.zeros(len(z))
    for d in dipole_info:
        i = d["layer"]
        rebuilt_dD += d["q_e_per_uc"] / d["uc_area"] * d["w_pos"] / dzn
        rebuilt_dA += d["q_e_per_uc"] / d["uc_area"] * d["w_neg"] / dzn
        qp, qm = d["q_e_per_uc"] * np.sum(d["w_pos"]), d["q_e_per_uc"] * np.sum(d["w_neg"])
        ok = contained(d["w_pos"], i) and contained(d["w_neg"], i) and abs(qp - qm) < 1e-12
        lines.append(f"   dipole #{d['entry']+1}    : +/-{d['q_e_per_uc']:.3f} e/uc, offset {d['offset']} -> layer {i} "
                     f"{d['mat_name']}:{d['name']} ({d['profile']}), separation {d['separation']*1e9:.4f} nm | "
                     f"+q {qp:.6f} / -q {qm:.6f} e/uc, net 0, both inside layer {flag(ok)}")
    for dip in config.get("dipole_configs", []):
        if not get_layer_indices_by_offset_range(fl, dip["interface"], dip["direction"], dip["layer_offset"]):
            lines.append(f"   dipole       : interface {tuple(dip['interface'])} not found in the stack -> NOT placed  {flag(False)}")

    # arrays == sum of the individual charges
    dev = lambda a, b: np.max(np.abs(a - b)) / max(np.max(np.abs(b)), 1.0)
    err = max(dev(rebuilt_D, device["N_D_array"]), dev(rebuilt_A, device["N_A_array"]),
              dev(rebuilt_dD, device["N_D_dipole"]), dev(rebuilt_dA, device["N_A_dipole"]))
    lines.append(f"   charge arrays rebuilt from the individual charges above: max relative difference {err:.1e}  {flag(err < 1e-12)}")

    # totals and what neutrality requires
    got_sheet = np.sum((device["N_D_array"] - device["N_A_array"]) * dzn)
    A_ref, mat_ref = reference_area(device, dipole_info)
    lines.append(f"   net fixed charge (excluding dipoles): expected {expected_sheet:+.4e} m^-2, on mesh {got_sheet:+.4e} m^-2 "
                 f"{flag(abs(got_sheet - expected_sheet) <= 1e-9 * max(abs(expected_sheet), 1e10))}")
    lines.append(f"   -> charge neutrality needs a net mobile charge of {-expected_sheet:+.4e} m^-2 "
                 f"= {-expected_sheet*A_ref:+.4f} e/uc of {mat_ref} (electrons if negative)")

    # unit-cell areas used for e/uc conversions
    areas = []
    for name in dict.fromkeys(L["mat_name"] for L in fl):
        a = M[name]["uc_area"] / ANG2
        areas.append(f"{name} {a:.2f} A^2" + (" (unusually small/large - check)" if not 5.0 <= a <= 500.0 else ""))
    lines.append(" [Areas ] uc_area per material (only used to convert e/uc): " + ", ".join(areas))

    print("==========================================================================================")
    print(" DEVICE CHECK: stack, mesh, layer map and fixed charges")
    print("------------------------------------------------------------------------------------------")
    for l in lines:
        print(l)
    print(f" Overall: {'all checks passed' if ok_all else 'SOME CHECKS NEED ATTENTION (see <-- CHECK)'}")
    print("==========================================================================================\n")
    return ok_all


# =============================================================================
# 2. Polarization / screening balance (printed every loop)
# =============================================================================
def polarization_report(config, device, dipole_info, U, eps_used, ND_plus_z, NA_minus_z, N_3D, P_3D):
    """Return (lines, numbers). U = Vout (electron potential energy in eV), eps_used = eps given to Poisson."""
    z = device["z"]; B = np.asarray(device["layer_boundaries"])
    z_ref, side = reference_side(config, z)
    A_ref, mat_ref = reference_area(device, dipole_info)

    p_ins = E_CHARGE * fe_moment(device["N_D_dipole"] - device["N_A_dipole"], z, z_ref)
    p_fix = E_CHARGE * fe_moment(ND_plus_z - NA_minus_z, z, z_ref)
    p_mob = E_CHARGE * fe_moment(P_3D - N_3D, z, z_ref)

    dU = np.diff(U)
    eps_el = 0.5 * (eps_used[1:] + eps_used[:-1])
    p_latt = EPS0 * np.sum((eps_el - 1.0) * dU)
    p_vac = EPS0 * (U[-1] - U[0])
    residual = p_ins + p_fix + p_mob + p_latt + p_vac
    scale = max(abs(p_ins), abs(p_fix), abs(p_mob), abs(p_latt), 1e-30)

    f = lambda p: e_nm_per_uc(p, A_ref)
    lines = []
    # --- one line per inserted dipole ---
    dzn = device["dz_node_arr"]
    for d in dipole_info:
        i = d["layer"]
        c_uc = uc_thickness(config, d["mat_name"])
        P_layer = d["p_Cm"] / d["t_layer"] * 100.0          # C/m^2 -> uC/cm^2
        P_uc = d["p_Cm"] / c_uc * 100.0
        # the same dipole as Poisson sees it (piecewise-linear charge; differs slightly from the nominal
        # value when a delta charge sits next to a very thin element)
        sheet = d["q_e_per_uc"] / d["uc_area"]
        p_fe = E_CHARGE * fe_moment(sheet * (d["w_pos"] - d["w_neg"]) / dzn, z, z_ref)
        inl = (z[:-1] >= B[i] - 1e-18) & (z[1:] <= B[i + 1] + 1e-18)        # elements inside the host layer
        eps_layer = np.sum(eps_el[inl] * np.diff(z)[inl]) / d["t_layer"]
        step_own = -p_fe / (EPS0 * eps_layer)                               # step from this dipole alone
        k0, k1 = np.argmin(np.abs(z - B[i])), np.argmin(np.abs(z - B[i + 1]))
        lines.append(f" Dipole #{d['entry']+1} | offset {d['offset']} -> {d['mat_name']}:{d['name']} | q = {d['q_e_per_uc']:.3f} e/uc | "
                     f"d = {d['separation']*1e9:.4f} nm | p = {e_nm_per_uc(d['p_Cm'], d['uc_area']):+.4f} e*nm/uc "
                     f"(on mesh {e_nm_per_uc(p_fe, d['uc_area']):+.4f}) | P = {P_layer:+.1f} uC/cm^2 over layer "
                     f"({d['t_layer']*1e9:.3f} nm) = {P_uc:+.1f} uC/cm^2 over uc ({c_uc*1e9:.3f} nm)")
        lines.append(f"            eps in layer {eps_layer:.1f} -> lattice screens {(1-1/eps_layer)*100:.2f}% of its own field | "
                     f"step from this dipole alone {step_own:+.4f} eV | band step across layer {U[k1]-U[k0]:+.4f} eV "
                     f"(dipole + background field)")
    # --- balance ---
    p_free = p_fix + p_mob
    p_scr = p_free + p_latt
    explicit = p_ins + p_free
    lines.append(f" Polarization balance (dipole moment per area, e*nm/uc of {mat_ref}; origin = {side} Dirichlet end; "
                 f"- = pointing towards z=0):")
    lines.append(f"   inserted dipoles {f(p_ins):+.5f} | free charges {f(p_free):+.5f} "
                 f"[ions/doping {f(p_fix):+.5f}, electrons/holes {f(p_mob):+.5f}] | lattice (eps-1) {f(p_latt):+.5f}")
    lines.append(f"   screening = net - inserted = {f(p_scr):+.5f} | net = {f(p_ins + p_scr):+.5f} = -eps0*dU -> "
                 f"potential step across stack dU = {U[-1]-U[0]:+.4f} eV | sum check {residual/scale:.1e}")
    if explicit != 0.0:
        lines.append(f"   -> lattice screens {-p_latt/explicit*100:.2f}% of all explicit dipoles (inserted + free charges); "
                     f"{-p_vac/explicit*100:.2f}% is left unscreened (= the band step across the stack)")
    nums = {"p_ins": p_ins, "p_fix": p_fix, "p_mob": p_mob, "p_latt": p_latt, "p_vac": p_vac,
            "residual": residual, "rel_residual": residual / scale, "A_ref": A_ref}
    return lines, nums
