# plotting.py
"""
All visualization for the heterostructure solver.

All densities plotted are 3D volume densities (m^-3) on log axes; settings in config["plot_settings"].

- add_layer_background : coloured atomic-layer bands / material interface markers
- plot_1 ... plot_7    : the individual figures of a single run
- plot_all_results     : draws plots 1-7 and shows them
- plot_thickness_sweep : critical-thickness curve of a thickness sweep
"""
import numpy as np
import matplotlib.pyplot as plt

ScF = 0.05  # Scaling factor for visualization heights


def _plot_settings(config):
    return (config or {}).get("plot_settings", {})


def add_layer_background(ax, flat_layers, show_colors=True, config=None):
    """Layer shading behind a plot.

    Colors: config["plot_settings"]["material_colors"] = {"STO": "#f4a261", ...} gives every atomic
    layer of that material the same color; materials not listed fall back to one tab20 color per
    atomic-layer name. Opacity of the shading: config["plot_settings"]["background_alpha"] (0-1).
    """
    ps = _plot_settings(config)
    mat_colors = ps.get("material_colors", {}) or {}
    alpha = float(ps.get("background_alpha", 0.25))
    if show_colors:
        # Map each unique atomic layer name to a distinct color dynamically (fallback)
        unique_layer_names = list(dict.fromkeys(L["name"] for L in flat_layers))
        cmap = plt.colormaps.get_cmap("tab20")
        color_map = {name: cmap(i % 20) for i, name in enumerate(unique_layer_names)}

    curr_z = 0.0
    prev_mat = None

    for L in flat_layers:
        t = L["t_layer"] * 1e9
        mat_name = L["mat_name"]

        if show_colors:
            layer_color = mat_colors.get(mat_name, color_map[L["name"]])
            # The blurred thickness spanning exactly the layer's width
            ax.axvspan(curr_z, curr_z + t, facecolor=layer_color, alpha=alpha, edgecolor='none', zorder=-100)
            # The central spike at the atomic center
            z_center = curr_z + (t / 2.0)
            ax.axvline(z_center, color=layer_color, alpha=min(1.0, 2.4 * alpha), lw=1.5, zorder=-99)

        # Material interface boundary marker (drawn in both modes)
        if prev_mat is not None and mat_name != prev_mat:
            ax.axvline(curr_z, color='black', lw=2, linestyle='--', zorder=-50)

        curr_z += t
        prev_mat = mat_name


def band_ylim(res, config, margin=0.3):
    """Energy window of the band plots: the bent and flat bands of the material with the largest
    band gap (plus Ef), so a vacuum / barrier layer does not stretch the axis.

    Materials listed in plot_settings["band_ylim_exclude"] are ignored; by default every material
    whose gap exceeds plot_settings["max_real_gap_eV"] (20 eV) is treated as vacuum and ignored.
    plot_settings["band_ylim"] = (ymin, ymax) overrides everything.
    """
    ps = _plot_settings(config)
    if ps.get("band_ylim"):
        return tuple(ps["band_ylim"])
    z, fl = res["z"], res["flat_layers"]
    B = np.cumsum([0.0] + [L["t_layer"] for L in fl])
    lmap = np.clip(np.searchsorted(B, z, side="right") - 1, 0, len(fl) - 1)
    mats = np.array([fl[i]["mat_name"] for i in lmap])
    exclude = set(ps.get("band_ylim_exclude", []))
    max_gap = float(ps.get("max_real_gap_eV", 20.0))
    gaps = {}
    for L in fl:
        g = L["Ec"] - L["Ev"]
        if L["mat_name"] not in exclude and g <= max_gap:
            gaps[L["mat_name"]] = max(gaps.get(L["mat_name"], 0.0), g)
    if not gaps:          # nothing left: fall back to the whole device
        sel = np.ones(len(z), bool)
    else:
        sel = mats == max(gaps, key=gaps.get)
    lo = min(np.min(res["Vtot_v"][sel]), np.min(res["V0_v"][sel]), res["Ef"])
    hi = max(np.max(res["Vtot_c"][sel]), np.max(res["V0_c"][sel]), res["Ef"])
    return lo - margin, hi + margin


def density_window(config, *arrays):
    """(floor, top) for log-scale density axes: top = 3x the largest value; floor =
    plot_settings["density_floor"] (m^-3) if given, else top * plot_settings["density_decades"]
    (default 1e-10, i.e. 10 decades), never below 1e10 m^-3."""
    ps = _plot_settings(config)
    peak = max([np.max(np.abs(a)) for a in arrays if a is not None and np.size(a)] + [1e10])
    top = 3.0 * peak
    floor = ps.get("density_floor") or max(top * float(ps.get("density_decades", 1e-10)), 1e10)
    return floor, top


def _log_curve(ax, x, y, floor, **kw):
    """Plot positive values of y on a log axis; values below the floor are hidden (no zeros)."""
    ax.plot(x, np.where(y > floor, y, np.nan), **kw)


def _fixed_charge_deltas(res):
    """Fixed + / - delta charges (doping/polar planes + dipoles) used by the charge plots."""
    fixed_pos = res["ND_plus_z"] + res["N_D_dipole"]
    fixed_neg = res["NA_minus_z"] + res["N_A_dipole"]

    pos_idx = np.where(fixed_pos > 1e10)[0]
    neg_idx = np.where(fixed_neg > 1e10)[0]
    return fixed_pos, fixed_neg, pos_idx, neg_idx


def _charge_distribution(ax, res, config, title):
    """Charge distribution with two y-axes (all volume densities, m^-3):

    left  (linear, signed): fixed charges -- polar planes, doping and dipoles (+ red / - blue) -- and
                            the net charge density rho_net, so signs and dopings stay visible
    right (log): electron density N_3D and hole density P_3D
    plot_settings["charge_linear_ylim"] = (ymin, ymax) fixes the left axis (e.g. to see a doping step
    next to much taller delta-dipole spikes); plot_settings["charge_xlim"] = (zmin, zmax) in nm zooms in.
    """
    z, N_3D, P_3D, rho_net = res["z"], res["N_3D"], res["P_3D"], res["rho_net"]
    fixed_pos, fixed_neg, pos_idx, neg_idx = _fixed_charge_deltas(res)
    x = z * 1e9

    # --- left: linear, signed ---
    if len(pos_idx) > 0:
        ax.fill_between(x, 0, fixed_pos, color='darkred', alpha=0.35, lw=0)
        ax.plot(x, fixed_pos, color='darkred', lw=1.5, label=r'Fixed $+$ ($N_D^+$ + dipole $+$)')
    if len(neg_idx) > 0:
        ax.fill_between(x, 0, -fixed_neg, color='darkblue', alpha=0.35, lw=0)
        ax.plot(x, -fixed_neg, color='darkblue', lw=1.5, label=r'Fixed $-$ ($N_A^-$ + dipole $-$)')
    ax.plot(x, rho_net, color='black', lw=1.2, alpha=0.8, label=r'$\rho_{net}$ (signed)')
    ax.axhline(0, color='gray', ls='--', lw=1)
    lim = _plot_settings(config).get("charge_linear_ylim")
    if lim:
        ax.set_ylim(*lim)
    xlim = _plot_settings(config).get("charge_xlim") or (0, x[-1])
    ax.set(xlim=tuple(xlim), xlabel='z (nm)', ylabel=r'Fixed / net charge density (m$^{-3}$, linear)', title=title)
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    ax.grid(alpha=0.3)

    # --- right: log, magnitudes ---
    axr = ax.twinx()
    floor, top = density_window(config, N_3D, P_3D)
    _log_curve(axr, x, N_3D, floor, lw=2, color='blue', label=r'$N_{3D}$ (electrons)')
    _log_curve(axr, x, P_3D, floor, lw=2, color='red', label=r'$P_{3D}$ (holes)')
    axr.set_yscale('log')
    axr.set_ylim(floor, top)
    axr.set_ylabel(r'Carrier density (m$^{-3}$, log)')

    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = axr.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, loc='upper left', bbox_to_anchor=(1.10, 1.0), fontsize=9)
    return axr


# =============================================================
# PLOTS WITH LAYER SHADING
# =============================================================

def plot_1_band_diagram_wavefunctions(res, config):
    """Band Diagram and Wavefunctions (config key: band_diagram_wavefunctions)"""
    z, Vtot_c, Vtot_v = res["z"], res["Vtot_c"], res["Vtot_v"]
    V0_c, V0_v, Ef = res["V0_c"], res["V0_v"], res["Ef"]
    Ec, Ev, psic, psiv = res["Ec"], res["Ev"], res["psic"], res["psiv"]
    flat_layers = res["flat_layers"]

    fig1, ax1 = plt.subplots(figsize=(10, 6))
    add_layer_background(ax1, flat_layers, show_colors=True, config=config)

    ax1.plot(z*1e9, V0_c, 'b--', alpha=0.5, lw=2, label='Ec (Flat)', zorder=10)
    ax1.plot(z*1e9, V0_v, 'r--', alpha=0.5, lw=2, label='Ev (Flat)', zorder=10)
    ax1.plot(z*1e9, Vtot_c, 'b', lw=2, label='Ec (Bent)', zorder=10)
    ax1.plot(z*1e9, Vtot_v, 'r', lw=2, label='Ev (Bent)', zorder=10)
    ax1.axhline(Ef, color='lime', ls='--', lw=2, label=r'$E_f$', zorder=10)

    if len(Ec) > 0:
        for i in range(min(config["n"], len(Ec))):
            intensity = np.abs(psic[:, i])**2
            if np.max(intensity) > 1e-15:
                norm_intensity = intensity / np.max(intensity)
                ax1.plot(z*1e9, np.full(len(z), Ec[i]) + (ScF * norm_intensity),
                         color='blue', linestyle=':', lw=1.5, alpha=0.8, zorder=5)
                ax1.axhline(Ec[i], color='blue', linestyle='-', lw=0.5, alpha=0.3, zorder=4)

    if len(Ev) > 0:
        for i in range(min(config["n"], len(Ev))):
            intensity = np.abs(psiv[:, i])**2
            if np.max(intensity) > 1e-15:
                norm_intensity = intensity / np.max(intensity)
                ax1.plot(z*1e9, np.full(len(z), Ev[i]) - (ScF * norm_intensity),
                         color='red', linestyle=':', lw=1.5, alpha=0.8, zorder=5)
                ax1.axhline(Ev[i], color='red', linestyle='-', lw=0.5, alpha=0.3, zorder=4)

    ax1.set(xlim=(0, z[-1]*1e9), ylim=band_ylim(res, config),
            xlabel='z (nm)', ylabel='Energy (eV)', title='Band Diagram and Wavefunctions')
    ax1.grid(alpha=0.3)
    ax1.legend(loc='center right')
    fig1.tight_layout()
    return fig1


def plot_2_density_heatmap(res, config):
    """Carrier Distribution Heatmap (config key: carrier_distribution_heatmap)"""
    z, Vtot_c, Vtot_v, Ef = res["z"], res["Vtot_c"], res["Vtot_v"], res["Ef"]
    N_3D, P_3D = res["N_3D"], res["P_3D"]
    flat_layers = res["flat_layers"]

    fig2, ax2 = plt.subplots(figsize=(10, 6))
    add_layer_background(ax2, flat_layers, show_colors=True, config=config)

    ax2.plot(z*1e9, Vtot_c, 'b', lw=2, label='Ec (Bent)', zorder=10)
    ax2.plot(z*1e9, Vtot_v, 'r', lw=2, label='Ev (Bent)', zorder=10)
    ax2.axhline(Ef, color='lime', ls='--', lw=2, label=r'$E_f$', zorder=10)

    if np.max(N_3D) > 1e10:
        norm_N3D = N_3D / np.max(N_3D)
        mask_e = norm_N3D > 1e-3
        ax2.scatter(z[mask_e]*1e9, Vtot_c[mask_e] + (ScF * norm_N3D[mask_e]),
                    c=norm_N3D[mask_e], cmap='Blues', s=25, alpha=0.9, edgecolor='none', zorder=5, label='Electron Heatmap')

    if np.max(P_3D) > 1e10:
        norm_P3D = P_3D / np.max(P_3D)
        mask_h = norm_P3D > 1e-3
        ax2.scatter(z[mask_h]*1e9, Vtot_v[mask_h] - (ScF * norm_P3D[mask_h]),
                    c=norm_P3D[mask_h], cmap='Reds', s=25, alpha=0.9, edgecolor='none', zorder=5, label='Hole Heatmap')

    ax2.set(xlim=(0, z[-1]*1e9), ylim=band_ylim(res, config),
            xlabel='z (nm)', ylabel='Energy (eV)', title='Carrier Distribution Heatmap')
    ax2.grid(alpha=0.3)
    handles, labels = ax2.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax2.legend(by_label.values(), by_label.keys(), loc='center right')
    fig2.tight_layout()
    return fig2


def plot_3_charge_distribution(res, config):
    """Charge Distribution with Layer Shading (config key: charge_distribution_shaded)"""
    z, N_3D, P_3D, rho_net = res["z"], res["N_3D"], res["P_3D"], res["rho_net"]
    flat_layers = res["flat_layers"]

    fig3, ax3 = plt.subplots(figsize=(12, 5))
    add_layer_background(ax3, flat_layers, show_colors=True, config=config)

    _charge_distribution(ax3, res, config, 'Charge Distribution with Layer Shading')
    fig3.tight_layout()
    return fig3


# =============================================================
# PLOTS WITH PLAIN BACKGROUND
# =============================================================

def plot_4_conduction_band_wavefunctions(res, config):
    """Conduction Band and Electron Wavefunctions (config key: conduction_band_wavefunctions)"""
    z, Vtot_c, V0_c = res["z"], res["Vtot_c"], res["V0_c"]
    Ec, psic = res["Ec"], res["psic"]
    flat_layers = res["flat_layers"]

    fig4, ax4 = plt.subplots(figsize=(10, 6))
    add_layer_background(ax4, flat_layers, show_colors=False, config=config)

    # Conduction band
    ax4.plot(
        z*1e9, V0_c,
        'b--', alpha=0.5, lw=2,
        label='Ec (Flat)', zorder=10
    )

    ax4.plot(
        z*1e9, Vtot_c,
        'b', lw=2,
        label='Ec (Bent)', zorder=10
    )

    # Electron wavefunctions
    if len(Ec) > 0:
        for i in range(min(config["n"], len(Ec))):
            intensity = np.abs(psic[:, i])**2

            if np.max(intensity) > 1e-15:
                norm_intensity = intensity / np.max(intensity)

                ax4.plot(
                    z*1e9,
                    np.full(len(z), Ec[i]) + (ScF * norm_intensity),
                    color='blue',
                    linestyle='-',
                    lw=1.5,
                    alpha=0.8,
                    zorder=5
                )

                ax4.axhline(
                    Ec[i],
                    color='blue',
                    linestyle='-',
                    lw=0.5,
                    alpha=0.3,
                    zorder=4
                )

    # Fixed viewing window
    ax4.set(
        xlim=(20, 40),
        ylim=(-4, 0),
        xlabel='z (nm)',
        ylabel='Energy (eV)',
        title='Conduction Band and Electron Wavefunctions'
    )

    ax4.grid(alpha=0.3)
    ax4.legend(loc='center right')
    fig4.tight_layout()
    return fig4


def plot_5_band_diagram_carrier_density(res, config):
    """Band Diagram and Carrier Density (config key: band_diagram_carrier_density)"""
    z, Vtot_c, Vtot_v, Ef = res["z"], res["Vtot_c"], res["Vtot_v"], res["Ef"]
    N_3D, P_3D = res["N_3D"], res["P_3D"]
    flat_layers = res["flat_layers"]

    fig5, ax5 = plt.subplots(figsize=(10, 6))
    add_layer_background(ax5, flat_layers, show_colors=False, config=config)

    ax5.plot(z*1e9, Vtot_c, 'b', lw=2, label='Ec (Bent)', zorder=10)
    ax5.plot(z*1e9, Vtot_v, 'r', lw=2, label='Ev (Bent)', zorder=10)
    ax5.axhline(Ef, color='lime', ls='--', lw=2, label=r'$E_f$', zorder=10)

    ax5.set(xlim=(0, z[-1]*1e9), ylim=band_ylim(res, config),
            xlabel='z (nm)', ylabel='Energy (eV)', title='Band Diagram and Carrier Density')
    ax5.grid(alpha=0.3)

    # Secondary Y-axis for 1D Density Plot
    ax5_twin = ax5.twinx()

    floor, top = density_window(config, N_3D, P_3D)
    if np.max(N_3D) > floor:
        _log_curve(ax5_twin, z*1e9, N_3D, floor, color='darkblue', lw=2, label=r'Electron Density ($N_{3D}$)', zorder=15)
    if np.max(P_3D) > floor:
        _log_curve(ax5_twin, z*1e9, P_3D, floor, color='darkred', lw=2, label=r'Hole Density ($P_{3D}$)', zorder=15)

    ax5_twin.set_ylabel(r'Carrier Density (m$^{-3}$, 3D)')
    ax5_twin.set_yscale('log')
    ax5_twin.set_ylim(floor, top)

    handles1, labels1 = ax5.get_legend_handles_labels()
    handles2, labels2 = ax5_twin.get_legend_handles_labels()
    by_label5 = dict(zip(labels1 + labels2, handles1 + handles2))
    ax5.legend(by_label5.values(), by_label5.keys(), loc='center right')
    fig5.tight_layout()
    return fig5


def plot_6_charge_distribution_clean(res, config):
    """Charge Distribution (config key: charge_distribution)"""
    z, N_3D, P_3D, rho_net = res["z"], res["N_3D"], res["P_3D"], res["rho_net"]
    flat_layers = res["flat_layers"]

    fig6, ax6 = plt.subplots(figsize=(12, 5))
    add_layer_background(ax6, flat_layers, show_colors=False, config=config)

    _charge_distribution(ax6, res, config, 'Charge Distribution')
    fig6.tight_layout()
    return fig6


def plot_7_field_and_dielectric(res, config):
    """Electric Field and Dielectric Constant (config key: electric_field_dielectric_constant)"""
    z = res["z"]
    flat_layers = res["flat_layers"]
    Epsi_final = res["Epsi_arr"]
    # Convert E_field from V/m to MV/cm for easier reading on the plot
    E_final_MV_cm = res["E_field"] * 1e-8

    fig7, ax7 = plt.subplots(figsize=(10, 6))
    add_layer_background(ax7, flat_layers, show_colors=False, config=config)

    # Plot Electric Field on the primary y-axis
    ax7.plot(z*1e9, E_final_MV_cm, color='purple', lw=2, label='Electric Field (MV/cm)', zorder=10)
    ax7.set(xlabel='z (nm)', ylabel='Electric Field (MV/cm)', title='Electric Field and Dielectric Constant')
    ax7.grid(alpha=0.3)

    # Plot Dielectric Constant on a secondary y-axis
    ax7_twin = ax7.twinx()
    ax7_twin.plot(z*1e9, Epsi_final, color='orange', lw=2, linestyle='--', label=r'Dielectric Constant ($\epsilon_r$)', zorder=15)
    ax7_twin.set_ylabel(r'Relative Permittivity ($\epsilon_r$)')

    # Combine legends
    handles1, labels1 = ax7.get_legend_handles_labels()
    handles2, labels2 = ax7_twin.get_legend_handles_labels()
    by_label7 = dict(zip(labels1 + labels2, handles1 + handles2))
    ax7.legend(by_label7.values(), by_label7.keys(), loc='upper right')

    fig7.tight_layout()
    return fig7


# Plot keys used in config["plots"] (True = draw, False = skip); each key matches the figure title
PLOTS = {
    "band_diagram_wavefunctions":         plot_1_band_diagram_wavefunctions,     # Band Diagram and Wavefunctions
    "carrier_distribution_heatmap":       plot_2_density_heatmap,                # Carrier Distribution Heatmap
    "charge_distribution_shaded":         plot_3_charge_distribution,            # Charge Distribution with Layer Shading
    "conduction_band_wavefunctions":      plot_4_conduction_band_wavefunctions,  # Conduction Band and Electron Wavefunctions
    "band_diagram_carrier_density":       plot_5_band_diagram_carrier_density,   # Band Diagram and Carrier Density
    "charge_distribution":                plot_6_charge_distribution_clean,      # Charge Distribution
    "electric_field_dielectric_constant": plot_7_field_and_dielectric,           # Electric Field and Dielectric Constant
}
DEFAULT_PLOTS = {k: True for k in PLOTS}
DEFAULT_PLOTS["band_diagram_wavefunctions"] = False

# Earlier key names, still accepted
PLOT_KEY_ALIASES = {
    "carrier_heatmap": "carrier_distribution_heatmap",
    "cb_wavefunctions": "conduction_band_wavefunctions",
    "band_carrier_density": "band_diagram_carrier_density",
    "charge_distribution_clean": "charge_distribution",
    "field_dielectric": "electric_field_dielectric_constant",
}


def add_minor_ticks(fig):
    """Minor ticks on every axis (incl. twin axes): linear axes get evenly spaced minor ticks,
    log axes get the 2..9 sub-decade ticks, so linear and log axes are easy to tell apart."""
    from matplotlib.ticker import LogLocator, NullFormatter
    for ax in fig.axes:
        ax.minorticks_on()
        for axis, scale in ((ax.yaxis, ax.get_yscale()), (ax.xaxis, ax.get_xscale())):
            if scale == "log":       # keep 2..9 minor ticks even when the axis spans many decades
                axis.set_major_locator(LogLocator(base=10.0, numticks=30))
                axis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * 1.0, numticks=300))
                axis.set_minor_formatter(NullFormatter())
        ax.tick_params(which="major", length=6, width=1.0)
        ax.tick_params(which="minor", length=3, width=0.8)


def apply_uniform_linewidth(fig, lw=2.0):
    """Give every line in the figure (data curves, band edges, wavefunctions, Ef, subband levels,
    zero lines, interface markers) and every legend entry the same linewidth.
    The thin per-layer centre lines of the layer shading (zorder -99) are background and are kept."""
    for ax in fig.axes:
        for line in ax.get_lines():
            if line.get_zorder() > -99:
                line.set_linewidth(lw)
        leg = ax.get_legend()
        if leg is not None:
            handles = leg.legend_handles if hasattr(leg, "legend_handles") else leg.legendHandles
            for h in handles:
                if type(h).__name__ == "Line2D":
                    h.set_linewidth(lw)


def plot_all_results(res, config, show=True):
    """Draw the plots switched on in config["plots"] (missing keys use DEFAULT_PLOTS) and show them."""
    user = config.get("plots") or {}
    if "charge_distribution_clean" in user and "charge_distribution_shaded" not in user:
        # old-style config: there "charge_distribution" meant the shaded version
        user = dict(user)
        if "charge_distribution" in user:
            user["charge_distribution_shaded"] = user.pop("charge_distribution")
    selection = dict(DEFAULT_PLOTS)
    for k, v in user.items():
        selection[PLOT_KEY_ALIASES.get(k, k)] = v
    unknown = set(selection) - set(PLOTS)
    if unknown:
        print(f"[plots] unknown plot key(s) in config['plots']: {sorted(unknown)}; valid: {list(PLOTS)}")
    figs = [func(res, config) for key, func in PLOTS.items() if selection.get(key, False)]
    lw = float(_plot_settings(config).get("linewidth", 2.0))
    for fig in figs:
        add_minor_ticks(fig)
        apply_uniform_linewidth(fig, lw)
    if show:
        plt.show()
    return figs


def plot_thickness_sweep(uc_array, e_per_uc_array, fixed_stack, sweep_mat, show=True):
    """Critical Thickness Curve (electrons per unit cell vs. swept thickness)."""
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(uc_array, e_per_uc_array, marker='o', lw=2, color='blue', markersize=8, label='Simulated')

    # Add a horizontal dashed line to show the 0.5 e/uc theoretical limit
    ax.axhline(0.5, color='red', linestyle='--', label='Polar Catastrophe Limit (0.5 $e^-$/u.c.)')

    ax.set_xlabel(f'{sweep_mat} Thickness (Unit Cells)', fontsize=12)
    ax.set_ylabel(r'Electron Density ($e^- / u.c.$)', fontsize=12)
    ax.set_title(f'Critical Thickness Sweep: {fixed_stack} / {sweep_mat}', fontsize=12)
    add_minor_ticks(fig)
    ax.grid(True, alpha=0.4)
    ax.legend()

    apply_uniform_linewidth(fig)
    # e/uc is best viewed on a linear scale to see it asymptotically approach 0.5
    fig.tight_layout()
    if show:
        plt.show()
    return fig
