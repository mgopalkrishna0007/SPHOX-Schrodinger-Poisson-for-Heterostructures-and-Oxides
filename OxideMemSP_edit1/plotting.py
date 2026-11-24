# plotting.py
"""
All visualization for the heterostructure solver (plot logic unchanged from master_run.py).

- add_layer_background : coloured atomic-layer bands / material interface markers
- plot_1 ... plot_7    : the individual figures of a single run
- plot_all_results     : draws plots 1-7 and shows them
- plot_thickness_sweep : critical-thickness curve of a thickness sweep
"""
import numpy as np
import matplotlib.pyplot as plt

ScF = 0.05  # Scaling factor for visualization heights


def add_layer_background(ax, flat_layers, show_colors=True):
    if show_colors:
        # Map each unique atomic layer name to a distinct color dynamically
        unique_layer_names = list(dict.fromkeys(L["name"] for L in flat_layers))
        cmap = plt.colormaps.get_cmap("tab20")
        color_map = {name: cmap(i % 20) for i, name in enumerate(unique_layer_names)}

    curr_z = 0.0
    prev_mat = None

    for L in flat_layers:
        t = L["t_layer"] * 1e9
        mat_name = L["mat_name"]

        if show_colors:
            layer_color = color_map[L["name"]]
            # The blurred thickness spanning exactly the layer's width
            ax.axvspan(curr_z, curr_z + t, facecolor=layer_color, alpha=0.25, edgecolor='none', zorder=-100)
            # The central spike at the atomic center
            z_center = curr_z + (t / 2.0)
            ax.axvline(z_center, color=layer_color, alpha=0.6, lw=1.5, zorder=-99)

        # Material interface boundary marker (drawn in both modes)
        if prev_mat is not None and mat_name != prev_mat:
            ax.axvline(curr_z, color='black', lw=2, linestyle='--', zorder=-50)

        curr_z += t
        prev_mat = mat_name


def _fixed_charge_deltas(res):
    """Fixed + / - delta charges (doping/polar planes + dipoles) used by the charge plots."""
    fixed_pos = res["ND_plus_z"] + res["N_D_dipole"]
    fixed_neg = res["NA_minus_z"] + res["N_A_dipole"]

    pos_idx = np.where(fixed_pos > 1e10)[0]
    neg_idx = np.where(fixed_neg > 1e10)[0]
    return fixed_pos, fixed_neg, pos_idx, neg_idx


# =============================================================
# SET 1: ORIGINAL PLOTS (WITH COLORS & HEATMAPS)
# =============================================================

def plot_1_band_diagram_wavefunctions(res, config):
    """--- PLOT 1: Original Band Diagram + Wavefunctions ---"""
    z, Vtot_c, Vtot_v = res["z"], res["Vtot_c"], res["Vtot_v"]
    V0_c, V0_v, Ef = res["V0_c"], res["V0_v"], res["Ef"]
    Ec, Ev, psic, psiv = res["Ec"], res["Ev"], res["psic"], res["psiv"]
    flat_layers = res["flat_layers"]

    fig1, ax1 = plt.subplots(figsize=(10, 6))
    add_layer_background(ax1, flat_layers, show_colors=True)

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

    ax1.set(xlim=(0, z[-1]*1e9), ylim=(np.min(Vtot_v) - 0.2, np.max(Vtot_c) + 0.2),
            xlabel='z (nm)', ylabel='Energy (eV)', title='[Original] Band Diagram & Individual Wavefunctions')
    ax1.grid(alpha=0.3)
    ax1.legend(loc='center right')
    fig1.tight_layout()
    return fig1


def plot_2_density_heatmap(res, config):
    """--- PLOT 2: Original Band Diagram + Total Density Heatmap ---"""
    z, Vtot_c, Vtot_v, Ef = res["z"], res["Vtot_c"], res["Vtot_v"], res["Ef"]
    N_3D, P_3D = res["N_3D"], res["P_3D"]
    flat_layers = res["flat_layers"]

    fig2, ax2 = plt.subplots(figsize=(10, 6))
    add_layer_background(ax2, flat_layers, show_colors=True)

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

    ax2.set(xlim=(0, z[-1]*1e9), ylim=(np.min(Vtot_v) - 0.2, np.max(Vtot_c) + 0.2),
            xlabel='z (nm)', ylabel='Energy (eV)', title='[Original] Carrier Spatial Distribution Heatmap')
    ax2.grid(alpha=0.3)
    handles, labels = ax2.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax2.legend(by_label.values(), by_label.keys(), loc='center right')
    fig2.tight_layout()
    return fig2


def plot_3_charge_distribution(res, config):
    """--- PLOT 3: Original Charge Distribution ---"""
    z, N_3D, P_3D, rho_net = res["z"], res["N_3D"], res["P_3D"], res["rho_net"]
    flat_layers = res["flat_layers"]

    fig3, ax3 = plt.subplots(figsize=(10, 5))
    add_layer_background(ax3, flat_layers, show_colors=True)

    ax3.plot(z*1e9, -N_3D, lw=2, color='blue', label=r'$-N_{3D}$ (Electrons)')
    ax3.plot(z*1e9, P_3D, lw=2, color='red', label=r'$P_{3D}$ (Holes)')
    ax3.plot(z*1e9, rho_net, lw=2, color='black', alpha=0.8, label=r'$\rho_{net}$')

    fixed_pos, fixed_neg, pos_idx, neg_idx = _fixed_charge_deltas(res)

    if len(pos_idx) > 0:
        ax3.vlines(z[pos_idx]*1e9, 0, fixed_pos[pos_idx], color='darkred', linestyle='-', lw=2, label=r'Fixed $+$ (Delta)')
    if len(neg_idx) > 0:
        ax3.vlines(z[neg_idx]*1e9, 0, -fixed_neg[neg_idx], color='darkblue', linestyle='-', lw=2, label=r'Fixed $-$ (Delta)')

    ax3.set_yscale('symlog', linthresh=1e25)
    ax3.axhline(0, color='gray', ls='--', lw=1)
    ax3.set(xlabel='z (nm)', ylabel=r'Charge Density (m$^{-3}$)', title='[Original] Charge Distribution Profile')
    ax3.grid(alpha=0.3)
    ax3.legend(loc='upper right', bbox_to_anchor=(1.25, 1.0))
    fig3.tight_layout()
    return fig3


# =============================================================
# SET 2: CLEAN PLOTS (NO COLORS, 1D DENSITY INSTEAD OF HEATMAP)
# =============================================================

def plot_4_conduction_band_wavefunctions(res, config):
    """--- PLOT 4: Clean Conduction Band + Electron Wavefunctions ---"""
    z, Vtot_c, V0_c = res["z"], res["Vtot_c"], res["V0_c"]
    Ec, psic = res["Ec"], res["psic"]
    flat_layers = res["flat_layers"]

    fig4, ax4 = plt.subplots(figsize=(10, 6))
    add_layer_background(ax4, flat_layers, show_colors=False)

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
    """--- PLOT 5: Clean Band Diagram + 1D Carrier Density (Replaces Heatmap) ---"""
    z, Vtot_c, Vtot_v, Ef = res["z"], res["Vtot_c"], res["Vtot_v"], res["Ef"]
    N_3D, P_3D = res["N_3D"], res["P_3D"]
    flat_layers = res["flat_layers"]

    fig5, ax5 = plt.subplots(figsize=(10, 6))
    add_layer_background(ax5, flat_layers, show_colors=False)

    ax5.plot(z*1e9, Vtot_c, 'b', lw=2, label='Ec (Bent)', zorder=10)
    ax5.plot(z*1e9, Vtot_v, 'r', lw=2, label='Ev (Bent)', zorder=10)
    ax5.axhline(Ef, color='lime', ls='--', lw=2, label=r'$E_f$', zorder=10)

    ax5.set(xlim=(0, z[-1]*1e9), ylim=(np.min(Vtot_v) - 0.2, np.max(Vtot_c) + 0.2),
            xlabel='z (nm)', ylabel='Energy (eV)', title='[Clean] Band Diagram & 1D Carrier Density')
    ax5.grid(alpha=0.3)

    # Secondary Y-axis for 1D Density Plot
    ax5_twin = ax5.twinx()

    if np.max(N_3D) > 1e10:
        ax5_twin.plot(z*1e9, N_3D, color='darkblue', lw=2, label=r'Electron Density ($N_{3D}$)', zorder=15)
    if np.max(P_3D) > 1e10:
        ax5_twin.plot(z*1e9, P_3D, color='darkred', lw=2, label=r'Hole Density ($P_{3D}$)', zorder=15)

    ax5_twin.set_ylabel(r'Carrier Density (m$^{-3}$)')
    ax5_twin.set_yscale('symlog', linthresh=1e18)

    handles1, labels1 = ax5.get_legend_handles_labels()
    handles2, labels2 = ax5_twin.get_legend_handles_labels()
    by_label5 = dict(zip(labels1 + labels2, handles1 + handles2))
    ax5.legend(by_label5.values(), by_label5.keys(), loc='center right')
    fig5.tight_layout()
    return fig5


def plot_6_charge_distribution_clean(res, config):
    """--- PLOT 6: Clean Charge Distribution ---"""
    z, N_3D, P_3D, rho_net = res["z"], res["N_3D"], res["P_3D"], res["rho_net"]
    flat_layers = res["flat_layers"]

    fig6, ax6 = plt.subplots(figsize=(10, 5))
    add_layer_background(ax6, flat_layers, show_colors=False)

    ax6.plot(z*1e9, -N_3D, lw=2, color='blue', label=r'$-N_{3D}$ (Electrons)')
    ax6.plot(z*1e9, P_3D, lw=2, color='red', label=r'$P_{3D}$ (Holes)')
    ax6.plot(z*1e9, rho_net, lw=2, color='black', alpha=0.8, label=r'$\rho_{net}$')

    fixed_pos, fixed_neg, pos_idx, neg_idx = _fixed_charge_deltas(res)

    if len(pos_idx) > 0:
        ax6.vlines(z[pos_idx]*1e9, 0, fixed_pos[pos_idx], color='darkred', linestyle='-', lw=2, label=r'Fixed $+$ (Delta)')
    if len(neg_idx) > 0:
        ax6.vlines(z[neg_idx]*1e9, 0, -fixed_neg[neg_idx], color='darkblue', linestyle='-', lw=2, label=r'Fixed $-$ (Delta)')

    ax6.set_yscale('symlog', linthresh=1e25)
    ax6.axhline(0, color='gray', ls='--', lw=1)
    ax6.set(xlabel='z (nm)', ylabel=r'Charge Density (m$^{-3}$)', title='[Clean] Charge Distribution Profile')
    ax6.grid(alpha=0.3)
    ax6.legend(loc='upper right', bbox_to_anchor=(1.25, 1.0))
    fig6.tight_layout()
    return fig6


def plot_7_field_and_dielectric(res, config):
    """--- PLOT 7: Electric Field and Dynamic Dielectric Constant ---"""
    z = res["z"]
    flat_layers = res["flat_layers"]
    Epsi_final = res["Epsi_arr"]
    # Convert E_field from V/m to MV/cm for easier reading on the plot
    E_final_MV_cm = res["E_field"] * 1e-8

    fig7, ax7 = plt.subplots(figsize=(10, 6))
    add_layer_background(ax7, flat_layers, show_colors=False)

    # Plot Electric Field on the primary y-axis
    ax7.plot(z*1e9, E_final_MV_cm, color='purple', lw=2, label='Electric Field (MV/cm)', zorder=10)
    ax7.set(xlabel='z (nm)', ylabel='Electric Field (MV/cm)', title='Electric Field & Dielectric Profile vs Space')
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


def plot_all_results(res, config, show=True):
    """Draw plots 1-7 of a single run (and show them)."""
    figs = [
        plot_1_band_diagram_wavefunctions(res, config),
        plot_2_density_heatmap(res, config),
        plot_3_charge_distribution(res, config),
        plot_4_conduction_band_wavefunctions(res, config),
        plot_5_band_diagram_carrier_density(res, config),
        plot_6_charge_distribution_clean(res, config),
        plot_7_field_and_dielectric(res, config),
    ]
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
    ax.grid(True, alpha=0.4)
    ax.legend()

    # e/uc is best viewed on a linear scale to see it asymptotically approach 0.5
    fig.tight_layout()
    if show:
        plt.show()
    return fig
