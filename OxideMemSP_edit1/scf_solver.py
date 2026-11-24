# scf_solver.py
"""
FEniCS initialization and the self-consistent Poisson-Schrödinger loop.

run_self_consistent_loop(config, device) takes the device built by
simulation.build_device(...) and returns the results dictionary used for plotting.
"""
import numpy as np
import dolfin as df

from master_helper import find_Ef_general, Schroed1D_dolfin, Poisson_dolfin
from material_properties import compute_permittivity
from mixers import create_mixing_scheme
from stack_geometry import E_CHARGE


def create_function_space(z):
    """--- 7. FEniCS Initialization ---"""
    Nz = len(z)
    mesh = df.IntervalMesh(Nz - 1, z[0], z[-1])
    mesh.coordinates()[:, 0] = z
    V_space = df.FunctionSpace(mesh, 'CG', 1)
    return V_space


def print_loop_report(nloop, Err, Ef, config, Q_e, Q_h, Q_ND, Q_NA, Q_net, net_integral_charge,
                      eps_old_int, eps_new_int, max_eps_change, rho_net, z, interface_idx):
    print(f"==========================================================================================")
    print(f" LOOP {nloop:03d} | Error: {Err:.3e} (Tol: {config['tol']:.1e}) | Ef: {Ef:.5f} eV")
    print(f"------------------------------------------------------------------------------------------")
    print(f" 2D Electron Dens (np.sum N_n) : {-Q_e:.4e} m^-2")
    print(f" 2D Hole Dens (np.sum P_m) : {+Q_h:.4e} m^-2")
    print(f" 2D Donor Dens (Trapz ND+) : {+Q_ND:.4e} m^-2")
    print(f" 2D Acceptor Dens (Trapz NA-) : {-Q_NA:.4e} m^-2")
    print(f" Net Charge to Poisson (Trapz) : {Q_net:.4e} m^-2")
    print(f" Q Residual (Convergence check): {net_integral_charge:.4e} m^-2")
    print(f" Interface Dielectric (Old -> New): {eps_old_int:.2f} -> {eps_new_int:.2f}")
    print(f" Max Dielectric Change across mesh: {max_eps_change:.4f}")
    print(f"------------------------------------------------------------------------------------------")
    print(f" Charge Dens @ Interface Node : {rho_net[interface_idx]:.4e} m^-3 (z = {z[interface_idx]*1e9:.2f} nm)")
    print(f"==========================================================================================\n")


def run_self_consistent_loop(config, device):
    # --- Unpack the device (mesh, properties, fixed charges) ---
    z = device["z"]
    flat_layers = device["flat_layers"]
    layer_idx_map = device["layer_idx_map"]
    interface_idx = device["interface_idx"]
    eps_frozen = device.get("eps_frozen")
    V0_c, V0_v = device["V0_c"], device["V0_v"]
    meff_e_arr, meff_h_arr = device["meff_e_arr"], device["meff_h_arr"]
    Epsi_arr = device["Epsi_arr"].copy()
    N_D_array, N_A_array = device["N_D_array"], device["N_A_array"]
    N_D_dipole, N_A_dipole = device["N_D_dipole"], device["N_A_dipole"]
    P_layers = device["P_layers"]

    # --- 7. FEniCS Initialization ---
    V_space = create_function_space(z)

    scheme = create_mixing_scheme(config)
    # How the field-dependent eps is updated:
    #   "damped": eps from the previous OUTPUT field, with the mixer's own eps damping (original behaviour)
    #   "input" : eps recomputed from the INPUT potential's field every loop (after eps_start), so each
    #             loop is a pure function of Vin and the mixer only has to converge V
    eps_update = str(config.get("eps_update", "damped")).lower()
    if eps_update not in ("damped", "input"):
        raise ValueError(f"Unknown eps_update '{eps_update}'. Choose 'damped' or 'input'.")
    eps_start = config.get("eps_start", 0)
    Vin = np.zeros_like(z)
    convergence_met = False

    print(f"\nStarting Self-Consistent Poisson-Schrödinger Solver... (mixing: {scheme.name}, eps update: {eps_update})")

    try:
        # --- 8. Self-Consistent Loop ---
        for nloop in range(1, config["Nloops"] + 1):
            Vtot_c = V0_c + Vin
            Vtot_v = V0_v + Vin

            if eps_update == "input" and nloop > eps_start:
                Epsi_arr = compute_permittivity(flat_layers, layer_idx_map, np.abs(np.gradient(Vin, z)), eps_frozen)

            Ec, psic = Schroed1D_dolfin(V_space, Vtot_c, meff_e_arr, config["n"])
            Ev_inv, psiv = Schroed1D_dolfin(V_space, -Vtot_v, meff_h_arr, config["n"])
            Ev = -np.array(Ev_inv) if len(Ev_inv) > 0 else np.array([])

            if len(Ec) == 0:
                break

            Ef, N_n, P_m, ND_plus_z, NA_minus_z, net_integral_charge = find_Ef_general(
                Ec, Ev, psic, psiv, N_D_array, N_A_array, Vtot_c, Vtot_v, z, config["T"],
                meff_e_arr, meff_h_arr, config["E_bind_D"], config["E_bind_A"],
                config["fully_ionized_flag"]
            )

            N_3D = np.sum(np.tile(N_n, (len(z), 1)) * np.abs(psic)**2, axis=1) if len(Ec) > 0 else np.zeros_like(z)
            P_3D = np.sum(np.tile(P_m, (len(z), 1)) * np.abs(psiv)**2, axis=1) if len(Ev) > 0 else np.zeros_like(z)

            rho_net = P_3D - N_3D + ND_plus_z - NA_minus_z + N_D_dipole - N_A_dipole

            Vout = Poisson_dolfin(
                V_space, rho_net, Epsi_arr,
                config.get("bc_left_type", "neumann"), config.get("V_left_applied", 0.0),
                config.get("bc_right_type", "neumann"), config.get("V_right_applied", 0.0)
            )

            # --- Field-dependent dielectric constant ---
            E_Vm = np.abs(np.gradient(Vout, z))  # Absolute electric field in V/m
            Epsi_new = compute_permittivity(flat_layers, layer_idx_map, E_Vm, eps_frozen)

            if eps_update == "input":
                # eps already follows Vin; report how far eps(E_out) still is from it
                eps_err = np.max(np.abs(Epsi_new - Epsi_arr) / Epsi_arr) if nloop > eps_start else 0.0
            else:
                Epsi_arr, eps_err = scheme.update_epsilon(nloop, Epsi_new, Epsi_arr)

            Err = np.max(np.abs(Vout - Vin))

            # Output Formatting
            Q_e = np.sum(N_n)
            Q_h = np.sum(P_m)
            Q_ND = np.trapz(ND_plus_z, z)
            Q_NA = np.trapz(NA_minus_z, z)
            Q_net = np.trapz(rho_net, z)

            L_z = z[-1] - z[0]
            z0 = 0.5 * (z[0] + z[-1])
            P_elec = -E_CHARGE * np.trapz((z - z0) * N_3D, x=z) / L_z
            P_hole = E_CHARGE * np.trapz((z - z0) * P_3D, x=z) / L_z
            P_dipole_inserted = sum(P_layers)

            eps_old_int = Epsi_arr[interface_idx]
            eps_new_int = Epsi_new[interface_idx]
            max_eps_change = np.max(np.abs(Epsi_new - Epsi_arr))

            print_loop_report(nloop, Err, Ef, config, Q_e, Q_h, Q_ND, Q_NA, Q_net, net_integral_charge,
                              eps_old_int, eps_new_int, max_eps_change, rho_net, z, interface_idx)

            if scheme.is_converged(Err, eps_err):
                convergence_met = True
                break

            Vin = scheme.next_potential(Vin, Vout, Err)

    except KeyboardInterrupt:
        print("\n\n[!] SIMULATION INTERRUPTED BY USER (Ctrl+C)")
        print(f"[!] Halting at Loop {nloop}. Proceeding to plot the partially converged state...\n")

    return {
            "z": z, "Vtot_c": V0_c + Vout, "Vtot_v": V0_v + Vout,
            "V0_c": V0_c, "V0_v": V0_v, "Ef": Ef,
            "N_3D": N_3D, "P_3D": P_3D, "Ec": Ec, "Ev": Ev,
            "psic": psic, "psiv": psiv,
            "ND_plus_z": ND_plus_z, "NA_minus_z": NA_minus_z,
            "N_D_dipole": N_D_dipole, "N_A_dipole": N_A_dipole,
            "rho_net": rho_net, "flat_layers": flat_layers,
            "Epsi_arr": Epsi_arr,
            "E_field": E_Vm,
            "converged": convergence_met,
        }
