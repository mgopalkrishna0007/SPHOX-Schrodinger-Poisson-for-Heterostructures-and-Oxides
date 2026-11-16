# import numpy as np
# import dolfin as df
# import warnings
# import matplotlib.pyplot as plt

# # Import from your local files
# from master_helper import find_Ef_general, Schroed1D_dolfin, Poisson_dolfin

# # from config_BBOiLNB import config
# # from config_SIrOiLNB import config
# from config_LAOiSTO import config
# # from config_STO100iSTO111 import config
# # from config_LLOiSBO import config
# # from config_STOiSAP import config
# # from config_STOiGTO import config

# warnings.filterwarnings("ignore")
# df.set_log_level(df.LogLevel.WARNING)


# def run_heterostructure_simulation(config):
#     e = 1.602176487E-19  # e = 1.602176487×10^−19 C

#     # --- Helper Function Moved Here ---
#     def get_layer_indices_by_offset_range(flat_layers, interface_pair, direction, offset_range):
#         mat_left, mat_right = interface_pair
#         idx = -1
#         for i in range(len(flat_layers) - 1):
#             if flat_layers[i]["mat_name"] == mat_left and flat_layers[i+1]["mat_name"] == mat_right:
#                 idx = i
#                 break
#         if idx == -1: return []

#         start_off, end_off = offset_range
#         indices = []
#         if direction == 'right':
#             max_off = len(flat_layers) - 1 - idx
#             if end_off == -1 or end_off > max_off: end_off = max_off
#             for off in range(start_off, end_off + 1): indices.append(idx + off)
#         elif direction == 'left':
#             max_off = idx + 1
#             if end_off == -1 or end_off > max_off: end_off = max_off
#             for off in range(start_off, end_off + 1): indices.append(idx - off + 1)
#         return indices

#     # --- 1. Expand Stack Definition ---
#     flat_layers = []  # will contain every atomic layer like ababababcdcdcdcdcd
#     MATERIALS = config["MATERIALS"]
    
#     for item in config["stack_definition"]:
#         mat_name, count = item.split('*')
#         count = int(count)
#         mat = MATERIALS[mat_name]
#         for _ in range(count):
#             for layer in mat["layers"]:  # adds the below extra information to all the layers
#                 layer_copy = layer.copy()
#                 layer_copy["mat_name"] = mat_name
#                 layer_copy["uc_area"] = mat["uc_area"]
#                 flat_layers.append(layer_copy)

#     # --- Pre-calculate Layer Boundaries ---
#     layer_boundaries = [0.0]
#     curr_z = 0.0
#     for L in flat_layers:
#         curr_z += L["t_layer"]
#         layer_boundaries.append(curr_z)
#     z_max = layer_boundaries[-1]

#     # --- 2. Exact Mesh Generation ---
    
#     # Collect exact geometric points to prevent snapping offsets
#     exact_z = []
#     for i, L in enumerate(flat_layers):
#         exact_z.append(layer_boundaries[i])                       # Exact boundary
#         exact_z.append(layer_boundaries[i] + (L["t_layer"]/2.0))  # Exact center
#     exact_z.append(z_max)

#     # Pre-calculate explicit dipole offset points
#     if "dipole_configs" in config:
#         for dip in config["dipole_configs"]:
#             indices = get_layer_indices_by_offset_range(
#                 flat_layers, dip["interface"], dip["direction"], dip["layer_offset"]
#             )
#             for idx in indices:
#                 L = flat_layers[idx]
#                 z_center = layer_boundaries[idx] + (L["t_layer"] / 2.0)
#                 d_plus = dip.get("delta_plus_0", 0.0)
#                 d_minus = dip.get("delta_minus_0", 0.0)
                
#                 if (d_plus + d_minus) > L["t_layer"]:
#                     scale = (L["t_layer"] * 0.99) / (d_plus + d_minus)
#                     d_plus *= scale
#                     d_minus *= scale
                    
#                 exact_z.append(z_center - d_plus)
#                 exact_z.append(z_center + d_minus)

#     # Generate base uniform mesh
#     min_t = min(L["t_layer"] for L in flat_layers)
#     dz_fine = min_t / config["points_per_layer"]
#     Nz_uniform = int(round(z_max / dz_fine)) + 1
#     z_uniform = np.linspace(0, z_max, Nz_uniform)

#     # --- TOLERANCE-BASED MESH MERGING ---
#     # Start with the exact geometric points
#     z_exact = np.unique(np.array(exact_z))
#     z_exact.sort()

#     # Add uniform points only if they are not already sufficiently close to an exact point
#     tol = 1e-15  # meters

#     z_list = list(z_exact)

#     for zu in z_uniform:
#         if np.min(np.abs(z_exact - zu)) > tol:
#             z_list.append(zu)

#     z = np.sort(np.array(z_list))
#     Nz = len(z)

#     # Check minimum element size
#     dz = np.diff(z)
#     print(f"Mesh generated with {Nz} nodes")
#     print(f"Minimum dz = {np.min(dz):.6e} m")
#     print(f"Maximum dz = {np.max(dz):.6e} m")
#     print(f"Target dz  = {dz_fine:.6e} m")

#     # Calculate local effective nodal volumes for 1D density mapping
#     dz_node_arr = np.zeros(Nz)
#     dz_node_arr[0] = (z[1] - z[0]) / 2.0
#     dz_node_arr[-1] = (z[-1] - z[-2]) / 2.0
#     dz_node_arr[1:-1] = (z[2:] - z[:-2]) / 2.0

#     # --- 3. Property Array Initialization ---
#     V0_c, V0_v = np.zeros(Nz), np.zeros(Nz)
#     meff_e_arr, meff_h_arr = np.zeros(Nz), np.zeros(Nz)
#     Epsi_arr = np.zeros(Nz)
#     N_D_array, N_A_array = np.zeros(Nz), np.zeros(Nz)
#     N_D_dipole, N_A_dipole = np.zeros(Nz), np.zeros(Nz)

#     layer_idx_map = np.searchsorted(layer_boundaries, z, side='right') - 1

#     # Now every mesh point knows which physical layer it belongs to.
#     layer_idx_map[layer_idx_map >= len(flat_layers)] = len(flat_layers) - 1

#     # Find the primary interface node for tracking
#     interface_idx = 0

#     # Fill material properties
#     for i in range(Nz - 1):
#         if flat_layers[layer_idx_map[i]]["mat_name"] != flat_layers[layer_idx_map[i + 1]]["mat_name"]:
#             # "Does the material change between these two neighbouring nodes?"
#             interface_idx = i
#             break

#     for k in range(Nz):
#             L = flat_layers[layer_idx_map[k]]
#             V0_c[k] = L["Ec"]
#             V0_v[k] = L["Ev"]
#             meff_e_arr[k] = L["meff_e"]
#             meff_h_arr[k] = L["meff_h"]
#             if L.get("is_eps_dynamic", False) and "eps_dynamic_func" in L:
#                 # Initialize with an electric field of 0.0 V/m
#                 Epsi_arr[k] = L["eps_dynamic_func"](0.0) 
#             else:
#                 Epsi_arr[k] = L["eps_r"]

#     # --- 4. Intrinsic Layer-Based Charges (Polar Planes) ---
#     for idx, L in enumerate(flat_layers):
#         if "layer_charge" in L and L["layer_charge"] != 0.0:
#             val_2d = L["layer_charge"] / L["uc_area"]
#             z_center = layer_boundaries[idx] + (L["t_layer"] / 2.0)
#             idx_center = np.argmin(np.abs(z - z_center))
            
#             val_3d_delta = val_2d / dz_node_arr[idx_center]
            
#             if L["layer_charge"] > 0:
#                 N_D_array[idx_center] += val_3d_delta
#             else:
#                 N_A_array[idx_center] += abs(val_3d_delta)

#     # --- 5. Interface-Relative Doping Injection (Array Support) ---
#     if "doping_profiles" in config:
#         for dop in config["doping_profiles"]:
#             target_arr = N_D_array if dop["type"] == "donor" else N_A_array
#             val_e_per_uc = dop["value"]

#             indices = get_layer_indices_by_offset_range(
#                 flat_layers, dop["interface"], dop["direction"], dop["layer_offset"]
#             )

#             for step, idx in enumerate(indices):
#                 L = flat_layers[idx]
#                 z_center = layer_boundaries[idx] + (L["t_layer"] / 2.0)
#                 idx_center = np.argmin(np.abs(z - z_center))

#                 # Apply the exact 3D mesh density strictly for this node
#                 target_arr[idx_center] += (val_e_per_uc / L["uc_area"]) / dz_node_arr[idx_center]
                
#     # --- 6. Split Charge / Dipole Polarization (Array Support) ---
#     P_layers = []
#     if "dipole_configs" in config:
#         for dip in config["dipole_configs"]:
#             base_val_e_per_uc = dip["value"]

#             indices = get_layer_indices_by_offset_range(
#                 flat_layers, dip["interface"], dip["direction"], dip["layer_offset"]
#             )

#             for step, idx in enumerate(indices):
#                 L = flat_layers[idx]
#                 z_center = layer_boundaries[idx] + (L["t_layer"] / 2.0)

#                 if dip.get("decay_mode") == 'exponential':
#                     decay = np.exp(-step / dip.get("decay_rate", 1.0))
#                 else:
#                     decay = 1.0
                
#                 current_val_e_per_uc = base_val_e_per_uc * decay

#                 d_plus = dip.get("delta_plus_0", 0.0)
#                 d_minus = dip.get("delta_minus_0", 0.0)

#                 if (d_plus + d_minus) > L["t_layer"]:
#                     scale = (L["t_layer"] * 0.99) / (d_plus + d_minus)
#                     d_plus *= scale
#                     d_minus *= scale

#                 z_pos = z_center - d_plus
#                 z_neg = z_center + d_minus

#                 mesh_idx_pos = np.argmin(np.abs(z - z_pos))
#                 mesh_idx_neg = np.argmin(np.abs(z - z_neg))

#                 # Deposit Deltas using local nodal widths
#                 N_D_dipole[mesh_idx_pos] += (current_val_e_per_uc / L["uc_area"]) / dz_node_arr[mesh_idx_pos]
#                 N_A_dipole[mesh_idx_neg] += (current_val_e_per_uc / L["uc_area"]) / dz_node_arr[mesh_idx_neg]

#                 dipole_moment = current_val_e_per_uc * e * (z_pos - z_neg)
#                 P_layers.append(dipole_moment / (L["uc_area"] * L["t_layer"]))

#     # --- 7. FEniCS Initialization ---
#     mesh = df.IntervalMesh(Nz - 1, z[0], z[-1])
#     mesh.coordinates()[:, 0] = z
#     V_space = df.FunctionSpace(mesh, 'CG', 1)

#     alpha = config["alpha"]
#     Vin = np.zeros_like(z)
#     convergence_met = False
#     prev_Err = None

#     print("\nStarting Self-Consistent Poisson-Schrödinger Solver...")

#     try:
#         # --- 8. Self-Consistent Loop ---
#         for nloop in range(1, config["Nloops"] + 1):
#             Vtot_c = V0_c + Vin
#             Vtot_v = V0_v + Vin

#             Ec, psic = Schroed1D_dolfin(V_space, Vtot_c, meff_e_arr, config["n"])
#             Ev_inv, psiv = Schroed1D_dolfin(V_space, -Vtot_v, meff_h_arr, config["n"])
#             Ev = -np.array(Ev_inv) if len(Ev_inv) > 0 else np.array([])

#             if len(Ec) == 0:
#                 break

#             Ef, N_n, P_m, ND_plus_z, NA_minus_z, net_integral_charge = find_Ef_general(
#                 Ec, Ev, psic, psiv, N_D_array, N_A_array, Vtot_c, Vtot_v, z, config["T"],
#                 meff_e_arr, meff_h_arr, config["E_bind_D"], config["E_bind_A"], 
#                 config["fully_ionized_flag"]
#             )

#             N_3D = np.sum(np.tile(N_n, (len(z), 1)) * np.abs(psic)**2, axis=1) if len(Ec) > 0 else np.zeros_like(z)
#             P_3D = np.sum(np.tile(P_m, (len(z), 1)) * np.abs(psiv)**2, axis=1) if len(Ev) > 0 else np.zeros_like(z)

#             rho_net = P_3D - N_3D + ND_plus_z - NA_minus_z + N_D_dipole - N_A_dipole

#             Vout = Poisson_dolfin(
#                 V_space, rho_net, Epsi_arr,
#                 config.get("bc_left_type", "neumann"), config.get("V_left_applied", 0.0),
#                 config.get("bc_right_type", "neumann"), config.get("V_right_applied", 0.0)
#             )

#             E_Vm = np.abs(np.gradient(Vout, z))  # Absolute electric field in V/m
#             Epsi_new = np.zeros_like(Epsi_arr)
#             for k in range(Nz):
#                 L = flat_layers[layer_idx_map[k]]
#                 if L.get("is_eps_dynamic", False) and "eps_dynamic_func" in L:
#                     # The main script doesn't know the physics; it just calls the function from config
#                     Epsi_new[k] = L["eps_dynamic_func"](E_Vm[k]) 
#                 else:
#                     Epsi_new[k] = L["eps_r"]

#             Epsi_arr = alpha * Epsi_new + (1 - alpha) * Epsi_arr
#             Err = np.max(np.abs(Vout - Vin))

#             # Output Formatting
#             Q_e = np.sum(N_n)
#             Q_h = np.sum(P_m)
#             Q_ND = np.trapz(ND_plus_z, z)
#             Q_NA = np.trapz(NA_minus_z, z)
#             Q_net = np.trapz(rho_net, z)

#             L_z = z[-1] - z[0]
#             z0 = 0.5 * (z[0] + z[-1])
#             P_elec = -e * np.trapz((z - z0) * N_3D, x=z) / L_z
#             P_hole = e * np.trapz((z - z0) * P_3D, x=z) / L_z
#             P_dipole_inserted = sum(P_layers)

# # --- ADD THESE LINES RIGHT ABOVE THE PRINTS ---
#             eps_old_int = Epsi_arr[interface_idx]
#             eps_new_int = Epsi_new[interface_idx]
#             max_eps_change = np.max(np.abs(Epsi_new - Epsi_arr))

#             print(f"==========================================================================================")
#             print(f" LOOP {nloop:03d} | Error: {Err:.3e} (Tol: {config['tol']:.1e}) | Ef: {Ef:.5f} eV")
#             print(f"------------------------------------------------------------------------------------------")
#             print(f" 2D Electron Dens (np.sum N_n) : {-Q_e:.4e} m^-2")
#             print(f" 2D Hole Dens (np.sum P_m) : {+Q_h:.4e} m^-2")
#             print(f" 2D Donor Dens (Trapz ND+) : {+Q_ND:.4e} m^-2")
#             print(f" 2D Acceptor Dens (Trapz NA-) : {-Q_NA:.4e} m^-2")
#             print(f" Net Charge to Poisson (Trapz) : {Q_net:.4e} m^-2")
#             print(f" Q Residual (Convergence check): {net_integral_charge:.4e} m^-2")
#             print(f" Interface Dielectric (Old -> New): {eps_old_int:.2f} -> {eps_new_int:.2f}")
#             print(f" Max Dielectric Change across mesh: {max_eps_change:.4f}")
#             print(f"------------------------------------------------------------------------------------------")
#             print(f" Charge Dens @ Interface Node : {rho_net[interface_idx]:.4e} m^-3 (z = {z[interface_idx]*1e9:.2f} nm)")
#             print(f"==========================================================================================\n")

#             if Err < config["tol"]:
#                 convergence_met = True
#                 break

#             if prev_Err is not None:
#                 if Err > prev_Err:
#                     alpha = max(alpha * config["shrink_factor"], config["alpha_min"])
#                 else:
#                     alpha = min(alpha * config["grow_factor"], config["alpha_max"])
#             prev_Err = Err
#             Vin = alpha * Vout + (1 - alpha) * Vin

#     except KeyboardInterrupt:
#         print("\n\n[!] SIMULATION INTERRUPTED BY USER (Ctrl+C)")
#         print(f"[!] Halting at Loop {nloop}. Proceeding to plot the partially converged state...\n")

#     return {
#             "z": z, "Vtot_c": V0_c + Vout, "Vtot_v": V0_v + Vout,
#             "V0_c": V0_c, "V0_v": V0_v, "Ef": Ef,
#             "N_3D": N_3D, "P_3D": P_3D, "Ec": Ec, "Ev": Ev,
#             "psic": psic, "psiv": psiv,
#             "ND_plus_z": ND_plus_z, "NA_minus_z": NA_minus_z,
#             "N_D_dipole": N_D_dipole, "N_A_dipole": N_A_dipole,
#             "rho_net": rho_net, "flat_layers": flat_layers,
#             # ADD THESE TWO LINES:
#             "Epsi_arr": Epsi_arr, 
#             "E_field": E_Vm 
#         }

# def add_layer_background(ax, flat_layers, show_colors=True):
#     import matplotlib.pyplot as plt

#     if show_colors:
#         # Map each unique atomic layer name to a distinct color dynamically
#         unique_layer_names = list(dict.fromkeys(L["name"] for L in flat_layers))
#         cmap = plt.colormaps.get_cmap("tab20")
#         color_map = {name: cmap(i % 20) for i, name in enumerate(unique_layer_names)}

#     curr_z = 0.0
#     prev_mat = None

#     for L in flat_layers:
#         t = L["t_layer"] * 1e9
#         mat_name = L["mat_name"]

#         if show_colors:
#             layer_color = color_map[L["name"]]
#             # The blurred thickness spanning exactly the layer's width
#             ax.axvspan(curr_z, curr_z + t, facecolor=layer_color, alpha=0.25, edgecolor='none', zorder=-100)
#             # The central spike at the atomic center
#             z_center = curr_z + (t / 2.0)
#             ax.axvline(z_center, color=layer_color, alpha=0.6, lw=1.5, zorder=-99)

#         # Material interface boundary marker (drawn in both modes)
#         if prev_mat is not None and mat_name != prev_mat:
#             ax.axvline(curr_z, color='black', lw=2, linestyle='--', zorder=-50)

#         curr_z += t
#         prev_mat = mat_name

# def run_thickness_sweep(base_config):
#     import copy
#     import matplotlib.pyplot as plt
#     import numpy as np
    
#     sweep_cfg = base_config.get("thickness_sweep", {})
#     min_uc, max_uc = sweep_cfg["uc_range"]
#     fixed_stack = sweep_cfg["fixed_stack"]
#     sweep_mat = sweep_cfg["sweep_mat"]
    
#     # STO Lattice constant in cm (approx 3.905 Angstroms)
#     a_sto_cm = 3.905e-8 
#     # In-plane unit cell area at the interface in cm^2
#     uc_area_cm2 = a_sto_cm ** 2 
    
#     uc_array = np.arange(min_uc, max_uc + 1)
#     e_per_uc_array = []
    
#     for uc in uc_array:
#         print(f"\n==================================================")
#         print(f" SWEEP RUN: {sweep_mat} * {uc} Unit Cells")
#         print(f"==================================================")
        
#         # Localize config for this iteration
#         cfg_run = copy.deepcopy(base_config)
#         cfg_run["stack_definition"] = [fixed_stack, f"{sweep_mat}*{uc}"]
        
#         # Suppress FEniCS plotting/outputs for intermediate loops if needed
#         res = run_heterostructure_simulation(cfg_run)
        
#         # Integrate N_3D (m^-3) over z (m) to get 2D sheet density (m^-2)
#         n2d_m2 = np.trapz(res["N_3D"], res["z"])
        
#         # Convert m^-2 to cm^-2
#         n2d_cm2 = n2d_m2 * 1e-4
        
#         # Convert cm^-2 to electrons per unit cell (e-/u.c.)
#         e_per_uc = n2d_cm2 * uc_area_cm2
#         e_per_uc_array.append(e_per_uc)
        
#     # Plotting the Critical Thickness Curve
#     fig, ax = plt.subplots(figsize=(8, 5))
#     ax.plot(uc_array, e_per_uc_array, marker='o', lw=2, color='blue', markersize=8, label='Simulated')
    
#     # Add a horizontal dashed line to show the 0.5 e/uc theoretical limit
#     ax.axhline(0.5, color='red', linestyle='--', label='Polar Catastrophe Limit (0.5 $e^-$/u.c.)')
    
#     ax.set_xlabel(f'{sweep_mat} Thickness (Unit Cells)', fontsize=12)
#     ax.set_ylabel(r'Electron Density ($e^- / u.c.$)', fontsize=12)
#     ax.set_title(f'Critical Thickness Sweep: {fixed_stack} / {sweep_mat}', fontsize=12)
#     ax.grid(True, alpha=0.4)
#     ax.legend()
    
#     # Removed the log scale logic since e/uc is best viewed on a linear scale 
#     # to see it asymptotically approach 0.5
        
#     fig.tight_layout()
#     plt.show()

# # def run_thickness_sweep(base_config):
# #     import copy
# #     import matplotlib.pyplot as plt
    
# #     sweep_cfg = base_config.get("thickness_sweep", {})
# #     min_uc, max_uc = sweep_cfg["uc_range"]
# #     fixed_stack = sweep_cfg["fixed_stack"]
# #     sweep_mat = sweep_cfg["sweep_mat"]
    
# #     uc_array = np.arange(min_uc, max_uc + 1)
# #     n2d_array = []
    
# #     for uc in uc_array:
# #         print(f"\n==================================================")
# #         print(f" SWEEP RUN: {sweep_mat} * {uc} Unit Cells")
# #         print(f"==================================================")
        
# #         # Localize config for this iteration
# #         cfg_run = copy.deepcopy(base_config)
# #         cfg_run["stack_definition"] = [fixed_stack, f"{sweep_mat}*{uc}"]
        
# #         # Suppress FEniCS plotting/outputs for intermediate loops if needed, 
# #         # though your solver prints are handled internally.
# #         res = run_heterostructure_simulation(cfg_run)
        
# #         # Integrate N_3D (m^-3) over z (m) to get 2D sheet density (m^-2)
# #         n2d_m2 = np.trapz(res["N_3D"], res["z"])
        
# #         # Convert m^-2 to cm^-2
# #         n2d_cm2 = n2d_m2 * 1e-4
# #         n2d_array.append(n2d_cm2)
        
# #     # Plotting the Critical Thickness Curve
# #     fig, ax = plt.subplots(figsize=(8, 5))
# #     ax.plot(uc_array, n2d_array, marker='o', lw=2, color='blue', markersize=8)
    
# #     ax.set_xlabel(f'{sweep_mat} Thickness (Unit Cells)', fontsize=12)
# #     ax.set_ylabel(r'2D Electron Density ($cm^{-2}$)', fontsize=12)
# #     ax.set_title(f'Critical Thickness Sweep: {fixed_stack} / {sweep_mat}', fontsize=12)
# #     ax.grid(True, alpha=0.4)
    
# #     # Optional: Use log scale if density spans multiple orders of magnitude
# #     if max(n2d_array) > 0 and min(n2d_array) > 0 and (max(n2d_array) / min(n2d_array) > 100):
# #         ax.set_yscale('log')
        
# #     fig.tight_layout()
# #     plt.show()


# if __name__ == "__main__":
#     res = run_heterostructure_simulation(config)

#     z, Vtot_c, Vtot_v = res["z"], res["Vtot_c"], res["Vtot_v"]
#     V0_c, V0_v, Ef = res["V0_c"], res["V0_v"], res["Ef"]
#     Ec, Ev, psic, psiv = res["Ec"], res["Ev"], res["psic"], res["psiv"]
#     N_3D, P_3D, rho_net = res["N_3D"], res["P_3D"], res["rho_net"]
#     ND_plus_z, NA_minus_z = res["ND_plus_z"], res["NA_minus_z"]
#     flat_layers = res["flat_layers"]

#     ScF = 0.05 # Scaling factor for visualization heights


#     # Check if sweep mode is enabled
#     if config.get("thickness_sweep", {}).get("enabled", False):
#         run_thickness_sweep(config)
        
#     # else:
#     #     # --- EXISTING SINGLE RUN LOGIC ---
#     #     res = run_heterostructure_simulation(config)
    
#     #     z, Vtot_c, Vtot_v = res["z"], res["Vtot_c"], res["Vtot_v"]
#     #     V0_c, V0_v, Ef = res["V0_c"], res["V0_v"], res["Ef"]
#     #     Ec, Ev, psic, psiv = res["Ec"], res["Ev"], res["psic"], res["psiv"]
#     #     N_3D, P_3D, rho_net = res["N_3D"], res["P_3D"], res["rho_net"]
#     #     ND_plus_z, NA_minus_z = res["ND_plus_z"], res["NA_minus_z"]
#     #     flat_layers = res["flat_layers"]
    
#     #     ScF = 0.05 # Scaling factor for visualization heights
        
#     #     # ... [Keep all of your existing Plot 1 through Plot 6 logic here] ...
        
#     #     plt.show()

#     # =============================================================
#     # SET 1: ORIGINAL PLOTS (WITH COLORS & HEATMAPS)
#     # =============================================================

#     # --- PLOT 1: Original Band Diagram + Wavefunctions ---
#     fig1, ax1 = plt.subplots(figsize=(10, 6))
#     add_layer_background(ax1, flat_layers, show_colors=True)

#     ax1.plot(z*1e9, V0_c, 'b--', alpha=0.5, lw=2, label='Ec (Flat)', zorder=10)
#     ax1.plot(z*1e9, V0_v, 'r--', alpha=0.5, lw=2, label='Ev (Flat)', zorder=10)
#     ax1.plot(z*1e9, Vtot_c, 'b', lw=2, label='Ec (Bent)', zorder=10)
#     ax1.plot(z*1e9, Vtot_v, 'r', lw=2, label='Ev (Bent)', zorder=10)
#     ax1.axhline(Ef, color='lime', ls='--', lw=2, label=r'$E_f$', zorder=10)

#     if len(Ec) > 0:
#         for i in range(min(config["n"], len(Ec))):
#             intensity = np.abs(psic[:, i])**2
#             if np.max(intensity) > 1e-15:
#                 norm_intensity = intensity / np.max(intensity)
#                 ax1.plot(z*1e9, np.full(len(z), Ec[i]) + (ScF * norm_intensity), 
#                          color='blue', linestyle=':', lw=1.5, alpha=0.8, zorder=5)
#                 ax1.axhline(Ec[i], color='blue', linestyle='-', lw=0.5, alpha=0.3, zorder=4)

#     if len(Ev) > 0:
#         for i in range(min(config["n"], len(Ev))):
#             intensity = np.abs(psiv[:, i])**2
#             if np.max(intensity) > 1e-15:
#                 norm_intensity = intensity / np.max(intensity)
#                 ax1.plot(z*1e9, np.full(len(z), Ev[i]) - (ScF * norm_intensity), 
#                          color='red', linestyle=':', lw=1.5, alpha=0.8, zorder=5)
#                 ax1.axhline(Ev[i], color='red', linestyle='-', lw=0.5, alpha=0.3, zorder=4)

#     ax1.set(xlim=(0, z[-1]*1e9), ylim=(np.min(Vtot_v) - 0.2, np.max(Vtot_c) + 0.2),
#             xlabel='z (nm)', ylabel='Energy (eV)', title='[Original] Band Diagram & Individual Wavefunctions')
#     ax1.grid(alpha=0.3)
#     ax1.legend(loc='center right')
#     fig1.tight_layout()

#     # --- PLOT 2: Original Band Diagram + Total Density Heatmap ---
#     fig2, ax2 = plt.subplots(figsize=(10, 6))
#     add_layer_background(ax2, flat_layers, show_colors=True)

#     ax2.plot(z*1e9, Vtot_c, 'b', lw=2, label='Ec (Bent)', zorder=10)
#     ax2.plot(z*1e9, Vtot_v, 'r', lw=2, label='Ev (Bent)', zorder=10)
#     ax2.axhline(Ef, color='lime', ls='--', lw=2, label=r'$E_f$', zorder=10)

#     if np.max(N_3D) > 1e10:
#         norm_N3D = N_3D / np.max(N_3D)
#         mask_e = norm_N3D > 1e-3
#         ax2.scatter(z[mask_e]*1e9, Vtot_c[mask_e] + (ScF * norm_N3D[mask_e]), 
#                     c=norm_N3D[mask_e], cmap='Blues', s=25, alpha=0.9, edgecolor='none', zorder=5, label='Electron Heatmap')

#     if np.max(P_3D) > 1e10:
#         norm_P3D = P_3D / np.max(P_3D)
#         mask_h = norm_P3D > 1e-3
#         ax2.scatter(z[mask_h]*1e9, Vtot_v[mask_h] - (ScF * norm_P3D[mask_h]), 
#                     c=norm_P3D[mask_h], cmap='Reds', s=25, alpha=0.9, edgecolor='none', zorder=5, label='Hole Heatmap')

#     ax2.set(xlim=(0, z[-1]*1e9), ylim=(np.min(Vtot_v) - 0.2, np.max(Vtot_c) + 0.2),
#             xlabel='z (nm)', ylabel='Energy (eV)', title='[Original] Carrier Spatial Distribution Heatmap')
#     ax2.grid(alpha=0.3)
#     handles, labels = ax2.get_legend_handles_labels()
#     by_label = dict(zip(labels, handles))
#     ax2.legend(by_label.values(), by_label.keys(), loc='center right')
#     fig2.tight_layout()

#     # --- PLOT 3: Original Charge Distribution ---
#     fig3, ax3 = plt.subplots(figsize=(10, 5))
#     add_layer_background(ax3, flat_layers, show_colors=True)

#     ax3.plot(z*1e9, -N_3D, lw=2, color='blue', label=r'$-N_{3D}$ (Electrons)')
#     ax3.plot(z*1e9, P_3D, lw=2, color='red', label=r'$P_{3D}$ (Holes)')
#     ax3.plot(z*1e9, rho_net, lw=2, color='black', alpha=0.8, label=r'$\rho_{net}$')

#     fixed_pos = ND_plus_z + res["N_D_dipole"]
#     fixed_neg = NA_minus_z + res["N_A_dipole"]
    
#     pos_idx = np.where(fixed_pos > 1e10)[0] 
#     neg_idx = np.where(fixed_neg > 1e10)[0]

#     if len(pos_idx) > 0:
#         ax3.vlines(z[pos_idx]*1e9, 0, fixed_pos[pos_idx], color='darkred', linestyle='-', lw=2, label=r'Fixed $+$ (Delta)')
#     if len(neg_idx) > 0:
#         ax3.vlines(z[neg_idx]*1e9, 0, -fixed_neg[neg_idx], color='darkblue', linestyle='-', lw=2, label=r'Fixed $-$ (Delta)')

#     ax3.set_yscale('symlog', linthresh=1e25)
#     ax3.axhline(0, color='gray', ls='--', lw=1)
#     ax3.set(xlabel='z (nm)', ylabel=r'Charge Density (m$^{-3}$)', title='[Original] Charge Distribution Profile')
#     ax3.grid(alpha=0.3)
#     ax3.legend(loc='upper right', bbox_to_anchor=(1.25, 1.0))
#     fig3.tight_layout()


#     # =============================================================
#     # SET 2: CLEAN PLOTS (NO COLORS, 1D DENSITY INSTEAD OF HEATMAP)
#     # =============================================================

# ############################################################################################################################


#     # --- PLOT 4: Clean Conduction Band + Electron Wavefunctions ---
#     fig4, ax4 = plt.subplots(figsize=(10, 6))
#     add_layer_background(ax4, flat_layers, show_colors=False)

#     # Conduction band
#     ax4.plot(
#         z*1e9, V0_c,
#         'b--', alpha=0.5, lw=2,
#         label='Ec (Flat)', zorder=10
#     )

#     ax4.plot(
#         z*1e9, Vtot_c,
#         'b', lw=2,
#         label='Ec (Bent)', zorder=10
#     )

#     # Electron wavefunctions
#     if len(Ec) > 0:
#         for i in range(min(config["n"], len(Ec))):
#             intensity = np.abs(psic[:, i])**2

#             if np.max(intensity) > 1e-15:
#                 norm_intensity = intensity / np.max(intensity)

#                 ax4.plot(
#                     z*1e9,
#                     np.full(len(z), Ec[i]) + (ScF * norm_intensity),
#                     color='blue',
#                     linestyle='-',
#                     lw=1.5,
#                     alpha=0.8,
#                     zorder=5
#                 )

#                 ax4.axhline(
#                     Ec[i],
#                     color='blue',
#                     linestyle='-',
#                     lw=0.5,
#                     alpha=0.3,
#                     zorder=4
#                 )

#     # Fixed viewing window
#     ax4.set(
#         xlim=(20, 40),
#         ylim=(-4, 0),
#         xlabel='z (nm)',
#         ylabel='Energy (eV)',
#         title='Conduction Band and Electron Wavefunctions'
#     )

#     ax4.grid(alpha=0.3)
#     ax4.legend(loc='center right')
#     fig4.tight_layout()


#     # # --- PLOT 4: Clean Band Diagram + Wavefunctions ---
#     # fig4, ax4 = plt.subplots(figsize=(10, 6))
#     # add_layer_background(ax4, flat_layers, show_colors=False)

#     # ax4.plot(z*1e9, V0_c, 'b--', alpha=0.5, lw=2, label='Ec (Flat)', zorder=10)
#     # ax4.plot(z*1e9, V0_v, 'r--', alpha=0.5, lw=2, label='Ev (Flat)', zorder=10)
#     # ax4.plot(z*1e9, Vtot_c, 'b', lw=2, label='Ec (Bent)', zorder=10)
#     # ax4.plot(z*1e9, Vtot_v, 'r', lw=2, label='Ev (Bent)', zorder=10)
#     # ax4.axhline(Ef, color='lime', ls='--', lw=2, label=r'$E_f$', zorder=10)

#     # if len(Ec) > 0:
#     #     for i in range(min(config["n"], len(Ec))):
#     #         intensity = np.abs(psic[:, i])**2
#     #         if np.max(intensity) > 1e-15:
#     #             norm_intensity = intensity / np.max(intensity)
#     #             ax4.plot(z*1e9, np.full(len(z), Ec[i]) + (ScF * norm_intensity), 
#     #                      color='blue', linestyle='-', lw=1.5, alpha=0.8, zorder=5)
#     #             ax4.axhline(Ec[i], color='blue', linestyle='-', lw=0.5, alpha=0.3, zorder=4)

#     # if len(Ev) > 0:
#     #     for i in range(min(config["n"], len(Ev))):
#     #         intensity = np.abs(psiv[:, i])**2
#     #         if np.max(intensity) > 1e-15:
#     #             norm_intensity = intensity / np.max(intensity)
#     #             ax4.plot(z*1e9, np.full(len(z), Ev[i]) - (ScF * norm_intensity), 
#     #                      color='red', linestyle='-', lw=1.5, alpha=0.8, zorder=5)
#     #             ax4.axhline(Ev[i], color='red', linestyle='-', lw=0.5, alpha=0.3, zorder=4)

#     # # ax4.set(xlim=(0, z[-1]*1e9), ylim=(np.min(Vtot_v) - 0.2, np.max(Vtot_c) + 0.2),
#     # #         xlabel='z (nm)', ylabel='Energy (eV)', title='[Clean] Band Diagram & Individual Wavefunctions')
#     # # Fixed viewing window
#     # ax4.set(
#     #     xlim=(20, 40),
#     #     ylim=(-4, 0),
#     #     xlabel='z (nm)',
#     #     ylabel='Energy (eV)',
#     #     title='[Clean] Band Diagram & Individual Wavefunctions'
#     # )

#     # ax4.grid(alpha=0.3)
#     # ax4.legend(loc='center right')
#     # fig4.tight_layout()

# ############################################################################################################################

#     # --- PLOT 5: Clean Band Diagram + 1D Carrier Density (Replaces Heatmap) ---
#     fig5, ax5 = plt.subplots(figsize=(10, 6))
#     add_layer_background(ax5, flat_layers, show_colors=False)

#     ax5.plot(z*1e9, Vtot_c, 'b', lw=2, label='Ec (Bent)', zorder=10)
#     ax5.plot(z*1e9, Vtot_v, 'r', lw=2, label='Ev (Bent)', zorder=10)
#     ax5.axhline(Ef, color='lime', ls='--', lw=2, label=r'$E_f$', zorder=10)

#     ax5.set(xlim=(0, z[-1]*1e9), ylim=(np.min(Vtot_v) - 0.2, np.max(Vtot_c) + 0.2),
#             xlabel='z (nm)', ylabel='Energy (eV)', title='[Clean] Band Diagram & 1D Carrier Density')
#     ax5.grid(alpha=0.3)
    
#     # Secondary Y-axis for 1D Density Plot
#     ax5_twin = ax5.twinx()
    
#     if np.max(N_3D) > 1e10:
#         ax5_twin.plot(z*1e9, N_3D, color='darkblue', lw=2, label=r'Electron Density ($N_{3D}$)', zorder=15)
#     if np.max(P_3D) > 1e10:
#         ax5_twin.plot(z*1e9, P_3D, color='darkred', lw=2, label=r'Hole Density ($P_{3D}$)', zorder=15)

#     ax5_twin.set_ylabel(r'Carrier Density (m$^{-3}$)')
#     ax5_twin.set_yscale('symlog', linthresh=1e18) 
    
#     handles1, labels1 = ax5.get_legend_handles_labels()
#     handles2, labels2 = ax5_twin.get_legend_handles_labels()
#     by_label5 = dict(zip(labels1 + labels2, handles1 + handles2))
#     ax5.legend(by_label5.values(), by_label5.keys(), loc='center right')
#     fig5.tight_layout()

#     # --- PLOT 6: Clean Charge Distribution ---
#     fig6, ax6 = plt.subplots(figsize=(10, 5))
#     add_layer_background(ax6, flat_layers, show_colors=False)

#     ax6.plot(z*1e9, -N_3D, lw=2, color='blue', label=r'$-N_{3D}$ (Electrons)')
#     ax6.plot(z*1e9, P_3D, lw=2, color='red', label=r'$P_{3D}$ (Holes)')
#     ax6.plot(z*1e9, rho_net, lw=2, color='black', alpha=0.8, label=r'$\rho_{net}$')

#     if len(pos_idx) > 0:
#         ax6.vlines(z[pos_idx]*1e9, 0, fixed_pos[pos_idx], color='darkred', linestyle='-', lw=2, label=r'Fixed $+$ (Delta)')
#     if len(neg_idx) > 0:
#         ax6.vlines(z[neg_idx]*1e9, 0, -fixed_neg[neg_idx], color='darkblue', linestyle='-', lw=2, label=r'Fixed $-$ (Delta)')

#     ax6.set_yscale('symlog', linthresh=1e25)
#     ax6.axhline(0, color='gray', ls='--', lw=1)
#     ax6.set(xlabel='z (nm)', ylabel=r'Charge Density (m$^{-3}$)', title='[Clean] Charge Distribution Profile')
#     ax6.grid(alpha=0.3)
#     ax6.legend(loc='upper right', bbox_to_anchor=(1.25, 1.0))
#     fig6.tight_layout()

# # =============================================================
#     # --- PLOT 7: Electric Field and Dynamic Dielectric Constant ---
#     # =============================================================
#     Epsi_final = res["Epsi_arr"]
#     # Convert E_field from V/m to MV/cm for easier reading on the plot
#     E_final_MV_cm = res["E_field"] * 1e-8 

#     fig7, ax7 = plt.subplots(figsize=(10, 6))
#     add_layer_background(ax7, flat_layers, show_colors=False)
    
#     # Plot Electric Field on the primary y-axis
#     ax7.plot(z*1e9, E_final_MV_cm, color='purple', lw=2, label='Electric Field (MV/cm)', zorder=10)
#     ax7.set(xlabel='z (nm)', ylabel='Electric Field (MV/cm)', title='Electric Field & Dielectric Profile vs Space')
#     ax7.grid(alpha=0.3)
    
#     # Plot Dielectric Constant on a secondary y-axis
#     ax7_twin = ax7.twinx()
#     ax7_twin.plot(z*1e9, Epsi_final, color='orange', lw=2, linestyle='--', label=r'Dielectric Constant ($\epsilon_r$)', zorder=15)
#     ax7_twin.set_ylabel(r'Relative Permittivity ($\epsilon_r$)')
    
#     # Combine legends
#     handles1, labels1 = ax7.get_legend_handles_labels()
#     handles2, labels2 = ax7_twin.get_legend_handles_labels()
#     by_label7 = dict(zip(labels1 + labels2, handles1 + handles2))
#     ax7.legend(by_label7.values(), by_label7.keys(), loc='upper right')
    
#     fig7.tight_layout()

    
#     plt.show()

import numpy as np
import dolfin as df
import warnings
import matplotlib.pyplot as plt

# Import from your local files
from master_helper import find_Ef_general, Schroed1D_dolfin, Poisson_dolfin

# from config_BBOiLNB import config
# from config_SIrOiLNB import config
# from config_LAOiSTO import config
# from config_STO100iSTO111 import config
# from config_LLOiSBO import config
from config_STOiSAP import config
# from config_STOiGTO import config

warnings.filterwarnings("ignore")
df.set_log_level(df.LogLevel.WARNING)



class AndersonMixer:
    def __init__(self, m=6, beta=0.1, max_step=0.1, reg=1e-10):
        self.m, self.beta, self.max_step, self.reg = m, beta, max_step, reg
        self.X, self.F = [], []

    def reset(self):
        self.X, self.F = [], []

    def mix(self, x, fx):
        f = fx - x                                   # residual
        self.X.append(x.copy()); self.F.append(f.copy())
        if len(self.X) > self.m + 1:
            self.X.pop(0); self.F.pop(0)

        x_new = x + self.beta * f                    # plain step on first iteration
        if len(self.F) > 1:
            dF = np.column_stack([self.F[i+1] - self.F[i] for i in range(len(self.F)-1)])
            dX = np.column_stack([self.X[i+1] - self.X[i] for i in range(len(self.X)-1)])
            A = dF.T @ dF
            A += self.reg * np.trace(A) / A.shape[0] * np.eye(A.shape[0])   # ridge
            gamma = np.linalg.solve(A, dF.T @ f)
            x_new = x + self.beta * f - (dX + self.beta * dF) @ gamma

        step = x_new - x                             # cap the largest potential change
        smax = np.max(np.abs(step))
        if smax > self.max_step:
            x_new = x + step * (self.max_step / smax)
        return x_new

def run_heterostructure_simulation(config):
    e = 1.602176487E-19  # e = 1.602176487×10^−19 C

    # --- Helper Function Moved Here ---
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

    # --- 1. Expand Stack Definition ---
    flat_layers = []  # will contain every atomic layer like ababababcdcdcdcdcd
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

    # --- Pre-calculate Layer Boundaries ---
    layer_boundaries = [0.0]
    curr_z = 0.0
    for L in flat_layers:
        curr_z += L["t_layer"]
        layer_boundaries.append(curr_z)
    z_max = layer_boundaries[-1]

    # --- 2. Exact Mesh Generation ---
    
    # Collect exact geometric points to prevent snapping offsets
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
                z_center = layer_boundaries[idx] + (L["t_layer"] / 2.0)
                d_plus = dip.get("delta_plus_0", 0.0)
                d_minus = dip.get("delta_minus_0", 0.0)
                
                if (d_plus + d_minus) > L["t_layer"]:
                    scale = (L["t_layer"] * 0.99) / (d_plus + d_minus)
                    d_plus *= scale
                    d_minus *= scale
                    
                exact_z.append(z_center - d_plus)
                exact_z.append(z_center + d_minus)

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
    tol = 1e-15  # meters

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

    # Calculate local effective nodal volumes for 1D density mapping
    dz_node_arr = np.zeros(Nz)
    dz_node_arr[0] = (z[1] - z[0]) / 2.0
    dz_node_arr[-1] = (z[-1] - z[-2]) / 2.0
    dz_node_arr[1:-1] = (z[2:] - z[:-2]) / 2.0

    # --- 3. Property Array Initialization ---
    V0_c, V0_v = np.zeros(Nz), np.zeros(Nz)
    meff_e_arr, meff_h_arr = np.zeros(Nz), np.zeros(Nz)
    Epsi_arr = np.zeros(Nz)
    N_D_array, N_A_array = np.zeros(Nz), np.zeros(Nz)
    N_D_dipole, N_A_dipole = np.zeros(Nz), np.zeros(Nz)

    layer_idx_map = np.searchsorted(layer_boundaries, z, side='right') - 1

    # Now every mesh point knows which physical layer it belongs to.
    layer_idx_map[layer_idx_map >= len(flat_layers)] = len(flat_layers) - 1

    # Find the primary interface node for tracking
    interface_idx = 0

    # Fill material properties
    for i in range(Nz - 1):
        if flat_layers[layer_idx_map[i]]["mat_name"] != flat_layers[layer_idx_map[i + 1]]["mat_name"]:
            # "Does the material change between these two neighbouring nodes?"
            interface_idx = i
            break

    for k in range(Nz):
            L = flat_layers[layer_idx_map[k]]
            V0_c[k] = L["Ec"]
            V0_v[k] = L["Ev"]
            meff_e_arr[k] = L["meff_e"]
            meff_h_arr[k] = L["meff_h"]
            if L.get("is_eps_dynamic", False) and "eps_dynamic_func" in L:
                # Initialize with an electric field of 0.0 V/m
                Epsi_arr[k] = L["eps_dynamic_func"](0.0) 
            else:
                Epsi_arr[k] = L["eps_r"]

    # --- 4. Intrinsic Layer-Based Charges (Polar Planes) ---
    for idx, L in enumerate(flat_layers):
        if "layer_charge" in L and L["layer_charge"] != 0.0:
            val_2d = L["layer_charge"] / L["uc_area"]
            z_center = layer_boundaries[idx] + (L["t_layer"] / 2.0)
            idx_center = np.argmin(np.abs(z - z_center))
            
            val_3d_delta = val_2d / dz_node_arr[idx_center]
            
            if L["layer_charge"] > 0:
                N_D_array[idx_center] += val_3d_delta
            else:
                N_A_array[idx_center] += abs(val_3d_delta)

    # --- 5. Interface-Relative Doping Injection (Array Support) ---
    if "doping_profiles" in config:
        for dop in config["doping_profiles"]:
            target_arr = N_D_array if dop["type"] == "donor" else N_A_array
            val_e_per_uc = dop["value"]

            indices = get_layer_indices_by_offset_range(
                flat_layers, dop["interface"], dop["direction"], dop["layer_offset"]
            )

            for step, idx in enumerate(indices):
                L = flat_layers[idx]
                z_center = layer_boundaries[idx] + (L["t_layer"] / 2.0)
                idx_center = np.argmin(np.abs(z - z_center))

                # Apply the exact 3D mesh density strictly for this node
                target_arr[idx_center] += (val_e_per_uc / L["uc_area"]) / dz_node_arr[idx_center]
                
    # --- 6. Split Charge / Dipole Polarization (Array Support) ---
    P_layers = []
    if "dipole_configs" in config:
        for dip in config["dipole_configs"]:
            base_val_e_per_uc = dip["value"]

            indices = get_layer_indices_by_offset_range(
                flat_layers, dip["interface"], dip["direction"], dip["layer_offset"]
            )

            for step, idx in enumerate(indices):
                L = flat_layers[idx]
                z_center = layer_boundaries[idx] + (L["t_layer"] / 2.0)

                if dip.get("decay_mode") == 'exponential':
                    decay = np.exp(-step / dip.get("decay_rate", 1.0))
                else:
                    decay = 1.0
                
                current_val_e_per_uc = base_val_e_per_uc * decay

                d_plus = dip.get("delta_plus_0", 0.0)
                d_minus = dip.get("delta_minus_0", 0.0)

                if (d_plus + d_minus) > L["t_layer"]:
                    scale = (L["t_layer"] * 0.99) / (d_plus + d_minus)
                    d_plus *= scale
                    d_minus *= scale

                z_pos = z_center - d_plus
                z_neg = z_center + d_minus

                mesh_idx_pos = np.argmin(np.abs(z - z_pos))
                mesh_idx_neg = np.argmin(np.abs(z - z_neg))

                # Deposit Deltas using local nodal widths
                N_D_dipole[mesh_idx_pos] += (current_val_e_per_uc / L["uc_area"]) / dz_node_arr[mesh_idx_pos]
                N_A_dipole[mesh_idx_neg] += (current_val_e_per_uc / L["uc_area"]) / dz_node_arr[mesh_idx_neg]

                dipole_moment = current_val_e_per_uc * e * (z_pos - z_neg)
                P_layers.append(dipole_moment / (L["uc_area"] * L["t_layer"]))

# --- 7. FEniCS Initialization ---
    mesh = df.IntervalMesh(Nz - 1, z[0], z[-1])
    mesh.coordinates()[:, 0] = z
    V_space = df.FunctionSpace(mesh, 'CG', 1)

    mixer = AndersonMixer(m=config["mix_m"], beta=config["mix_beta"], max_step=config["mix_max_step"])
    beta_eps = config["mix_beta_eps"]
    Vin = np.zeros_like(z)
    convergence_met = False
    prev_Err = None

    print("\nStarting Self-Consistent Poisson-Schrödinger Solver...")

    try:
        # --- 8. Self-Consistent Loop ---
        for nloop in range(1, config["Nloops"] + 1):
            Vtot_c = V0_c + Vin
            Vtot_v = V0_v + Vin

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

            E_Vm = np.abs(np.gradient(Vout, z))  # Absolute electric field in V/m
            Epsi_new = np.zeros_like(Epsi_arr)
            for k in range(Nz):
                L = flat_layers[layer_idx_map[k]]
                if L.get("is_eps_dynamic", False) and "eps_dynamic_func" in L:
                    # The main script doesn't know the physics; it just calls the function from config
                    Epsi_new[k] = L["eps_dynamic_func"](E_Vm[k]) 
                else:
                    Epsi_new[k] = L["eps_r"]

# ε update: its own fixed damping, delayed so V settles first
            if nloop > config["eps_start"]:
                eps_err = np.max(np.abs(Epsi_new - Epsi_arr) / Epsi_arr)
                Epsi_arr = beta_eps * Epsi_new + (1 - beta_eps) * Epsi_arr
            else:
                eps_err = 0.0

            Err = np.max(np.abs(Vout - Vin))
            # Output Formatting
            Q_e = np.sum(N_n)
            Q_h = np.sum(P_m)
            Q_ND = np.trapz(ND_plus_z, z)
            Q_NA = np.trapz(NA_minus_z, z)
            Q_net = np.trapz(rho_net, z)

            L_z = z[-1] - z[0]
            z0 = 0.5 * (z[0] + z[-1])
            P_elec = -e * np.trapz((z - z0) * N_3D, x=z) / L_z
            P_hole = e * np.trapz((z - z0) * P_3D, x=z) / L_z
            P_dipole_inserted = sum(P_layers)

# --- ADD THESE LINES RIGHT ABOVE THE PRINTS ---
            eps_old_int = Epsi_arr[interface_idx]
            eps_new_int = Epsi_new[interface_idx]
            max_eps_change = np.max(np.abs(Epsi_new - Epsi_arr))

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

            if Err < config["tol"] and eps_err < config["eps_tol"]:
                convergence_met = True
                break

            if prev_Err is not None and Err > 5.0 * prev_Err:    # blow-up -> drop stale history
                mixer.reset()
            prev_Err = Err

            Vin = mixer.mix(Vin, Vout)

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
            # ADD THESE TWO LINES:
            "Epsi_arr": Epsi_arr, 
            "E_field": E_Vm 
        }

def add_layer_background(ax, flat_layers, show_colors=True):
    import matplotlib.pyplot as plt

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

def run_thickness_sweep(base_config):
    import copy
    import matplotlib.pyplot as plt
    import numpy as np
    
    sweep_cfg = base_config.get("thickness_sweep", {})
    min_uc, max_uc = sweep_cfg["uc_range"]
    fixed_stack = sweep_cfg["fixed_stack"]
    sweep_mat = sweep_cfg["sweep_mat"]
    
    # STO Lattice constant in cm (approx 3.905 Angstroms)
    a_sto_cm = 3.905e-8 
    # In-plane unit cell area at the interface in cm^2
    uc_area_cm2 = a_sto_cm ** 2 
    
    uc_array = np.arange(min_uc, max_uc + 1)
    e_per_uc_array = []
    
    for uc in uc_array:
        print(f"\n==================================================")
        print(f" SWEEP RUN: {sweep_mat} * {uc} Unit Cells")
        print(f"==================================================")
        
        # Localize config for this iteration
        cfg_run = copy.deepcopy(base_config)
        cfg_run["stack_definition"] = [fixed_stack, f"{sweep_mat}*{uc}"]
        
        # Suppress FEniCS plotting/outputs for intermediate loops if needed
        res = run_heterostructure_simulation(cfg_run)
        
        # Integrate N_3D (m^-3) over z (m) to get 2D sheet density (m^-2)
        n2d_m2 = np.trapz(res["N_3D"], res["z"])
        
        # Convert m^-2 to cm^-2
        n2d_cm2 = n2d_m2 * 1e-4
        
        # Convert cm^-2 to electrons per unit cell (e-/u.c.)
        e_per_uc = n2d_cm2 * uc_area_cm2
        e_per_uc_array.append(e_per_uc)
        
    # Plotting the Critical Thickness Curve
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(uc_array, e_per_uc_array, marker='o', lw=2, color='blue', markersize=8, label='Simulated')
    
    # Add a horizontal dashed line to show the 0.5 e/uc theoretical limit
    ax.axhline(0.5, color='red', linestyle='--', label='Polar Catastrophe Limit (0.5 $e^-$/u.c.)')
    
    ax.set_xlabel(f'{sweep_mat} Thickness (Unit Cells)', fontsize=12)
    ax.set_ylabel(r'Electron Density ($e^- / u.c.$)', fontsize=12)
    ax.set_title(f'Critical Thickness Sweep: {fixed_stack} / {sweep_mat}', fontsize=12)
    ax.grid(True, alpha=0.4)
    ax.legend()
    
    # Removed the log scale logic since e/uc is best viewed on a linear scale 
    # to see it asymptotically approach 0.5
        
    fig.tight_layout()
    plt.show()

# def run_thickness_sweep(base_config):
#     import copy
#     import matplotlib.pyplot as plt
    
#     sweep_cfg = base_config.get("thickness_sweep", {})
#     min_uc, max_uc = sweep_cfg["uc_range"]
#     fixed_stack = sweep_cfg["fixed_stack"]
#     sweep_mat = sweep_cfg["sweep_mat"]
    
#     uc_array = np.arange(min_uc, max_uc + 1)
#     n2d_array = []
    
#     for uc in uc_array:
#         print(f"\n==================================================")
#         print(f" SWEEP RUN: {sweep_mat} * {uc} Unit Cells")
#         print(f"==================================================")
        
#         # Localize config for this iteration
#         cfg_run = copy.deepcopy(base_config)
#         cfg_run["stack_definition"] = [fixed_stack, f"{sweep_mat}*{uc}"]
        
#         # Suppress FEniCS plotting/outputs for intermediate loops if needed, 
#         # though your solver prints are handled internally.
#         res = run_heterostructure_simulation(cfg_run)
        
#         # Integrate N_3D (m^-3) over z (m) to get 2D sheet density (m^-2)
#         n2d_m2 = np.trapz(res["N_3D"], res["z"])
        
#         # Convert m^-2 to cm^-2
#         n2d_cm2 = n2d_m2 * 1e-4
#         n2d_array.append(n2d_cm2)
        
#     # Plotting the Critical Thickness Curve
#     fig, ax = plt.subplots(figsize=(8, 5))
#     ax.plot(uc_array, n2d_array, marker='o', lw=2, color='blue', markersize=8)
    
#     ax.set_xlabel(f'{sweep_mat} Thickness (Unit Cells)', fontsize=12)
#     ax.set_ylabel(r'2D Electron Density ($cm^{-2}$)', fontsize=12)
#     ax.set_title(f'Critical Thickness Sweep: {fixed_stack} / {sweep_mat}', fontsize=12)
#     ax.grid(True, alpha=0.4)
    
#     # Optional: Use log scale if density spans multiple orders of magnitude
#     if max(n2d_array) > 0 and min(n2d_array) > 0 and (max(n2d_array) / min(n2d_array) > 100):
#         ax.set_yscale('log')
        
#     fig.tight_layout()
#     plt.show()


if __name__ == "__main__":
    res = run_heterostructure_simulation(config)

    z, Vtot_c, Vtot_v = res["z"], res["Vtot_c"], res["Vtot_v"]
    V0_c, V0_v, Ef = res["V0_c"], res["V0_v"], res["Ef"]
    Ec, Ev, psic, psiv = res["Ec"], res["Ev"], res["psic"], res["psiv"]
    N_3D, P_3D, rho_net = res["N_3D"], res["P_3D"], res["rho_net"]
    ND_plus_z, NA_minus_z = res["ND_plus_z"], res["NA_minus_z"]
    flat_layers = res["flat_layers"]

    ScF = 0.05 # Scaling factor for visualization heights


    # Check if sweep mode is enabled
    if config.get("thickness_sweep", {}).get("enabled", False):
        run_thickness_sweep(config)
        
    # else:
    #     # --- EXISTING SINGLE RUN LOGIC ---
    #     res = run_heterostructure_simulation(config)
    
    #     z, Vtot_c, Vtot_v = res["z"], res["Vtot_c"], res["Vtot_v"]
    #     V0_c, V0_v, Ef = res["V0_c"], res["V0_v"], res["Ef"]
    #     Ec, Ev, psic, psiv = res["Ec"], res["Ev"], res["psic"], res["psiv"]
    #     N_3D, P_3D, rho_net = res["N_3D"], res["P_3D"], res["rho_net"]
    #     ND_plus_z, NA_minus_z = res["ND_plus_z"], res["NA_minus_z"]
    #     flat_layers = res["flat_layers"]
    
    #     ScF = 0.05 # Scaling factor for visualization heights
        
    #     # ... [Keep all of your existing Plot 1 through Plot 6 logic here] ...
        
    #     plt.show()

    # =============================================================
    # SET 1: ORIGINAL PLOTS (WITH COLORS & HEATMAPS)
    # =============================================================

    # --- PLOT 1: Original Band Diagram + Wavefunctions ---
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

    # --- PLOT 2: Original Band Diagram + Total Density Heatmap ---
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

    # --- PLOT 3: Original Charge Distribution ---
    fig3, ax3 = plt.subplots(figsize=(10, 5))
    add_layer_background(ax3, flat_layers, show_colors=True)

    ax3.plot(z*1e9, -N_3D, lw=2, color='blue', label=r'$-N_{3D}$ (Electrons)')
    ax3.plot(z*1e9, P_3D, lw=2, color='red', label=r'$P_{3D}$ (Holes)')
    ax3.plot(z*1e9, rho_net, lw=2, color='black', alpha=0.8, label=r'$\rho_{net}$')

    fixed_pos = ND_plus_z + res["N_D_dipole"]
    fixed_neg = NA_minus_z + res["N_A_dipole"]
    
    pos_idx = np.where(fixed_pos > 1e10)[0] 
    neg_idx = np.where(fixed_neg > 1e10)[0]

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


    # =============================================================
    # SET 2: CLEAN PLOTS (NO COLORS, 1D DENSITY INSTEAD OF HEATMAP)
    # =============================================================

############################################################################################################################


    # --- PLOT 4: Clean Conduction Band + Electron Wavefunctions ---
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


    # # --- PLOT 4: Clean Band Diagram + Wavefunctions ---
    # fig4, ax4 = plt.subplots(figsize=(10, 6))
    # add_layer_background(ax4, flat_layers, show_colors=False)

    # ax4.plot(z*1e9, V0_c, 'b--', alpha=0.5, lw=2, label='Ec (Flat)', zorder=10)
    # ax4.plot(z*1e9, V0_v, 'r--', alpha=0.5, lw=2, label='Ev (Flat)', zorder=10)
    # ax4.plot(z*1e9, Vtot_c, 'b', lw=2, label='Ec (Bent)', zorder=10)
    # ax4.plot(z*1e9, Vtot_v, 'r', lw=2, label='Ev (Bent)', zorder=10)
    # ax4.axhline(Ef, color='lime', ls='--', lw=2, label=r'$E_f$', zorder=10)

    # if len(Ec) > 0:
    #     for i in range(min(config["n"], len(Ec))):
    #         intensity = np.abs(psic[:, i])**2
    #         if np.max(intensity) > 1e-15:
    #             norm_intensity = intensity / np.max(intensity)
    #             ax4.plot(z*1e9, np.full(len(z), Ec[i]) + (ScF * norm_intensity), 
    #                      color='blue', linestyle='-', lw=1.5, alpha=0.8, zorder=5)
    #             ax4.axhline(Ec[i], color='blue', linestyle='-', lw=0.5, alpha=0.3, zorder=4)

    # if len(Ev) > 0:
    #     for i in range(min(config["n"], len(Ev))):
    #         intensity = np.abs(psiv[:, i])**2
    #         if np.max(intensity) > 1e-15:
    #             norm_intensity = intensity / np.max(intensity)
    #             ax4.plot(z*1e9, np.full(len(z), Ev[i]) - (ScF * norm_intensity), 
    #                      color='red', linestyle='-', lw=1.5, alpha=0.8, zorder=5)
    #             ax4.axhline(Ev[i], color='red', linestyle='-', lw=0.5, alpha=0.3, zorder=4)

    # # ax4.set(xlim=(0, z[-1]*1e9), ylim=(np.min(Vtot_v) - 0.2, np.max(Vtot_c) + 0.2),
    # #         xlabel='z (nm)', ylabel='Energy (eV)', title='[Clean] Band Diagram & Individual Wavefunctions')
    # # Fixed viewing window
    # ax4.set(
    #     xlim=(20, 40),
    #     ylim=(-4, 0),
    #     xlabel='z (nm)',
    #     ylabel='Energy (eV)',
    #     title='[Clean] Band Diagram & Individual Wavefunctions'
    # )

    # ax4.grid(alpha=0.3)
    # ax4.legend(loc='center right')
    # fig4.tight_layout()

############################################################################################################################

    # --- PLOT 5: Clean Band Diagram + 1D Carrier Density (Replaces Heatmap) ---
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

    # --- PLOT 6: Clean Charge Distribution ---
    fig6, ax6 = plt.subplots(figsize=(10, 5))
    add_layer_background(ax6, flat_layers, show_colors=False)

    ax6.plot(z*1e9, -N_3D, lw=2, color='blue', label=r'$-N_{3D}$ (Electrons)')
    ax6.plot(z*1e9, P_3D, lw=2, color='red', label=r'$P_{3D}$ (Holes)')
    ax6.plot(z*1e9, rho_net, lw=2, color='black', alpha=0.8, label=r'$\rho_{net}$')

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

# =============================================================
    # --- PLOT 7: Electric Field and Dynamic Dielectric Constant ---
    # =============================================================
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

    
    plt.show()