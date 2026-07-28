import numpy as np
import dolfin as df
import warnings
import matplotlib.pyplot as plt

# Import from your local files
from master_helper import find_Ef_general, Schroed1D_dolfin, Poisson_dolfin
# from config_LAOiSTO import config
from config_STO100iSTO111 import config

# from config_LLOiSBO import config

# from config_STOiSAP import config

warnings.filterwarnings("ignore")
df.set_log_level(df.LogLevel.WARNING)

def run_heterostructure_simulation(config):
    e = 1.602176487E-19 #e=1.602176487×10^−19 C
    
    # --- 1. Expand Stack Definition ---
    flat_layers = [] # will contain every atomic layer lik3 ababababcdcdcdcdcd
    MATERIALS = config["MATERIALS"]
    for item in config["stack_definition"]:
        mat_name, count = item.split('*')
        count = int(count)
        mat = MATERIALS[mat_name]
        for _ in range(count):
            for layer in mat["layers"]: # adds the below extra information to all the layers
                layer_copy = layer.copy()
                layer_copy["mat_name"] = mat_name
                layer_copy["uc_area"] = mat["uc_area"]
                flat_layers.append(layer_copy)

    # --- 2. Mesh Generation ---
    min_t = min(L["t_layer"] for L in flat_layers) # finds the thinnest layer 
    dz_fine = min_t / config["points_per_layer"] # thinnest layer divided by points giving spacing of the meshing
    z_max = sum(L["t_layer"] for L in flat_layers) # gives total thickness
    
    Nz = int(round(z_max / dz_fine)) + 1
    z = np.linspace(0, z_max, Nz)
    dz_node = z[1] - z[0]

    print("dz_fine, dz_node = ", dz_fine, dz_node)
    print("dz_fine - dz_node = ", dz_fine - dz_node)

    # --- 3. Property Array Initialization ---
    V0_c, V0_v = np.zeros(Nz), np.zeros(Nz)
    meff_e_arr, meff_h_arr = np.zeros(Nz), np.zeros(Nz)
    Epsi_arr = np.zeros(Nz)
    N_D_array, N_A_array = np.zeros(Nz), np.zeros(Nz)
    N_D_dipole, N_A_dipole = np.zeros(Nz), np.zeros(Nz)
    
    layer_boundaries = [0.0]
    curr_z = 0.0
    for L in flat_layers:
        curr_z += L["t_layer"]
        layer_boundaries.append(curr_z)

    layer_idx_map = np.searchsorted(layer_boundaries, z, side='right') - 1
# Now every mesh point knows which physical layer it belongs to.
    layer_idx_map[layer_idx_map >= len(flat_layers)] = len(flat_layers) - 1

    # Find the primary interface node for tracking
    interface_idx = 0

# Fill material properties
    for i in range(Nz - 1):
        if flat_layers[layer_idx_map[i]]["mat_name"] != flat_layers[layer_idx_map[i+1]]["mat_name"]:
# "Does the material change between these two neighbouring nodes?"
            interface_idx = i
            break

    for k in range(Nz):
        L = flat_layers[layer_idx_map[k]]
        V0_c[k] = L["Ec"]
        V0_v[k] = L["Ev"]
        meff_e_arr[k] = L["meff_e"]
        meff_h_arr[k] = L["meff_h"]
        if L.get("is_eps_dynamic", False):
            Epsi_arr[k] = 1.0 / L["eps_A"] 
        else:
            Epsi_arr[k] = L["eps_r"]

    # --- 4. Intrinsic Layer-Based Charges (Polar Planes) ---
    for idx, L in enumerate(flat_layers):
        if "layer_charge" in L and L["layer_charge"] != 0.0:
            # Convert 2D formal charge (e per uc_area) to 3D volume density (e/m^3)
            val_2d = L["layer_charge"] / L["uc_area"]
            val_3d = val_2d / L["t_layer"]
            
            z_start = layer_boundaries[idx]
            z_end = layer_boundaries[idx + 1]
            
            # Mask assigns the charge strictly within this atomic plane's bounds
            mask = (z >= z_start) & (z <= z_end)
            
            if L["layer_charge"] > 0:
                N_D_array[mask] += val_3d
            else:
                N_A_array[mask] += abs(val_3d)

#  After this section, every mesh node has a complete set of material properties
# for example 
#  Mesh node
# Material, Ec, Ev, me*, mh*, ε

    # Helper function for interface-relative targeting
    #  This helper function is designed to find specific atomic layers at 
    # specified distance from a material interface. 
    # Instead of using absolute layer numbers (e.g., layer 153), 
    # it lets you say things like:
    # "Find the first 5 TiO₂ layers in the STO film after the LLO/STO interface."
    # "Find the first 3 LaO layers inside LLO before the interface."


    def get_interface_relative_indices(
        flat_layers,          # Complete expanded stack
        interface_pair,       # (mat_left, mat_right)
        direction,            # "left" or "right" of interface
        target_layer,         # Atomic layer to count
        uc_start,             # First occurrence to include
        uc_end=-1             # Last occurrence to include (-1 = all)
    ):
        mat_left, mat_right = interface_pair

        idx = -1

        for i in range(len(flat_layers) - 1):
            if flat_layers[i]["mat_name"] == mat_left and flat_layers[i+1]["mat_name"] == mat_right:
                idx = i
                break
        if idx == -1: return []     

        indices = []
        count = 0
        if direction == 'right': 
            for i in range(idx + 1, len(flat_layers)):
                if flat_layers[i]["name"] == target_layer:
                    count += 1
                    if uc_start <= count <= uc_end: indices.append(i)
                    if count >= uc_end: break
        elif direction == 'left': 
            for i in range(idx, -1, -1):
                if flat_layers[i]["name"] == target_layer:
                    count += 1
                    if uc_start <= count <= uc_end: indices.append(i)
                    if count >= uc_end: break
        return indices
    
    # --- 5. Interface-Relative Doping Injection ---
    if "doping_profiles" in config:
        for dop in config["doping_profiles"]:
            target_arr = N_D_array if dop["type"] == "donor" else N_A_array # the units here are m-3
            # Chooses where to place the dopants.
            val = dop["value"] # untis again are m-3 (val number of donors per cubic meter)
            
            target_indices = get_interface_relative_indices(
                flat_layers, dop["interface"], dop["direction"], 
                dop["target_layer"], dop["uc_range"][0], dop["uc_range"][1]
            )
            
            for idx in target_indices:
                z_start_layer = layer_boundaries[idx]
                z_end_layer = layer_boundaries[idx + 1]
                
                if dop["profile"] == "square":
                    mask = (z >= z_start_layer) & (z <= z_end_layer)
                    target_arr[mask] += val
                elif dop["profile"] == "delta":
                    z_center = (z_start_layer + z_end_layer) / 2.0
                    idx_center = np.argmin(np.abs(z - z_center))
                    target_arr[idx_center] += val / dz_node

    # --- 6. Split Charge / Dipole Polarization ---
    P_layers = []
    if "dipole_configs" in config:
        for dip in config["dipole_configs"]:
            target_indices = get_interface_relative_indices(
                flat_layers, dip["interface"], dip["direction"], 
                dip["target_layer"], dip["uc_range"][0], dip["uc_range"][1]
            )
            
            m_charge = dip["split_m_charge"]
            
            for step, idx in enumerate(target_indices):
                L = flat_layers[idx]
                decay = np.exp(-step / dip["decay_rate"]) if dip["decay_mode"] == 'exponential' else 1.0
                
                d_plus = dip["delta_plus_0"] * (decay if dip["fixed_delta"] != 'plus' else 1.0)
                d_minus = dip["delta_minus_0"] * (decay if dip["fixed_delta"] != 'minus' else 1.0)
                
                if (d_plus + d_minus) > L["t_layer"]:
                    scale = (L["t_layer"] * 0.99) / (d_plus + d_minus)
                    d_plus *= scale; d_minus *= scale
                
                z_center = layer_boundaries[idx] + (L["t_layer"] / 2.0)
                z_pos = z_center - d_plus
                z_neg = z_center + d_minus
                
                sigma_0 = 1.0 / L["uc_area"]
                idx_pos = np.argmin(np.abs(z - z_pos))
                idx_neg = np.argmin(np.abs(z - z_neg))
                
                N_D_dipole[idx_pos] += (m_charge * sigma_0) / dz_node
                N_A_dipole[idx_neg] += (m_charge * sigma_0) / dz_node
                
                dipole_moment = m_charge * e * (z_pos - z_neg)
                P_layers.append(dipole_moment / (L["uc_area"] * L["t_layer"]))

    # --- 7. FEniCS Initialization ---
    mesh = df.IntervalMesh(Nz - 1, z[0], z[-1])
    mesh.coordinates()[:, 0] = z
    V_space = df.FunctionSpace(mesh, 'CG', 1)
    
    alpha = config["alpha"]
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

            if len(Ec) == 0: break

            Ef, N_n, P_m, ND_plus_z, NA_minus_z, net_integral_charge = find_Ef_general(
                Ec, Ev, psic, psiv, N_D_array, N_A_array, Vtot_c, Vtot_v, z, config["T"],
                meff_e_arr, meff_h_arr, config["E_bind_D"], config["E_bind_A"], config["fully_ionized_flag"]
            )

            N_3D = np.sum(np.tile(N_n, (len(z), 1)) * np.abs(psic)**2, axis=1) if len(Ec) > 0 else np.zeros_like(z)
            P_3D = np.sum(np.tile(P_m, (len(z), 1)) * np.abs(psiv)**2, axis=1) if len(Ev) > 0 else np.zeros_like(z)

            rho_net = P_3D - N_3D + ND_plus_z - NA_minus_z + N_D_dipole - N_A_dipole

            Vout = Poisson_dolfin(
                V_space, rho_net, Epsi_arr, 
                config.get("bc_left_type", "neumann"), config.get("V_left_applied", 0.0),
                config.get("bc_right_type", "neumann"), config.get("V_right_applied", 0.0)
            )

            E_cm = np.abs(np.gradient(Vout, z)) / 100.0
            Epsi_new = np.zeros_like(Epsi_arr)
            for k in range(Nz):
                L = flat_layers[layer_idx_map[k]]
                if L.get("is_eps_dynamic", False):
                    Epsi_new[k] = 1 / (L["eps_A"] + L["eps_B"] * E_cm[k])
                else:
                    Epsi_new[k] = L["eps_r"]
                    
            Epsi_arr = alpha * Epsi_new + (1 - alpha) * Epsi_arr
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
            P_hole =  e * np.trapz((z - z0) * P_3D, x=z) / L_z
            P_dipole_inserted = sum(P_layers)

            print(f"==========================================================================================")
            print(f" LOOP {nloop:03d} | Error: {Err:.3e} (Tol: {config['tol']:.1e}) | Ef: {Ef:.5f} eV")
            print(f"------------------------------------------------------------------------------------------")
            print(f" 2D Electron Dens (np.sum N_n) : {-Q_e:.4e} m^-2")
            print(f" 2D Hole Dens (np.sum P_m)     : {+Q_h:.4e} m^-2")
            print(f" 2D Donor Dens (Trapz ND+)     : {+Q_ND:.4e} m^-2")
            print(f" 2D Acceptor Dens (Trapz NA-)  : {-Q_NA:.4e} m^-2")
            print(f" Net Charge to Poisson (Trapz) : {Q_net:.4e} m^-2")
            print(f" Q Residual (Convergence check): {net_integral_charge:.4e} m^-2")
            print(f"------------------------------------------------------------------------------------------")
            print(f" Charge Dens @ Interface Node  : {rho_net[interface_idx]:.4e} m^-3 (z = {z[interface_idx]*1e9:.2f} nm)")
            print(f"==========================================================================================\n")

            if Err < config["tol"]:
                convergence_met = True
                break
                            
            if prev_Err is not None:
                alpha = max(alpha * config["shrink_factor"], config["alpha_min"]) if Err > prev_Err else min(alpha * config["grow_factor"], config["alpha_max"])
            prev_Err = Err
            Vin = alpha * Vout + (1 - alpha) * Vin

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
        "rho_net": rho_net, "flat_layers": flat_layers
    }

def add_layer_background(ax, flat_layers):
    colors = {
        "TiO2": "white",         
        "SrO": "lightgray",     
        "LaO": "#bcdcff",       
        "AlO2": "#ffcccb"
    }

    curr_z = 0.0
    prev_mat = None

    for L in flat_layers:
        t = L["t_layer"] * 1e9
        mat_name = L["mat_name"]
        
        c = colors.get(L["name"], "whitesmoke")
        ax.axvspan(curr_z, curr_z + t, facecolor=c, alpha=0.30, edgecolor='none', zorder=-100)

        if prev_mat is not None and mat_name != prev_mat:
            ax.axvline(curr_z, color='grey', lw=2, zorder=-50)

        curr_z += t
        prev_mat = mat_name

if __name__ == "__main__":
    res = run_heterostructure_simulation(config)

    z, Vtot_c, Vtot_v = res["z"], res["Vtot_c"], res["Vtot_v"]
    V0_c, V0_v, Ef = res["V0_c"], res["V0_v"], res["Ef"]
    Ec, Ev, psic, psiv = res["Ec"], res["Ev"], res["psic"], res["psiv"]
    N_3D, P_3D, rho_net = res["N_3D"], res["P_3D"], res["rho_net"]
    ND_plus_z, NA_minus_z = res["ND_plus_z"], res["NA_minus_z"]
    flat_layers = res["flat_layers"]

    if len(Ec) > 0:
        fig, ax = plt.subplots(figsize=(10, 6))
        add_layer_background(ax, flat_layers)

 # 1. Bring the lines to the FRONT by adding zorder=10
        ax.plot(z*1e9, V0_c, 'b--', alpha=0.5, lw=2, label='Ec (Flat)', zorder=10)
        ax.plot(z*1e9, V0_v, 'r--', alpha=0.5, lw=2, label='Ev (Flat)', zorder=10)
        ax.plot(z*1e9, Vtot_c, 'b', lw=2, label='Ec (Bent)', zorder=10)
        ax.plot(z*1e9, Vtot_v, 'r', lw=2, label='Ev (Bent)', zorder=10)
        ax.axhline(Ef, color='lime', ls='--', lw=2, label=r'$E_f$', zorder=10)

        ScF = 0.05
        for i in range(min(config["n"], len(Ec))):
            intensity = np.abs(psic[:, i])**2
            norm_intensity = intensity / np.max(intensity)
            mask = norm_intensity > 1e-3
            
            # 2. Send the heatmap to the BACK by changing zorder to 2
            ax.scatter(z[mask]*1e9, np.full(np.sum(mask), Ec[i]) + (ScF * norm_intensity[mask]), 
                       c=norm_intensity[mask], cmap='Blues', s=15, alpha=0.8, edgecolor='none', zorder=2)



# Plot Hole Wavefunctions (Reds) - ADD THIS BACK
        if len(Ev) > 0:
            for i in range(min(config["n"], len(Ev))):
                intensity = np.abs(psiv[:, i])**2
                norm_intensity = intensity / np.max(intensity)
                mask = norm_intensity > 1e-3
                # Only plot if there's actually a non-zero probability of finding a hole
                if np.max(intensity) > 1e-15: 
                    ax.scatter(z[mask]*1e9, np.full(np.sum(mask), Ev[i]) - (ScF * norm_intensity[mask]), 
                               c=norm_intensity[mask], cmap='Reds', s=15, alpha=0.8, edgecolor='none', zorder=2)


                    
        ax.set(xlim=(0, z[-1]*1e9), ylim=(np.min(Vtot_v) - 0.2, np.max(Vtot_c) + 0.2),
               xlabel='z (nm)', ylabel='Energy (eV)', title='Heterostructure Band Diagram')
        ax.grid(alpha=0.3)
        ax.legend(loc='center right')
        plt.tight_layout()

        # Charge Distribution
        fig, ax = plt.subplots(figsize=(10, 5))
        add_layer_background(ax, flat_layers)

        ax.plot(z*1e9, -N_3D, lw=2, color='blue', label=r'$-N_{3D}$ (Electrons)')
        ax.plot(z*1e9, P_3D, lw=2, color='red', label=r'$P_{3D}$ (Holes)')
        ax.plot(z*1e9, rho_net, lw=2, color='black', alpha=0.8, label=r'$\rho_{net}$')

        pos_idx = np.where(ND_plus_z > 0)[0]
        neg_idx = np.where(NA_minus_z > 0)[0]
        if len(pos_idx) > 0:
            ax.vlines(z[pos_idx]*1e9, 0, ND_plus_z[pos_idx], color='red', linestyle=':', lw=2, label=r'Polar $+$')
        if len(neg_idx) > 0:
            ax.vlines(z[neg_idx]*1e9, 0, -NA_minus_z[neg_idx], color='blue', linestyle=':', lw=2, label=r'Polar $-$')

        ax.set_yscale('symlog', linthresh=1e25) 
        ax.axhline(0, color='gray', ls='--', lw=1)
        ax.set(xlabel='z (nm)', ylabel=r'Charge Density (m$^{-3}$)', title='Charge Distribution Profile (Logarithmic)')
        ax.grid(alpha=0.3)
        ax.legend(loc='upper right', bbox_to_anchor=(1.25, 1.0))
        plt.tight_layout()
        plt.show()