# OxideMemSP_edit1

Modular version of `OxideMemSP_backup`. The physics and the plots are the same; the code is split into modules, and each config can choose its mixing scheme.

## Run

1. Pick a config by uncommenting its import in `master_run.py`.
2. Set `"mixing_method"` in that config to `"anderson"` or `"alpha"`. To override it without editing the config, set `MIXING_METHOD_OVERRIDE` in `master_run.py`.
3. Run `python master_run.py`.

## Mixing settings (in every config)

| `mixing_method` | Keys used | Behaviour |
|---|---|---|
| `"anderson"` | `mix_m`, `mix_beta`, `mix_max_step`, `mix_beta_eps`, `eps_start`, `eps_tol` | Anderson mixing of V. The field-dependent ε gets its own fixed damping, which starts after `eps_start` loops. Convergence needs `Err < tol` and a relative ε change below `eps_tol`. |
| `"alpha"` | `alpha`, `alpha_min`, `alpha_max`, `grow_factor`, `shrink_factor` | Adaptive linear mixing of V and ε. α shrinks when the error grows and grows when the error falls. Convergence needs `Err < tol`. |

## Field-dependent permittivity ε(E)

Layers with `"is_eps_dynamic": True` and an `"eps_dynamic_func"` get ε from the local |E| in every loop.

- **The model has a ceiling on D.** The STO model ε(E)·ε₀·E can never exceed about 0.65 C/m², which is about 0.62 e/uc. Any fixed sheet charge larger than that inside a field-dependent layer has no self-consistent solution: ε falls towards 0 and the loop diverges. Before the loop starts, the solver now prints an `[eps(E) WARNING]` for every charge like this.
- **Dipole layers keep a constant ε.** By default, a dipole layer uses its `eps_r` even if the layer is field dependent. To let ε(E) act inside one, set `"eps_dynamic_inside": True` on that dipole entry.
- **`"eps_update"` sets how ε is updated.**
  - `"damped"` (original): ε comes from the previous output field, damped by the mixer.
  - `"input"` (the default in the configs): ε is recomputed from the input potential each loop after `eps_start`. Anderson converges much more reliably this way. In tests, STOiGTO converged in 87 loops instead of never.
- **Convergence waits for ε.** If any layer is field dependent, the loop cannot declare convergence before `eps_start`. Before this change, Anderson could stop early with ε never updated.
- **ε(E) must be continuous.** The fit has to join smoothly at 0.27 MV/cm. With `C = 736.30` it does; with `C = 719.307` it jumps from 300 to 237, and the loop oscillates.

## Fixed-charge profiles (in every config)

Every fixed charge (polar plane, doping, dipole) has a shape inside its layer:

| Shape | Polar plane / doping | Dipole |
|---|---|---|
| `"delta"` (default, original behaviour) | whole charge on the layer-centre node (atomic plane) | +q on `centre − delta_plus_0`, −q on `centre + delta_minus_0` |
| `"step"` | uniform over the whole layer | +q uniform over the left half, −q over the right half |
| `"gaussian"` | Gaussian on the layer centre, truncated to the layer | Gaussians on the + and − points, truncated to the layer |

- The defaults per kind are set in `"charge_profile_defaults"`. Gaussian width is `sigma = sigma_frac × t_layer`.
- To override a single doping or dipole entry, add `"profile"` to it, optionally with `"sigma_frac"` or `"sigma"` (in m).
- To override a polar plane, add `"charge_profile"` to that layer's dict in `MATERIALS`, optionally with `"sigma_frac"` or `"sigma"`.
- Charge is spread by integrating the shape over each node's owned width (`charge_profiles.py`). The total sheet charge is therefore exact, and no charge leaves its layer.
- The dipole moment depends on the shape. Step gives a separation of `t_layer/2`, and gaussian gives one slightly smaller than the delta separation because of the truncation.

`charge_profiles_example.png` compares the three shapes at the STOiSAP interface.

## Printed diagnostics

Two switches, in every config, control the extra output:

| Switch | When it prints | What it shows |
|---|---|---|
| `"print_device_check"` | once, before the loop starts | Stack and mesh: layer count, thickness, elements per layer, whether every boundary, centre and dipole point is a mesh node, and whether every node carries its own layer's properties. Charges: every polar plane, doping and dipole in e/uc of its host layer, with its sheet density, the charge actually placed, and whether it stays inside its layer. Arrays: whether the charge arrays equal the sum of the individual charges. Totals: the net fixed charge, the electrons that charge neutrality needs, and any suspicious `uc_area`. |
| `"print_polarization"` | every loop | One line per inserted dipole: offset, q, separation, dipole moment and P in µC/cm². Then the polarization balance. |

**The polarization balance.** Every term is a dipole moment per area, printed in e·nm/uc. The origin is the end of the stack where the potential is fixed (the Dirichlet end).

```
inserted + free charges (ions/doping + electrons/holes) + lattice + eps0*dU = 0
```

- **Lattice** is ∫ ε₀(ε−1)E dz, the dielectric's own polarization.
- **ε₀·ΔU** is what is left over, i.e. the potential step across the whole stack.
- **"Screening"** is net − inserted, which equals free charges + lattice.
- The sum check is evaluated with the same piecewise-linear integrals FEniCS assembles, so it should print about 1e-12.

## Choosing plots (in every config, `"plots"`)

Set each plot to `True` (draw it) or `False` (skip it):

Each key matches the figure title:

| Key | Figure title |
|---|---|
| `band_diagram_wavefunctions` | Band Diagram and Wavefunctions (**off by default**) |
| `carrier_distribution_heatmap` | Carrier Distribution Heatmap |
| `charge_distribution_shaded` | Charge Distribution with Layer Shading |
| `conduction_band_wavefunctions` | Conduction Band and Electron Wavefunctions (fixed window) |
| `band_diagram_carrier_density` | Band Diagram and Carrier Density (log right axis) |
| `charge_distribution` | Charge Distribution (plain background) |
| `electric_field_dielectric_constant` | Electric Field and Dielectric Constant |

The earlier key names (`carrier_heatmap`, `cb_wavefunctions`, …) still work. Every line in every plot has the same width (`plot_settings["linewidth"]`, default 2.0; only the thin layer-centre lines of the shading stay thin). Every axis has minor ticks: evenly spaced on linear axes, 2–9 within each decade on log axes.

**The two charge distribution plots have two y-axes:**
- **Left axis (linear, signed):** fixed charges (polar planes, doping, dipoles; + red / − blue) and ρ_net.
- **Right axis (log):** N_3D and P_3D.
- `charge_xlim = (zmin, zmax)` zooms in, in nm. `charge_linear_ylim = (ymin, ymax)` fixes the left axis in m⁻³, for example to see a small doping next to tall dipole spikes.

## Plot settings (in every config, `"plot_settings"`)

- **What the density plots show.** Every density plot is a 3D volume density in m⁻³ (N_3D, P_3D, ρ_net, fixed charge), not a sheet density. Carrier densities use a **log** axis; the charge distribution plots add a linear signed axis (see above).
  - Signed quantities are shown as magnitudes. |ρ_net| is solid where ρ > 0 and dashed where ρ < 0.
  - Fixed + charges are dark red and fixed − charges are dark blue.
  - The lower limit is `density_floor` in m⁻³, or `top × density_decades` (10 decades by default) when `density_floor` is None.
- **Band-plot energy window** (plots 1, 2 and 5). It spans the bent and flat bands of the material with the **largest real band gap**, plus Ef.
  - Gaps above `max_real_gap_eV` (20 eV) count as vacuum and are ignored, as are materials listed in `band_ylim_exclude`.
  - `band_ylim = (ymin, ymax)` fixes the window by hand.
- **Layer shading.** `material_colors` gives each material one color (any matplotlib color). `background_alpha` sets the shading opacity, from 0 (invisible) to 1 (solid).

## Electrodes as boundary conditions (`config_W_YHO_LSMO_STO.py`)

**How the electrodes are modelled.**
- Metal electrodes don't need their bands solved. Each one is a fixed potential (a Dirichlet boundary) behind a thin, charge-free **screening gap** of thickness λ (the Thomas–Fermi screening length) and permittivity ε_e.
- Outside the metal this is electrostatically equivalent to a Thomas–Fermi electrode (Zhuravlev *et al.*, PRL 94, 246802 (2005)).
- Without the gap (ideal electrodes), a uniform polarization would be screened completely, and the barriers wouldn't depend on P.

**Electrodes inside the stack.** The STO substrate is kept below the LSMO electrode, so LSMO is a fixed-potential region inside the stack rather than an end boundary.
- Mark such layers with `"is_electrode": True` and `"electrode_potential": U`. The potential is then held at U on every node inside them, via `Poisson_dolfin(..., fixed_regions=...)`.
- This is the only change in `master_helper.py`. The previous file is kept as `master_helper_original.py`, and with no fixed regions the new one behaves exactly like it.

**Fermi level.** `"fermi_level": -Phi_LSMO` fixes Ef from the electrodes. Without it, the code uses global charge neutrality as before.

**Electrode potentials.** U_left = 0 and U_right = Φ_W − Φ_LSMO − V_bias. With no polarization this gives the work-function (Schottky–Mott) barriers Φ_M − χ.

**Ferroelectric.** Every YHO cell carries the same fixed dipole: one `dipole_configs` entry with offsets 1..N, a `"step"` profile and q = 2 × 0.254 e. That gives P = 36 μC/cm² in every cell.
- `POL = +1` points P toward LSMO, `−1` toward W, and `0` turns it off.

**Printed diagnostics.** The polarization balance gains an "electrode charges" term and prints the screening charge on each electrode. Identical dipoles of one entry are printed as a single line.

## Files

| File | Contents |
|---|---|
| `master_run.py` | Entry point. Calls simulation, then sweep, then plots. |
| `master_helper.py` | `find_Ef_general`, `Schroed1D_dolfin`, `Poisson_dolfin` (+ optional `fixed_regions` for electrodes inside the stack) |
| `master_helper_original.py` | The helper before that change |
| `config_*.py` | Materials, stack, BCs, doping/dipoles, sweep, solver and mixing settings |
| `stack_geometry.py` | Stack expansion, layer boundaries, interface-relative layer lookup |
| `mesh_generation.py` | Exact points, uniform mesh, tolerance-based merging, nodal widths |
| `material_properties.py` | Layer index map, interface node, band edges, masses, ε(E) |
| `fixed_charges.py` | Polar-plane charges, doping profiles, split-charge dipoles |
| `charge_profiles.py` | Delta / step / gaussian weights for distributing a sheet charge over nodes |
| `mixers.py` | `AndersonMixer`, `AndersonMixingScheme`, `AlphaMixingScheme` |
| `scf_solver.py` | FEniCS initialization and the self-consistent loop |
| `simulation.py` | `build_device` and `run_heterostructure_simulation` |
| `thickness_sweep.py` | Critical-thickness sweep |
| `diagnostics.py` | Device check at start-up, dipole lines and polarization balance every loop |
| `plotting.py` | `add_layer_background`, plots 1–7, sweep plot |
