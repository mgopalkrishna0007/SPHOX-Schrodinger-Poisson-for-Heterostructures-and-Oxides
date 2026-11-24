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

## Files

| File | Contents |
|---|---|
| `master_run.py` | Entry point. Calls simulation, then sweep, then plots. |
| `master_helper.py` | `find_Ef_general`, `Schroed1D_dolfin`, `Poisson_dolfin` (unchanged) |
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
| `plotting.py` | `add_layer_background`, plots 1–7, sweep plot |
