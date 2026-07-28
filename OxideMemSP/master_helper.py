import numpy as np
import scipy.integrate as integrate
import scipy.optimize as optimize
import matplotlib.pyplot as plt
import dolfin as df
import warnings

warnings.filterwarnings("ignore")
df.set_log_level(df.LogLevel.WARNING)

# =============================================================================
# Helper Function: find_Ef_general
# =============================================================================
def find_Ef_general(Ec, Ev, psic, psiv, N_D_array, N_A_array, Vc_local, Vv_local, z, T, meff_e_array, meff_h_array, E_bind_D=0.0, E_bind_A=0.0, fully_ionized=True):
    e = 1.602176487E-19
    kB = 1.3806488E-23
    m0 = 9.10938188E-31
    hbar = 6.62606896E-34 / (2 * np.pi)
    kT_eV = (kB * T) / e if T > 0 else 1e-10

    m_eff_n = np.array([np.trapz(meff_e_array * np.abs(psic[:, i])**2, x=z) for i in range(len(Ec))])
    m_eff_m = np.array([np.trapz(meff_h_array * np.abs(psiv[:, i])**2, x=z) for i in range(len(Ev))])

    twoD_DOS_e = (m_eff_n * m0 * kB * T) / (np.pi * hbar**2)
    twoD_DOS_h = (m_eff_m * m0 * kB * T) / (np.pi * hbar**2)

    def Q_residual(Ef):
        N_n = twoD_DOS_e * np.logaddexp(0, (Ef - Ec) / kT_eV)
        P_m = twoD_DOS_h * np.logaddexp(0, (Ev - Ef) / kT_eV)
        if fully_ionized:
            ND_plus = N_D_array
            NA_minus = N_A_array
        else:
            ND_plus = N_D_array / (1 + 2 * np.exp((Ef - (Vc_local - E_bind_D)) / kT_eV))
            NA_minus = N_A_array / (1 + 0.5 * np.exp((Vv_local + E_bind_A - Ef) / kT_eV))
        return np.sum(P_m) - np.sum(N_n) + np.trapz(ND_plus, x=z) - np.trapz(NA_minus, x=z)

    min_Ef = np.min(Vv_local) - 1.0
    max_Ef = np.max(Vc_local) + 1.0

    try:
        Ef_sol = optimize.brentq(Q_residual, min_Ef, max_Ef)
    except ValueError:
        print("doing bisection")
        Ef_sol = optimize.bisect(Q_residual, min_Ef - 5.0, max_Ef + 5.0)

    N_n = twoD_DOS_e * np.logaddexp(0, (Ef_sol - Ec) / kT_eV)
    P_m = twoD_DOS_h * np.logaddexp(0, (Ev - Ef_sol) / kT_eV)

    if fully_ionized:
        ND_plus_z = N_D_array
        NA_minus_z = N_A_array
    else:
        ND_plus_z = N_D_array / (1 + 2 * np.exp((Ef_sol - (Vc_local - E_bind_D)) / kT_eV))
        NA_minus_z = N_A_array / (1 + 0.5 * np.exp((Vv_local + E_bind_A - Ef_sol) / kT_eV))

    net_integral_charge = Q_residual(Ef_sol)
    return Ef_sol, N_n, P_m, ND_plus_z, NA_minus_z, net_integral_charge

# =============================================================================
# FEniCS Solver: Schrodinger Equation
# =============================================================================
def Schroed1D_dolfin(V_space, Vtot_array, Mass_array, n_eig):
    mesh = V_space.mesh()
    u = df.TrialFunction(V_space)
    v = df.TestFunction(V_space)

    h = 6.62606896E-34
    hbar = h / (2 * np.pi)
    e = 1.602176487E-19
    m0 = 9.10938188E-31

    v2d = df.vertex_to_dof_map(V_space)
    Vtot_func = df.Function(V_space)
    Vtot_func.vector()[v2d] = Vtot_array

    Mass_func = df.Function(V_space)
    Mass_func.vector()[v2d] = Mass_array

    scale = df.Constant(1e9)
    coeff_val = hbar**2 / (2 * m0 * e)

    a = scale * df.Constant(coeff_val) * (1.0 / Mass_func) * df.inner(df.grad(u), df.grad(v)) * df.dx + scale * Vtot_func * u * v * df.dx
    m_form = scale * u * v * df.dx

    A = df.PETScMatrix()
    M_mat = df.PETScMatrix()
    df.assemble(a, tensor=A)
    df.assemble(m_form, tensor=M_mat)

    solver = df.SLEPcEigenSolver(A, M_mat)
    solver.parameters['problem_type'] = 'gen_hermitian'
    solver.parameters['spectral_transform'] = 'shift-and-invert'
    solver.parameters['spectral_shift'] = np.min(Vtot_array) - 0.01

    solver.solve(n_eig)

    nconv = solver.get_number_converged()
    if nconv < n_eig:
        n_eig = nconv

    E = []
    psi_array = np.zeros((len(Vtot_array), n_eig))

    for i in range(n_eig):
        r, _, rx, _ = solver.get_eigenpair(i)
        E.append(r)

        eig_func = df.Function(V_space)
        eig_func.vector()[:] = rx
        psi_val = eig_func.vector().get_local()[v2d]

        z_array = mesh.coordinates()[:, 0]
        norm = np.sqrt(np.trapz(np.abs(psi_val)**2, x=z_array))
        psi_array[:, i] = psi_val / norm

    if len(E) > 0:
        idx = np.argsort(E)
        E = np.array(E)[idx]
        psi_array = psi_array[:, idx]

    return E, psi_array

# =============================================================================
# FEniCS Solver: Poisson Equation
# =============================================================================
# def Poisson_dolfin(V_space, rho_net_array, Epsi_array, V_left=0.0, V_right=0.0):
def Poisson_dolfin(V_space, rho_net_array, Epsi_array, 
                   bc_left_type='dirichlet', V_left=0.0, 
                   bc_right_type='neumann', V_right=0.0):
    
    u = df.TrialFunction(V_space)
    v = df.TestFunction(V_space)

    e = 1.602176487E-19
    Epsi0 = 8.854187817620E-12

    v2d = df.vertex_to_dof_map(V_space)
    rho_func = df.Function(V_space)
    rho_func.vector()[v2d] = rho_net_array

    Epsi_func = df.Function(V_space)
    Epsi_func.vector()[v2d] = Epsi_array

    coeff = df.Constant(e / Epsi0)

    a_form = Epsi_func * df.inner(df.grad(u), df.grad(v)) * df.dx
    L = -coeff * rho_func * v * df.dx

    z_local = V_space.mesh().coordinates()[:, 0]
    tol_bc = 1e-14

    bcs = []
    
    # Left Boundary Condition
    if bc_left_type.lower() == 'dirichlet':
        def left_boundary(x, on_boundary):
            return on_boundary and abs(x[0] - z_local[0]) < tol_bc
        bcs.append(df.DirichletBC(V_space, df.Constant(V_left), left_boundary))

    # Right Boundary Condition
    if bc_right_type.lower() == 'dirichlet':
        def right_boundary(x, on_boundary):
            return on_boundary and abs(x[0] - z_local[-1]) < tol_bc
        bcs.append(df.DirichletBC(V_space, df.Constant(V_right), right_boundary))

    Vs_func = df.Function(V_space)
    df.solve(a_form == L, Vs_func, bcs)

    return Vs_func.vector().get_local()[v2d]
    