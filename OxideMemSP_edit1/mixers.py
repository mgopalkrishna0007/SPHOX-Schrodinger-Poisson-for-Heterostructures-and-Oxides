# mixers.py
"""
Mixing schemes for the self-consistent Poisson-Schrödinger loop.

Choose the scheme in the config with:
    "mixing_method": "anderson"   # Anderson mixing  (uses mix_m, mix_beta, mix_max_step, mix_beta_eps, eps_start, eps_tol)
    "mixing_method": "alpha"      # adaptive alpha   (uses alpha, alpha_min, alpha_max, grow_factor, shrink_factor)

Each scheme exposes the same three steps that the SC loop calls every iteration:
    update_epsilon(nloop, Epsi_new, Epsi_arr) -> (Epsi_arr, eps_err)   # field-dependent eps damping
    is_converged(Err, eps_err)                -> bool
    next_potential(Vin, Vout, Err)            -> Vin for the next iteration
"""
import numpy as np


# =============================================================================
# Anderson mixer (core algorithm)
# =============================================================================
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
            try:
                gamma = np.linalg.solve(A, dF.T @ f)
                x_new = x + self.beta * f - (dX + self.beta * dF) @ gamma
            except np.linalg.LinAlgError:
                # residual history has collapsed (e.g. already converged, all dF = 0): keep the plain step
                pass

        step = x_new - x                             # cap the largest potential change
        smax = np.max(np.abs(step))
        if smax > self.max_step:
            x_new = x + step * (self.max_step / smax)
        return x_new


# =============================================================================
# Scheme 1: Anderson mixing for V, fixed delayed damping for eps
# =============================================================================
class AndersonMixingScheme:
    name = "anderson"

    def __init__(self, config):
        self.mixer = AndersonMixer(m=config["mix_m"], beta=config["mix_beta"], max_step=config["mix_max_step"])
        self.beta_eps = config["mix_beta_eps"]
        self.eps_start = config["eps_start"]
        self.tol = config["tol"]
        self.eps_tol = config["eps_tol"]
        self.prev_Err = None

    def update_epsilon(self, nloop, Epsi_new, Epsi_arr):
        # ε update: its own fixed damping, delayed so V settles first
        if nloop > self.eps_start:
            eps_err = np.max(np.abs(Epsi_new - Epsi_arr) / Epsi_arr)
            Epsi_arr = self.beta_eps * Epsi_new + (1 - self.beta_eps) * Epsi_arr
        else:
            eps_err = 0.0
        return Epsi_arr, eps_err

    def is_converged(self, Err, eps_err):
        return Err < self.tol and eps_err < self.eps_tol

    def next_potential(self, Vin, Vout, Err):
        if self.prev_Err is not None and Err > 5.0 * self.prev_Err:    # blow-up -> drop stale history
            self.mixer.reset()
        self.prev_Err = Err
        return self.mixer.mix(Vin, Vout)


# =============================================================================
# Scheme 2: adaptive alpha (linear) mixing for both V and eps
# =============================================================================
class AlphaMixingScheme:
    name = "alpha"

    def __init__(self, config):
        self.alpha = config["alpha"]
        self.alpha_min = config["alpha_min"]
        self.alpha_max = config["alpha_max"]
        self.grow_factor = config["grow_factor"]
        self.shrink_factor = config["shrink_factor"]
        self.tol = config["tol"]
        self.prev_Err = None

    def update_epsilon(self, nloop, Epsi_new, Epsi_arr):
        # ε is damped with the current alpha every loop; not part of the convergence check
        Epsi_arr = self.alpha * Epsi_new + (1 - self.alpha) * Epsi_arr
        return Epsi_arr, 0.0

    def is_converged(self, Err, eps_err):
        return Err < self.tol

    def next_potential(self, Vin, Vout, Err):
        if self.prev_Err is not None:
            if Err > self.prev_Err:
                self.alpha = max(self.alpha * self.shrink_factor, self.alpha_min)
            else:
                self.alpha = min(self.alpha * self.grow_factor, self.alpha_max)
        self.prev_Err = Err
        return self.alpha * Vout + (1 - self.alpha) * Vin


MIXING_SCHEMES = {
    "anderson": AndersonMixingScheme,
    "alpha": AlphaMixingScheme,
}


def create_mixing_scheme(config):
    """Build the mixing scheme selected by config["mixing_method"] (default: "anderson")."""
    method = str(config.get("mixing_method", "anderson")).lower()
    if method not in MIXING_SCHEMES:
        raise ValueError(
            f"Unknown mixing_method '{method}'. Choose one of: {', '.join(MIXING_SCHEMES)}"
        )
    return MIXING_SCHEMES[method](config)
