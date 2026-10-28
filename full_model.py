"""Unmodified restored-inertia class and closure from original combined diagnostic."""

from __future__ import annotations

import numpy as np

from scipy.integrate import simpson

import modal_model as n1

class FullerMultiModeModel(n1.MultiModeModel):
    """N1 multimode/direct-aero model with the N4 inertial terms restored."""

    def __init__(self, rw: float, n_modes: int = 3, nspan: int = n1.MAIN_NSPAN,
                 nx: int = n1.MAIN_NX):
        super().__init__(rw, n_modes, nspan=nspan, nx=nx)
        self.delta = self.chi*self.sigma

    def rhs(self, _t, y, mu):
        theta, theta_dot, q, q_dot, f, f_dot = self.split(y)
        a_theta, a_q = self.aerodynamic_loads_scalar(theta_dot, q_dot)

        generalized_rhs = np.empty(self.n_modes+1)
        generalized_rhs[0] = (
            mu*f-theta-a_theta
            - 2.0*theta_dot*np.dot(self.delta*q, q_dot)
        )
        generalized_rhs[1:] = (
            -2.0*n1.ZETA*self.r*q_dot
            -(self.r**2-self.sigma*theta_dot**2)*q-a_q
        )

        mass = self.mass_matrix.copy()
        mass[0, 0] += np.dot(self.delta, q*q)
        acc = np.linalg.solve(mass, generalized_rhs)

        f_ddot = -self.P["a2"]*f_dot-self.P["a3"]*f-self.P["a3"]*theta_dot
        out = np.empty_like(y)
        out[0], out[1] = theta_dot, acc[0]
        out[2:2+2*self.n_modes:2] = q_dot
        out[3:2+2*self.n_modes:2] = acc[1:]
        out[-2], out[-1] = f_dot, f_ddot
        return out

def energy_closure(result: dict, model: FullerMultiModeModel, mu: float) -> dict:
    """Check the fuller retained kinetic energy over the saved whole cycle."""
    t = result["t"]
    theta = result["theta"]
    theta_dot = result["theta_dot"]
    q = result["q"]
    q_dot = result["q_dot"]
    f = result["muscle_state"]
    beta = model.beta[:, None]
    chi = model.chi[:, None]
    delta = model.delta[:, None]

    energy = 0.5*(theta_dot**2+theta**2)
    energy += theta_dot*np.sum(beta*q_dot, axis=0)
    energy += 0.5*np.sum(chi*(q_dot**2+model.r[:, None]**2*q**2), axis=0)
    energy += 0.5*theta_dot**2*np.sum(delta*q**2, axis=0)

    p_muscle = mu*f*theta_dot
    a_theta, a_q = model.aerodynamic_loads_vector(theta_dot, q_dot)
    p_aero = a_theta*theta_dot + np.sum(chi*a_q*q_dot, axis=0)
    p_struct = np.sum(
        2.0*chi*n1.ZETA*model.r[:, None]*q_dot**2, axis=0)
    integrate = lambda value: float(simpson(value, x=t))
    delta_e = float(energy[-1]-energy[0])
    w_net = integrate(p_muscle)
    w_aero = integrate(p_aero)
    w_struct = integrate(p_struct)
    residual = w_net-delta_e-w_aero-w_struct
    return {
        "Delta_E": delta_e,
        "W_m_net_reintegrated": w_net,
        "W_aero_reintegrated": w_aero,
        "W_struct_reintegrated": w_struct,
        "energy_residual": residual,
        "energy_closure_relative": abs(residual)/max(abs(w_net), 1e-14),
    }
