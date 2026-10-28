"""Unmodified modal methods extracted from the original N1 diagnostic; corrected local BEE input."""

from __future__ import annotations

import numpy as np

import math

from scipy.integrate import simpson, solve_ivp

from parameters import BEE

from model import dimensionless_groups

THETA = float(BEE["Theta"])

ZETA = float(BEE["zeta"])

BETAS = np.array([
    1.875104068711961,
    4.694091132974175,
    7.854757438237612,
    10.99554073487547,
])

MAIN_NSPAN = 64

MAIN_NX = 401

T_END = 600.0

TRANSIENT_FRACTION = 0.75

N_CYCLE = 6001

Y0_STROKE = 0.02

def raw_mode(beta: float, xi: np.ndarray) -> np.ndarray:
    sigma = (np.cosh(beta)+np.cos(beta))/(np.sinh(beta)+np.sin(beta))
    return (np.cosh(beta*xi)-np.cos(beta*xi)
            - sigma*(np.sinh(beta*xi)-np.sin(beta*xi)))

def raw_mode_d1(beta: float, xi: np.ndarray) -> np.ndarray:
    sigma = (np.cosh(beta)+np.cos(beta))/(np.sinh(beta)+np.sin(beta))
    return beta*(np.sinh(beta*xi)+np.sin(beta*xi)
                 - sigma*(np.cosh(beta*xi)-np.cos(beta*xi)))

def raw_mode_d2(beta: float, xi: np.ndarray) -> np.ndarray:
    sigma = (np.cosh(beta)+np.cos(beta))/(np.sinh(beta)+np.sin(beta))
    return beta**2*(np.cosh(beta*xi)+np.cos(beta*xi)
                    - sigma*(np.sinh(beta*xi)+np.sin(beta*xi)))

def raw_mode_d3(beta: float, xi: np.ndarray) -> np.ndarray:
    sigma = (np.cosh(beta)+np.cos(beta))/(np.sinh(beta)+np.sin(beta))
    return beta**3*(np.sinh(beta*xi)-np.sin(beta*xi)
                    - sigma*(np.cosh(beta*xi)+np.cos(beta*xi)))

def mode_values(n_modes: int, xi: np.ndarray, derivative: int = 0) -> np.ndarray:
    rows = []
    for beta in BETAS[:n_modes]:
        tip = float(raw_mode(beta, np.array([1.0]))[0])
        if derivative == 0:
            value = raw_mode(beta, xi)
        elif derivative == 1:
            value = raw_mode_d1(beta, xi)
        elif derivative == 2:
            value = raw_mode_d2(beta, xi)
        elif derivative == 3:
            value = raw_mode_d3(beta, xi)
        else:
            raise ValueError("derivative must be 0, 1, 2, or 3")
        rows.append(value/tip)
    return np.asarray(rows)

def gauss_rule(n: int):
    x, w = np.polynomial.legendre.leggauss(int(n))
    return 0.5*(x+1.0), 0.5*w

def modal_data(n_modes: int, nquad: int = 400):
    xi, weights = gauss_rule(nquad)
    phi = mode_values(n_modes, xi, 0)
    dphi = mode_values(n_modes, xi, 1)
    ddphi = mode_values(n_modes, xi, 2)
    mass = (phi*weights) @ phi.T
    stiffness = (ddphi*weights) @ ddphi.T
    j0 = np.diag(mass)
    j1 = phi @ (weights*xi)
    ratios = (BETAS[:n_modes]/BETAS[0])**2
    return dict(xi=xi, weights=weights, phi=phi, dphi=dphi,
                mass=mass, stiffness=stiffness, J0=j0, J1=j1,
                frequency_ratios=ratios)

def modal_verification_rows():
    data = modal_data(4, 600)
    mass = data["mass"]
    stiff = data["stiffness"]
    j0 = data["J0"]
    mass_rel = mass/np.sqrt(np.outer(j0, j0))
    stiff_rel = stiff/np.sqrt(np.outer(np.diag(stiff), np.diag(stiff)))
    xends = np.array([0.0, 1.0])
    phi_ends = mode_values(4, xends, 0)
    dphi_ends = mode_values(4, xends, 1)
    ddphi_ends = mode_values(4, xends, 2)
    dddphi_ends = mode_values(4, xends, 3)
    rows = []
    for i in range(4):
        off_mass = max(abs(mass_rel[i, j]) for j in range(4) if j != i)
        off_stiff = max(abs(stiff_rel[i, j]) for j in range(4) if j != i)
        rows.append(dict(
            mode=i+1, beta=BETAS[i], tip_value=phi_ends[i, 1],
            frequency_ratio=data["frequency_ratios"][i],
            J0=j0[i], J1=data["J1"][i],
            stiffness_over_mass=stiff[i, i]/j0[i],
            beta_fourth=BETAS[i]**4,
            relative_frequency_error=(
                math.sqrt(stiff[i, i]/j0[i])/BETAS[i]**2-1.0),
            max_mass_offdiag_relative=off_mass,
            max_stiffness_offdiag_relative=off_stiff,
            root_displacement=phi_ends[i, 0],
            root_slope=dphi_ends[i, 0],
            tip_moment=ddphi_ends[i, 1],
            tip_shear=dddphi_ends[i, 1],
        ))
    summary = dict(
        max_mass_offdiag_relative=float(np.max(np.abs(mass_rel-np.eye(4)))),
        max_stiffness_offdiag_relative=float(np.max(np.abs(stiff_rel-np.eye(4)))),
        max_frequency_relative_error=float(max(abs(r["relative_frequency_error"])
                                               for r in rows)),
        max_boundary_condition_residual=float(max(
            np.max(np.abs(phi_ends[:, 0])), np.max(np.abs(dphi_ends[:, 0])),
            np.max(np.abs(ddphi_ends[:, 1])), np.max(np.abs(dddphi_ends[:, 1])))),
    )
    return rows, summary

class MultiModeModel:
    def __init__(self, rw: float, n_modes: int, nspan: int = MAIN_NSPAN,
                 nx: int = MAIN_NX):
        self.rw = float(rw)
        self.n_modes = int(n_modes)
        self.nspan = int(nspan)
        self.nx = int(nx)
        self.P = dimensionless_groups(BEE, self.rw)
        md = modal_data(self.n_modes)
        self.J0 = md["J0"]
        self.J1 = md["J1"]
        self.frequency_ratios = md["frequency_ratios"]
        self.r = self.rw*self.frequency_ratios
        lam = self.P["lambda_t"]
        self.beta = 3.0*self.J1/((1.0+lam)*THETA)
        self.gamma = self.J1*THETA/self.J0
        self.chi = self.beta/self.gamma
        self.sigma = THETA**2

        # The generalized mass matrix contains no bending-mode cross terms
        # because exact uniform-beam modes are mass orthogonal.
        mass = np.eye(self.n_modes+1)
        mass[0, 1:] = self.beta
        mass[1:, 0] = self.gamma
        self.mass_matrix = mass
        self.mass_matrix_inverse = np.linalg.inv(mass)

        self.xi, self.weights = gauss_rule(self.nspan)
        self.phi = mode_values(self.n_modes, self.xi, 0)
        self.x_eval = np.linspace(0.0, 1.0, self.nx)
        self.phi_eval = mode_values(self.n_modes, self.x_eval, 0)
        self.dphi_eval = mode_values(self.n_modes, self.x_eval, 1)

        # First-mode coefficient rewritten so every mode is normalized by its
        # own modal mass.  This exactly recovers the N5 one-mode formula.
        first_md = modal_data(1)
        j0_first = float(first_md["J0"][0])
        j2_first = float(np.dot(first_md["weights"],
                                first_md["xi"]**2*first_md["phi"][0]))
        k_first = self.P["eta_B"]/(THETA**2*j2_first)
        self.aero_q_scale = k_first*j0_first/self.J0

    @property
    def state_size(self):
        return 2+2*self.n_modes+2

    def initial_state(self):
        y = np.zeros(self.state_size)
        y[0] = Y0_STROKE
        return y

    def split(self, y):
        y = np.asarray(y)
        theta, theta_dot = y[0], y[1]
        q = y[2:2+2*self.n_modes:2]
        q_dot = y[3:2+2*self.n_modes:2]
        f, f_dot = y[-2], y[-1]
        return theta, theta_dot, q, q_dot, f, f_dot

    def aerodynamic_loads_scalar(self, theta_dot: float, q_dot: np.ndarray):
        u = THETA*self.xi*theta_dot + q_dot @ self.phi
        drag = u*np.abs(u)
        a_theta = (4.0/(self.P["Nw"]*THETA**2))*np.dot(
            self.weights, self.xi*drag)
        a_q = self.aero_q_scale*(self.phi @ (self.weights*drag))
        return float(a_theta), a_q

    def aerodynamic_loads_vector(self, theta_dot: np.ndarray,
                                  q_dot: np.ndarray):
        u = THETA*self.xi[:, None]*theta_dot[None, :]
        u = u + self.phi.T @ q_dot
        drag = u*np.abs(u)
        a_theta = (4.0/(self.P["Nw"]*THETA**2))*np.sum(
            self.weights[:, None]*self.xi[:, None]*drag, axis=0)
        a_q = self.aero_q_scale[:, None]*(
            self.phi @ (self.weights[:, None]*drag))
        return a_theta, a_q

    def rhs(self, _t, y, mu):
        theta, theta_dot, q, q_dot, f, f_dot = self.split(y)
        a_theta, a_q = self.aerodynamic_loads_scalar(theta_dot, q_dot)
        generalized_rhs = np.empty(self.n_modes+1)
        generalized_rhs[0] = mu*f-theta-a_theta
        generalized_rhs[1:] = (
            -2.0*ZETA*self.r*q_dot
            -(self.r**2-self.sigma*theta_dot**2)*q-a_q)
        acc = self.mass_matrix_inverse @ generalized_rhs
        f_ddot = -self.P["a2"]*f_dot-self.P["a3"]*f-self.P["a3"]*theta_dot
        out = np.empty_like(y)
        out[0], out[1] = theta_dot, acc[0]
        out[2:2+2*self.n_modes:2] = q_dot
        out[3:2+2*self.n_modes:2] = acc[1:]
        out[-2], out[-1] = f_dot, f_ddot
        return out

def upward_zero(_t, y, *_args):
    return float(y[0])

upward_zero.direction = 1

upward_zero.terminal = False

def phase_difference_deg(t, x, y, period):
    tau = t-t[0]
    kernel = np.exp(-1j*2.0*np.pi*tau/period)
    cx = simpson(x*kernel, x=t)
    cy = simpson(y*kernel, x=t)
    value = np.degrees(np.angle(cy)-np.angle(cx))
    return float((value+180.0) % 360.0-180.0)

def cycle_metrics(sol, ta, tb, model: MultiModeModel, mu: float,
                  n_cycle: int = N_CYCLE):
    t = np.linspace(ta, tb, int(n_cycle))
    state = sol.sol(t)
    theta, theta_dot = state[0], state[1]
    q = state[2:2+2*model.n_modes:2]
    q_dot = state[3:2+2*model.n_modes:2]
    f = state[-2]
    a_theta, a_q = model.aerodynamic_loads_vector(theta_dot, q_dot)
    p_aero = a_theta*theta_dot + np.sum(
        model.chi[:, None]*a_q*q_dot, axis=0)
    p_muscle = mu*f*theta_dot
    p_struct = np.sum(
        2.0*model.chi[:, None]*ZETA*model.r[:, None]*q_dot**2, axis=0)

    w_field = model.phi_eval.T @ q
    slope_field = model.dphi_eval.T @ q
    max_deflection = float(np.max(np.abs(w_field)))
    max_slope = float(np.max(np.abs(slope_field)))
    tip = np.sum(q, axis=0)
    period = float(tb-ta)
    integrate = lambda value: float(simpson(value, x=t))

    modal_energy = 0.5*model.chi[:, None]*(
        q_dot**2+model.r[:, None]**2*q**2
        + model.sigma*theta_dot[None, :]**2*q**2)
    modal_energy_mean = np.array([
        integrate(modal_energy[i])/period for i in range(model.n_modes)])
    modal_energy_fraction = modal_energy_mean/max(np.sum(modal_energy_mean), 1e-30)
    modal_mass_displacement = model.J0[:, None]*q**2
    modal_rms_share = np.mean(modal_mass_displacement, axis=1)
    modal_rms_fraction = modal_rms_share/max(np.sum(modal_rms_share), 1e-30)

    state_a, state_b = sol.sol(ta), sol.sol(tb)
    state_scale = max(float(np.max(np.abs(state))), 1e-12)
    result = dict(
        period=period, Omega=float(2.0*np.pi/period),
        stroke_amp=float(0.5*(np.max(theta)-np.min(theta))),
        max_deflection=max_deflection, max_slope=max_slope,
        rotation_deg=float(np.degrees(np.arctan(max_slope))),
        tip_amplitude=float(np.max(np.abs(tip))),
        phase_tip_minus_theta_deg=phase_difference_deg(t, theta, tip, period),
        peak_aero_stroke_load=float(np.max(np.abs(a_theta))),
        peak_aero_bend_load_scaled=float(np.max(np.sqrt(np.sum(
            (model.chi[:, None]*a_q)**2, axis=0)))),
        peak_aero_power=float(np.max(np.abs(p_aero))),
        min_aero_power=float(np.min(p_aero)), W_aero=integrate(p_aero),
        W_plus=integrate(np.maximum(p_muscle, 0.0)),
        W_minus=integrate(np.maximum(-p_muscle, 0.0)),
        W_m_net=integrate(p_muscle), W_struct=integrate(p_struct),
        peak_muscle_force=float(np.max(np.abs(mu*f))),
        peak_muscle_power=float(np.max(np.abs(p_muscle))),
        modal_amplitudes=np.max(np.abs(q), axis=1),
        modal_tip_phases_deg=np.array([
            phase_difference_deg(t, theta, q[i], period)
            for i in range(model.n_modes)]),
        modal_energy_mean=modal_energy_mean,
        modal_energy_fraction=modal_energy_fraction,
        modal_rms_deformation_fraction=modal_rms_fraction,
        periodic_state_abs=float(np.max(np.abs(state_b-state_a))),
        periodic_state_rel=float(np.max(np.abs(state_b-state_a))/state_scale),
        t=t, theta=theta, theta_dot=theta_dot, q=q, q_dot=q_dot,
        tip=tip, muscle_state=f, aero_stroke_load=a_theta,
        aero_bend_load=a_q, aero_power=p_aero,
    )
    return result

def solve_cycle(rw: float, n_modes: int, mu: float, nspan: int = MAIN_NSPAN,
                nx: int = MAIN_NX, t_end: float = T_END, rtol: float = 3e-9,
                atol: float = 3e-11, max_step: float = 0.045,
                n_cycle: int = N_CYCLE):
    model = MultiModeModel(rw, n_modes, nspan=nspan, nx=nx)
    sol = solve_ivp(
        model.rhs, (0.0, float(t_end)), model.initial_state(), args=(mu,),
        method="DOP853", rtol=rtol, atol=atol, max_step=max_step,
        dense_output=True, events=upward_zero,
    )
    if not sol.success:
        raise RuntimeError(f"{n_modes}-mode integration failed: {sol.message}")
    crossings = sol.t_events[0]
    crossings = crossings[crossings > TRANSIENT_FRACTION*t_end]
    if len(crossings) < 4:
        raise RuntimeError(f"{n_modes}-mode: insufficient steady cycles")
    result = cycle_metrics(sol, float(crossings[-2]), float(crossings[-1]),
                           model, mu, n_cycle=n_cycle)
    previous = cycle_metrics(sol, float(crossings[-3]), float(crossings[-2]),
                             model, mu, n_cycle=1601)
    result["cycle_amp_rel"] = abs(
        result["stroke_amp"]-previous["stroke_amp"])/max(result["stroke_amp"], 1e-12)
    result["cycle_work_rel"] = abs(
        result["W_m_net"]-previous["W_m_net"])/max(abs(result["W_m_net"]), 1e-12)
    result["solver_nfev"] = int(sol.nfev)
    result["solution"] = sol
    result["model"] = model
    return result
