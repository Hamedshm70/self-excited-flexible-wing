"""Retained structural equations, dSA dynamics, and aerodynamic operators.
The direct operator is the production model. The reduced operator and legacy
RK45 integration are retained ONLY for the Figure 4 provenance comparison.
I_theta returned by dimensionless_groups is HALF the bilateral inertia.
"""
from __future__ import annotations
import numpy as np
from scipy.integrate import simpson, solve_ivp
from scipy.optimize import brentq
from parameters import (BEE, SLOPE_FACTOR, NSPAN as MAIN_NSPAN, T_END,
                        TRANSIENT, N_CYCLE, Y0, THETA)
BETA1 = 1.8751040687
SIG1 = 0.7340955

def _phi_raw(xi):
    b = BETA1
    return np.cosh(b * xi) - np.cos(b * xi) - SIG1 * (np.sinh(b * xi) - np.sin(b * xi))
_PHI_TIP = _phi_raw(1.0)

def phi(xi):
    return _phi_raw(xi) / _PHI_TIP
_xi = np.linspace(0.0, 1.0, 20001)
_ph = phi(_xi)
I0 = np.trapz(_ph ** 2, _xi)
I1 = np.trapz(_xi * _ph, _xi)
I2 = np.trapz(_xi ** 2 * _ph, _xi)
I3 = np.trapz(_xi * _ph ** 2, _xi)
def dimensionless_groups(sp, r_w):
    """Build all dimensionless groups from physical species params + chosen r_w.
    kappa_a = aero blade-element prefactor (0.5 = textbook one-wing 1/2 factor).
    Added mass optionally folded into the effective wing mass via sp['added_mass'].
    Thoracic inertia folded in via sp['lambda_t'] = m*Gamma^2/(2 I_w)."""
    (L, c, mw) = (sp['L'], sp['c'], sp['mw'])
    (CD, rho, Theta) = (sp['CD'], sp['rho'], sp['Theta'])
    ka = sp.get('kappa_a', 0.5)
    lam = sp.get('lambda_t', 0.0)
    if sp.get('added_mass', False):
        ma = np.pi / 4.0 * rho * c ** 2
        mw_eff = mw + ma * L
    else:
        mw_eff = mw
    Iw = mw_eff * L ** 2 / 3.0
    I_theta = (1.0 + lam) * Iw
    A3 = c * L ** 4 / 4.0
    Nw = I_theta / (ka * CD * rho * A3 * Theta)
    b_tq = 3.0 * I1 / ((1.0 + lam) * Theta)
    g_qt = I1 * Theta / I0
    sig_cf = Theta ** 2
    eta_A = 8.0 * I2 / (Theta * Nw)
    eta_B = ka * CD * rho * c * L ** 2 * I2 * Theta ** 2 / (mw_eff * I0)
    eta_C = 2.0 * ka * CD * rho * c * L ** 2 * I3 * Theta / (mw_eff * I0)
    return dict(Nw=Nw, r=r_w, b_tq=b_tq, g_qt=g_qt, sig_cf=sig_cf, eta_A=eta_A, eta_B=eta_B, eta_C=eta_C, I_theta=I_theta, A3=A3, lambda_t=lam, a2=sp['a2'], a3=sp['a3'], zeta=sp['zeta'], zeta_theta=sp.get('zeta_theta', 0.0))

def _rhs(t, y, P, mu, flexible):
    (th, thd, q, qd, f, fd) = y
    (Nw, r, b, g) = (P['Nw'], P['r'], P['b_tq'], P['g_qt'])
    (sig, eA, eB, eC, z) = (P['sig_cf'], P['eta_A'], P['eta_B'], P['eta_C'], P['zeta'])
    (a2, a3) = (P['a2'], P['a3'])
    z_theta = P.get('zeta_theta', 0.0)
    RHS_th = mu * f - th - 2.0 * z_theta * thd - 1.0 / Nw * thd * abs(thd)
    if flexible:
        RHS_th -= eA * abs(thd) * qd
        RHS_q = -(2 * z * r + eC * abs(thd)) * qd - (r ** 2 - sig * thd ** 2) * q - eB * thd * abs(thd)
        det = 1.0 - b * g
        thdd = (RHS_th - b * RHS_q) / det
        qdd = (RHS_q - g * RHS_th) / det
    else:
        thdd = RHS_th
        (q, qd, qdd) = (0.0, 0.0, 0.0)
    fdd = -a2 * fd - a3 * f - a3 * thd
    return [thd, thdd, qd, qdd, fd, fdd]

def cycle_energy_metrics(t, th, thd, q, qd, f, P, mu, flexible, period):
    """Return cycle-averaged powers and per-cycle works on a whole-cycle segment.

    ``t`` must span an integer number of steady cycles.  ``W_minus`` is the
    positive magnitude of negative muscle work.  The aerodynamic channel is the
    complete reduced-order expression, including the eta_A stroke reaction.
    Legacy aliases are intentionally not created here; callers should use the
    explicit ``W_*`` and ``phi_*`` names.
    """
    (Nw, r) = (P['Nw'], P['r'])
    z_theta = P.get('zeta_theta', 0.0)
    if flexible:
        (beta, gamma) = (P['b_tq'], P['g_qt'])
        (sigma, eta_a) = (P['sig_cf'], P['eta_A'])
        (eta_b, eta_c, zeta) = (P['eta_B'], P['eta_C'], P['zeta'])
        rhs_th = mu * f - th - 2.0 * z_theta * thd - 1.0 / Nw * thd * np.abs(thd) - eta_a * np.abs(thd) * qd
        rhs_q = -(2.0 * zeta * r + eta_c * np.abs(thd)) * qd
        rhs_q -= (r ** 2 - sigma * thd ** 2) * q + eta_b * thd * np.abs(thd)
        det = 1.0 - beta * gamma
        thdd = (rhs_th - beta * rhs_q) / det
        qdd = (rhs_q - gamma * rhs_th) / det
        chi = beta / gamma
    else:
        beta = gamma = sigma = eta_a = eta_b = eta_c = zeta = chi = 0.0
        thdd = mu * f - th - 2.0 * z_theta * thd - 1.0 / Nw * thd * np.abs(thd)
        qdd = np.zeros_like(thdd)
    p_muscle = mu * f * thd
    p_thorax = 2.0 * z_theta * thd ** 2
    p_aero_primary = 1.0 / Nw * np.abs(thd) ** 3
    p_aero_eta_a = eta_a * thd * np.abs(thd) * qd
    p_aero_eta_b = chi * eta_b * thd * np.abs(thd) * qd
    p_aero_eta_c = chi * eta_c * np.abs(thd) * qd ** 2
    p_aero = p_aero_primary + p_aero_eta_a + p_aero_eta_b + p_aero_eta_c
    p_struct = 2.0 * chi * zeta * r * qd ** 2
    p_trunc = -chi * sigma * (q ** 2 * thd * thdd + 2.0 * q * qd * thd ** 2)
    duration = float(t[-1] - t[0])
    n_cycles = duration / period if period > 0 else np.nan
    integrate = lambda x: float(np.trapezoid(x, t) / n_cycles)
    w_plus = integrate(np.maximum(p_muscle, 0.0))
    w_minus = integrate(np.maximum(-p_muscle, 0.0))
    w_m_net = w_plus - w_minus
    w_abs = w_plus + w_minus
    w_thorax = integrate(p_thorax)
    w_aero_primary = integrate(p_aero_primary)
    w_aero_eta_a = integrate(p_aero_eta_a)
    w_aero_eta_b = integrate(p_aero_eta_b)
    w_aero_eta_c = integrate(p_aero_eta_c)
    w_aero = integrate(p_aero)
    w_struct = integrate(p_struct)
    w_trunc = integrate(p_trunc)
    energy = 0.5 * (thd ** 2 + 2.0 * beta * thd * qd + chi * qd ** 2 + th ** 2 + chi * r ** 2 * q ** 2 + chi * sigma * thd ** 2 * q ** 2)
    delta_e = float((energy[-1] - energy[0]) / n_cycles)
    w_num = w_m_net - delta_e - w_thorax - w_aero - w_struct - w_trunc
    p_flex_load = beta * qdd * thd + p_aero_eta_a
    mean_complete_aero = w_aero / period
    pi_flex = float(np.max(np.abs(p_flex_load))) / mean_complete_aero if mean_complete_aero > 0 else np.nan
    phi_flex_load = float(np.trapezoid(np.abs(p_flex_load), t) / duration) / mean_complete_aero if mean_complete_aero > 0 else np.nan
    modal_energy = 0.5 * chi * (qd ** 2 + r ** 2 * q ** 2 + sigma * thd ** 2 * q ** 2)
    return dict(W_plus=w_plus, W_minus=w_minus, W_abs=w_abs, W_m_net=w_m_net, W_aero=w_aero, W_aero_primary=w_aero_primary, W_aero_etaA=w_aero_eta_a, W_aero_etaB=w_aero_eta_b, W_aero_etaC=w_aero_eta_c, W_thorax=w_thorax, W_struct=w_struct, W_trunc=w_trunc, Delta_E=delta_e, W_num=w_num, phi_net=w_m_net / w_plus if w_plus > 0 else np.nan, phi_aero=w_aero / w_plus if w_plus > 0 else np.nan, mean_P_plus=w_plus / period, mean_P_minus=w_minus / period, mean_P_abs=w_abs / period, mean_P_m_net=w_m_net / period, mean_P_aero=w_aero / period, mean_P_thorax=w_thorax / period, mean_P_struct=w_struct / period, mean_P_trunc=w_trunc / period, peak_muscle_force=float(np.max(np.abs(mu * f))), peak_muscle_power=float(np.max(np.abs(p_muscle))), peak_positive_muscle_power=float(np.max(p_muscle)), peak_negative_muscle_power=float(-np.min(p_muscle)), peak_flex_load=float(np.max(np.abs(p_flex_load))), Pi_flex=pi_flex, Phi_flex_load=phi_flex_load, modal_energy_max=float(np.max(modal_energy)) if flexible else 0.0, closure_rel=abs(w_num) / max(abs(w_m_net), 1e-14))

def simulate(P, mu, flexible=True, t_end=400.0, n_per=4000, *, method='RK45', rtol=1e-06, atol=1e-08, max_step=np.inf):
    y0 = [0.02, 0.0, 0.0, 0.0, 0.0, 0.0]
    t_eval = np.linspace(0, t_end, n_per)
    sol = solve_ivp(_rhs, (0, t_end), y0, args=(P, mu, flexible), t_eval=t_eval, method=method, rtol=rtol, atol=atol, max_step=max_step)
    (t, Y) = (sol.t, sol.y)
    m = t > 0.7 * t_end
    (th, thd, q, qd, f) = (Y[0][m], Y[1][m], Y[2][m], Y[3][m], Y[4][m])
    ts = t[m]
    if th.size < 5 or not np.all(np.isfinite(th)):
        return dict(t=ts, th=th, thd=thd, q=q, qd=qd, f=f, amp=np.nan, Omega=np.nan, period=np.nan, phi_net=np.nan, phi_aero=np.nan, mu=mu, qamp=np.nan)
    amp = 0.5 * (th.max() - th.min())
    s = np.sign(th - th.mean())
    cross = np.where((s[:-1] < 0) & (s[1:] >= 0))[0]
    if len(cross) >= 2:
        (i0, i1) = (cross[0], cross[-1])
        period = (ts[i1] - ts[i0]) / (len(cross) - 1)
        Omega = 2 * np.pi / period if period > 0 else np.nan
        em = cycle_energy_metrics(ts[i0:i1 + 1], th[i0:i1 + 1], thd[i0:i1 + 1], q[i0:i1 + 1], qd[i0:i1 + 1], f[i0:i1 + 1], P, mu, flexible, period)
    else:
        (period, Omega, em) = (np.nan, np.nan, {})
    return dict(t=ts, th=th, thd=thd, q=q, qd=qd, f=f, amp=amp, Omega=Omega, period=period, mu=mu, qamp=0.5 * (q.max() - q.min()) if flexible else 0.0, **em)
def gauss_rule(n: int):
    (x, w) = np.polynomial.legendre.leggauss(n)
    xi = 0.5 * (x + 1.0)
    return (xi, 0.5 * w, phi(xi))

class Aerodynamics:

    def __init__(self, formulation: str, nspan: int=MAIN_NSPAN):
        self.formulation = formulation
        self.nspan = int(nspan)
        (self.xi, self.w, self.phi) = gauss_rule(self.nspan)

    def scalar(self, thd: float, qd: float, P: dict, flexible: bool):
        if self.formulation == 'reduced':
            ath = 1.0 / P['Nw'] * thd * abs(thd)
            aq = 0.0
            if flexible:
                ath += P['eta_A'] * abs(thd) * qd
                aq = P['eta_B'] * thd * abs(thd) + P['eta_C'] * abs(thd) * qd
            return (float(ath), float(aq))
        u = THETA * self.xi * thd
        if flexible:
            u = u + self.phi * qd
        drag = u * np.abs(u)
        ath = 4.0 / (P['Nw'] * THETA ** 2) * np.dot(self.w, self.xi * drag)
        aq = 0.0
        if flexible:
            kq = P['eta_B'] / (THETA ** 2 * I2)
            aq = kq * np.dot(self.w, self.phi * drag)
        return (float(ath), float(aq))

    def vector(self, thd: np.ndarray, qd: np.ndarray, P: dict, flexible: bool):
        thd = np.asarray(thd)
        qd = np.asarray(qd)
        if self.formulation == 'reduced':
            ath = 1.0 / P['Nw'] * thd * np.abs(thd)
            aq = np.zeros_like(thd)
            if flexible:
                ath = ath + P['eta_A'] * np.abs(thd) * qd
                aq = P['eta_B'] * thd * np.abs(thd) + P['eta_C'] * np.abs(thd) * qd
            return (ath, aq)
        u = THETA * self.xi[:, None] * thd[None, :]
        if flexible:
            u = u + self.phi[:, None] * qd[None, :]
        drag = u * np.abs(u)
        ath = 4.0 / (P['Nw'] * THETA ** 2) * np.sum(self.w[:, None] * self.xi[:, None] * drag, axis=0)
        if flexible:
            kq = P['eta_B'] / (THETA ** 2 * I2)
            aq = kq * np.sum(self.w[:, None] * self.phi[:, None] * drag, axis=0)
        else:
            aq = np.zeros_like(ath)
        return (ath, aq)

def rhs(t, y, P, mu, flexible, aero: Aerodynamics):
    (th, thd, q, qd, f, fd) = y
    (ath, aq) = aero.scalar(thd, qd, P, flexible)
    rhs_th = mu * f - th - 2.0 * P.get('zeta_theta', 0.0) * thd - ath
    if flexible:
        rhs_q = -2.0 * P['zeta'] * P['r'] * qd
        rhs_q -= (P['r'] ** 2 - P['sig_cf'] * thd ** 2) * q + aq
        det = 1.0 - P['b_tq'] * P['g_qt']
        thdd = (rhs_th - P['b_tq'] * rhs_q) / det
        qdd = (rhs_q - P['g_qt'] * rhs_th) / det
    else:
        (thdd, qdd) = (rhs_th, 0.0)
    fdd = -P['a2'] * fd - P['a3'] * f - P['a3'] * thd
    return np.array([thd, thdd, qd if flexible else 0.0, qdd, fd, fdd])

def upward_zero(_t, y, *_args):
    return float(y[0])

def phase_difference_deg(t, x, y, period):
    tau = t - t[0]
    kernel = np.exp(-1j * 2.0 * np.pi * tau / period)
    cx = simpson(x * kernel, x=t)
    cy = simpson(y * kernel, x=t)
    value = np.degrees(np.angle(cy) - np.angle(cx))
    return float((value + 180.0) % 360.0 - 180.0)

def cycle_metrics(sol, ta, tb, P, mu, flexible, aero, n_cycle=N_CYCLE):
    t = np.linspace(ta, tb, int(n_cycle))
    (th, thd, q, qd, f, fd) = sol.sol(t)
    (ath, aq) = aero.vector(thd, qd, P, flexible)
    chi = P['b_tq'] / P['g_qt'] if flexible else 0.0
    p_aero = ath * thd + chi * aq * qd
    p_muscle = mu * f * thd
    p_struct = 2.0 * chi * P['zeta'] * P['r'] * qd ** 2 if flexible else np.zeros_like(t)
    dy = np.column_stack([rhs(float(ti), yi, P, mu, flexible, aero) for (ti, yi) in zip(t, sol.sol(t).T)])
    thdd = dy[1]
    p_trunc = -chi * P['sig_cf'] * (q ** 2 * thd * thdd + 2.0 * q * qd * thd ** 2) if flexible else np.zeros_like(t)
    integ = lambda value: float(simpson(value, x=t))
    wplus = integ(np.maximum(p_muscle, 0.0))
    wminus = integ(np.maximum(-p_muscle, 0.0))
    wnet = wplus - wminus
    waero = integ(p_aero)
    wstruct = integ(p_struct)
    wtrunc = integ(p_trunc)
    (state_a, state_b) = (sol.sol(ta), sol.sol(tb))
    state_scale = max(np.max(np.abs(sol.sol(np.linspace(ta, tb, 1001)))), 1e-10)
    qmax = float(np.max(np.abs(q))) if flexible else 0.0
    slope = SLOPE_FACTOR * qmax
    return dict(period=float(tb - ta), Omega=float(2.0 * np.pi / (tb - ta)), stroke_amp=float(0.5 * (np.max(th) - np.min(th))), q_max=qmax, max_slope=float(slope), rotation_deg=float(np.degrees(np.arctan(slope))), phase_q_minus_theta_deg=phase_difference_deg(t, th, q, tb - ta) if flexible else 0.0, peak_aero_stroke_load=float(np.max(np.abs(ath))), peak_aero_bend_load_scaled=float(np.max(np.abs(chi * aq))) if flexible else 0.0, peak_aero_power=float(np.max(np.abs(p_aero))), min_aero_dissipation_power=float(np.min(p_aero)), W_aero=waero, W_plus=wplus, W_minus=wminus, W_m_net=wnet, W_struct=wstruct, W_trunc=wtrunc, peak_muscle_force=float(np.max(np.abs(mu * f))), peak_muscle_power=float(np.max(np.abs(p_muscle))), peak_positive_muscle_power=float(np.max(p_muscle)), peak_negative_muscle_power=float(-np.min(p_muscle)), periodic_state_abs=float(np.max(np.abs(state_b - state_a))), periodic_state_rel=float(np.max(np.abs(state_b - state_a)) / state_scale), t=t, theta=th, theta_dot=thd, q=q, q_dot=qd, muscle_state=f, aero_stroke_load=ath, aero_bend_load=aq, aero_power=p_aero)

def solve_cycle(P, mu, flexible, formulation, nspan=MAIN_NSPAN, t_end=T_END, rtol=3e-09, atol=3e-11, max_step=0.06, n_cycle=N_CYCLE):
    aero = Aerodynamics(formulation, nspan)
    sol = solve_ivp(rhs, (0.0, t_end), Y0, args=(P, mu, flexible, aero), method='DOP853', rtol=rtol, atol=atol, max_step=max_step, dense_output=True, events=upward_zero)
    if not sol.success:
        raise RuntimeError(f'{formulation} integration failed: {sol.message}')
    crossings = sol.t_events[0]
    crossings = crossings[crossings > TRANSIENT * t_end]
    if len(crossings) < 4:
        raise RuntimeError(f'{formulation}: insufficient steady cycles')
    last = cycle_metrics(sol, float(crossings[-2]), float(crossings[-1]), P, mu, flexible, aero, n_cycle=n_cycle)
    previous = cycle_metrics(sol, float(crossings[-3]), float(crossings[-2]), P, mu, flexible, aero, n_cycle=2001)
    last['cycle_amp_rel'] = abs(last['stroke_amp'] - previous['stroke_amp']) / max(last['stroke_amp'], 1e-12)
    last['cycle_work_rel'] = abs(last['W_m_net'] - previous['W_m_net']) / max(abs(last['W_m_net']), 1e-12)
    last['solver_nfev'] = int(sol.nfev)
    last['solution'] = sol
    return last

def quick_amplitude(P, mu, flexible, formulation, nspan):
    result = solve_cycle(P, mu, flexible, formulation, nspan=nspan, t_end=420.0, rtol=1e-08, atol=1e-10, max_step=0.09, n_cycle=1201)
    return result['stroke_amp']

class TargetUnreachable(RuntimeError):

    def __init__(self, formulation, target, max_amp, max_mu):
        super().__init__(f'{formulation} target {target} not reached; max {max_amp} at {max_mu}')
        self.max_amp = max_amp
        self.max_mu = max_mu

def tune_drive(P, flexible, formulation, target, initial, nspan=MAIN_NSPAN):
    cache = {}

    def fun(mu):
        key = round(float(mu), 8)
        if key not in cache:
            cache[key] = quick_amplitude(P, float(mu), flexible, formulation, nspan) - target
        return cache[key]
    grid = np.linspace(max(0.25, 0.3 * initial), 1.65 * initial, 19)
    previous = None
    sampled = []
    for mu in grid:
        try:
            value = fun(float(mu))
        except RuntimeError:
            break
        sampled.append((float(mu), value + target))
        if previous is not None and previous[1] <= 0 <= value:
            return float(brentq(fun, previous[0], float(mu), xtol=0.0002, rtol=2e-05, maxiter=20))
        previous = (float(mu), value)
    if not sampled:
        raise TargetUnreachable(formulation, target, float('nan'), float('nan'))
    (max_mu, max_amp) = max(sampled, key=lambda item: item[1])
    raise TargetUnreachable(formulation, target, max_amp, max_mu)

def coefficient_checks(P):
    rows = []
    for n in (8, 16, 32, 64, 128):
        (xi, w, ph) = gauss_rule(n)
        j2 = float(np.dot(w, xi ** 2 * ph))
        j3 = float(np.dot(w, xi * ph ** 2))
        primary_th = 4.0 / (P['Nw'] * THETA ** 2) * THETA ** 2 * np.dot(w, xi ** 3)
        cross_th = 4.0 / (P['Nw'] * THETA ** 2) * 2.0 * THETA * j2
        primary_q = P['eta_B'] / (THETA ** 2 * I2) * THETA ** 2 * j2
        cross_q = P['eta_B'] / (THETA ** 2 * I2) * 2.0 * THETA * j3
        rows.append(dict(nspan=n, J2=j2, J3=j3, primary_theta_rel_error=abs(primary_th - 1.0 / P['Nw']) / (1.0 / P['Nw']), cross_theta_rel_error=abs(cross_th - P['eta_A']) / P['eta_A'], primary_q_rel_error=abs(primary_q - P['eta_B']) / P['eta_B'], cross_q_rel_error=abs(cross_q - P['eta_C']) / P['eta_C']))
    (red, full) = (Aerodynamics('reduced', 64), Aerodynamics('direct', 64))
    rigid_err = []
    for thd in np.linspace(-3.0, 3.0, 61):
        ar = red.scalar(float(thd), 0.0, P, False)[0]
        af = full.scalar(float(thd), 0.0, P, False)[0]
        rigid_err.append(abs(af - ar))
    eps = np.logspace(-4, -1, 7)
    errors = []
    for value in eps:
        ar = np.array(red.scalar(1.0, float(value), P, True))
        af = np.array(full.scalar(1.0, float(value), P, True))
        errors.append(float(np.linalg.norm(af - ar)))
    order = float(np.polyfit(np.log(eps), np.log(errors), 1)[0])
    return (rows, dict(rigid_max_abs_error=max(rigid_err), small_qdot_error_order=order, eps=eps.tolist(), absolute_errors=errors))

upward_zero.direction = 1
upward_zero.terminal = False

def dimensional_metrics(P, result, omega_n=None):
    """Complete bilateral-system scales; never confuse code half inertia.

    Input work, power, and force metrics are nondimensional. The alternating
    generalized force is not individual muscle tensile force or stress.
    """
    from parameters import GAMMA
    wn = BEE['wn'] if omega_n is None else float(omega_n)
    inertia = 2.0*P['I_theta']
    force_scale = inertia*THETA*wn**2/GAMMA
    power_scale = inertia*THETA**2*wn**3
    work_scale = inertia*THETA**2*wn**2
    return dict(I_full_kg_m2=inertia,force_scale_N=force_scale,
                power_scale_W=power_scale,work_scale_J=work_scale,
                peak_effective_force_N=result['peak_muscle_force']*force_scale,
                peak_positive_power_W=result['peak_positive_muscle_power']*power_scale,
                mean_positive_power_W=result['W_plus']/result['period']*power_scale,
                positive_work_J=result['W_plus']*work_scale)
