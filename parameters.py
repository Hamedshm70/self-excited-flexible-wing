"""Physical parameters and exact per-figure protocols for the corrected manuscript.
All dimensional constants use SI units. Search seeds affect bracketing only;
reference results never replace an integration. See reference/provenance.json.
"""
from __future__ import annotations
import numpy as np
def muscle_a2_a3(t_o, kappa, omega_n):
    """dSA second-order filter coefficients, non-dimensionalised by omega_n.
    Corrected sum-of-rates damping; coefficients evaluated independently.
    The selected t_o=5.83 ms is a computational parameterization.
    """
    r3 = np.log(1.0 / kappa) / (t_o * (1.0 - kappa))
    r4 = kappa * r3
    alpha2 = r3 + r4
    alpha3 = r3 * r4
    return (alpha2 / omega_n, alpha3 / omega_n ** 2)


_fn = 94.9
_wn = 2*np.pi*_fn
T0 = 5.83e-3
KAPPA = 0.17
_a2, _a3 = muscle_a2_a3(T0, KAPPA, _wn)
BEE = dict(L=1.09e-2, c=3.20e-3, mw=3.44e-7, Theta=1.10,
           lambda_t=0.095, fn=_fn, wn=_wn, a2=_a2, a3=_a3,
           CD=3.0, rho=1.25, zeta=0.05, kappa_a=0.5, added_mass=True)
GAMMA = 3.82e-5
SLOPE_LIMIT = 0.15
SLOPE_FACTOR = 1.3765055104752342
METHOD = "DOP853"
RTOL = 3e-9
ATOL = 3e-11
MAX_STEP = 0.06
NSPAN = 64
T_END = 600.0
TRANSIENT = 0.75
N_CYCLE = 6001
Y0 = np.array([0.02, 0., 0., 0., 0., 0.])
THETA = float(BEE['Theta'])

# Figure-specific production settings. Lowercase prefixes match analysis methods.
f6_SLOPE_LIMIT = 0.15
f6_RW_A = [1.4, 2.0, 2.5, 3.0, 5.5, 7.5]
f6_MU_A = np.arange(2.0, 90.0 + 1e-12, 2.0)
f6_T0S_MS = [3.0, 4.4, 5.83, 8.0]
f6_RW_B = np.round(np.linspace(1.3, 5.5, 15), 1).tolist()
f6_NSPAN = 64
f6_SCALARS = ('Omega', 'stroke_amp', 'q_max', 'max_slope', 'rotation_deg', 'phase_q_minus_theta_deg', 'periodic_state_rel', 'cycle_amp_rel', 'cycle_work_rel', 'solver_nfev')
f7_RW = [round(float(value), 3) for value in np.linspace(1.35, 5.5, 18)]
f7_NWS = [0.3, 0.84, 1.3, 3.0]
f7_TARGET = 0.22
_wet_mass = BEE['mw'] + (np.pi/4.0)*BEE['rho']*BEE['c']**2*BEE['L'] if BEE['added_mass'] else BEE['mw']
_half_inertia = (1.0+BEE['lambda_t'])*(_wet_mass*BEE['L']**2/3.0)
_third_area_moment = BEE['c']*BEE['L']**4/4.0
f7_NW_REFERENCE = float(_half_inertia/(BEE['kappa_a']*BEE['CD']*BEE['rho']*_third_area_moment*BEE['Theta']))
f8_TARGETS = (0.2, 0.5, 0.8, 1.0)
f8_EXTRA_RWS = {0.2: (1.95, 2.0), 0.5: (3.6, 3.65), 0.8: (5.4, 5.45), 1.0: ()}
f8_H_FRACS = (0.04, 0.02, 0.01, 0.005)
f8_SELECTED_H = 0.01
f8_SLOPE_LIMIT = 0.15
f8_MARGINAL_LIMIT = 0.2
f8_SLOPE_FACTOR = SLOPE_FACTOR
f8_RTOL = 3e-09
f8_ATOL = 3e-11
f8_MAX_STEP = 0.06
f8_NSPAN = 64
f8_BASE_DURATION = 360.0
f8_WARM_DURATION = 150.0
f8_CYCLE_SAMPLES = 1601
f8_MAX_SETTLE_BLOCKS = 3
f9_TARGET = 0.2
f9_T0_MEASURED = T0  # selected computational baseline, not a new measurement
f9_T0_GRID = tuple(sorted(set(np.geomspace(0.0024, 0.009, 12).tolist() + [f9_T0_MEASURED])))
f9_H_FRACS = (0.04, 0.02, 0.01, 0.005)
f9_SELECTED_H = 0.01
f9_SLOPE_LIMIT = 0.15
f9_SLOPE_FACTOR = SLOPE_FACTOR
f9_RTOL = 3e-09
f9_ATOL = 3e-11
f9_MAX_STEP = 0.06
f9_NSPAN = 64
f9_CYCLE_SAMPLES = 1601
f10_TARGET = 0.2
f10_NSPAN = 64
f10_METRICS = ('Omega', 'stroke_amp', 'q_max', 'max_slope', 'rotation_deg', 'phase_q_minus_theta_deg', 'W_plus', 'W_minus', 'W_m_net', 'W_aero', 'W_struct', 'W_trunc', 'peak_muscle_force', 'peak_muscle_power', 'peak_aero_stroke_load', 'peak_aero_power')
f11_TARGET = 0.2
f11_SLOPE_LIMIT = 0.15
