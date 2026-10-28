"""Corrected scientific kernels and self-contained production workflows.
No reference result is used as a computed answer. See REPRODUCIBILITY_REPORT.md for execution evidence.
"""
from __future__ import annotations
import csv,json,math,time,hashlib
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp,simpson
from scipy.optimize import brentq,minimize_scalar
import model as n5
from model import Aerodynamics,dimensionless_groups,solve_cycle
from parameters import *
PACKAGE=Path(__file__).resolve().parent
OUTPUT=PACKAGE/'outputs'
REFERENCE=PACKAGE/'reference'
def f6_groups(t0_ms: float, r_w: float) -> dict:
    species = dict(n5.BEE)
    (species['a2'], species['a3']) = muscle_a2_a3(t0_ms / 1000.0, 0.17, species['wn'])
    return n5.dimensionless_groups(species, r_w)

def f6_scalar(result: dict) -> dict:
    return {key: int(result[key]) if key == 'solver_nfev' else float(result[key]) for key in f6_SCALARS}

def f6_boundary_task(task: tuple[float, float, float]) -> dict:
    """Classify and refine the first validity/rollover/numerical endpoint."""
    (t0_ms, r_w, old_mu) = task
    production: dict[float, dict] = {}

    def evaluate(mu: float) -> dict:
        key = round(float(mu), 8)
        if key not in production:
            production[key] = f6_direct_cycle(t0_ms, r_w, float(mu))
        return production[key]
    cap = min(95.0, max(20.0, 1.6 * old_mu + 8.0))
    quick = []
    first_failure_mu = None
    for mu in np.linspace(2.0, cap, 10):
        try:
            quick.append((float(mu), f6_quick_cycle(t0_ms, r_w, float(mu))))
        except RuntimeError:
            first_failure_mu = float(mu)
            break
    if len(quick) < 3:
        raise RuntimeError(f'insufficient direct branch samples: {task}')
    crossing = None
    for (left, right) in zip(quick[:-1], quick[1:]):
        if left[1]['max_slope'] <= f6_SLOPE_LIMIT < right[1]['max_slope']:
            crossing = (left[0], right[0])
            break
    quick_mu = np.asarray([item[0] for item in quick])
    quick_stroke = np.asarray([item[1]['stroke_amp'] for item in quick])
    quick_slope = np.asarray([item[1]['max_slope'] for item in quick])
    local_peaks = [index for index in range(1, len(quick_stroke) - 1) if quick_stroke[index] >= quick_stroke[index - 1] and quick_stroke[index + 1] < 0.99 * quick_stroke[index]]
    peak_index = local_peaks[0] if local_peaks else int(np.argmax(quick_stroke))
    rollover_seen = bool(peak_index < len(quick) - 1 and np.min(quick_stroke[peak_index + 1:]) < 0.98 * quick_stroke[peak_index])
    if rollover_seen and (crossing is None or quick_mu[peak_index] < crossing[0]):
        ilo = max(0, peak_index - 1)
        ihi = min(len(quick) - 1, peak_index + 1)
        result = minimize_scalar(lambda mu: -evaluate(float(mu))['stroke_amp'], bounds=(float(quick_mu[ilo]), float(quick_mu[ihi])), method='bounded', options={'xatol': 0.002, 'maxiter': 20})
        envelope_mu = float(result.x)
        envelope = evaluate(envelope_mu)
        endpoint_type = 'stable_rollover_before_slope_boundary'
        slope_crossing_found = False
        slope_crossing_mu = None
    elif crossing is not None:
        (lo, hi) = crossing
        flo = evaluate(lo)['max_slope'] - f6_SLOPE_LIMIT
        fhi = evaluate(hi)['max_slope'] - f6_SLOPE_LIMIT
        if flo > 0 or fhi < 0:
            candidates = sorted(set(quick_mu.tolist() + np.linspace(lo, hi, 5).tolist()))
            pairs = []
            for mu in candidates:
                try:
                    pairs.append((mu, evaluate(mu)))
                except RuntimeError:
                    continue
            crossing = None
            for (left, right) in zip(pairs[:-1], pairs[1:]):
                if left[1]['max_slope'] <= f6_SLOPE_LIMIT < right[1]['max_slope']:
                    crossing = (left[0], right[0])
                    break
            if crossing is None:
                raise RuntimeError(f'production bracket changed unexpectedly: {task}')
            (lo, hi) = crossing
        root = float(brentq(lambda mu: evaluate(float(mu))['max_slope'] - f6_SLOPE_LIMIT, lo, hi, xtol=0.0002, rtol=2e-05, maxiter=18))
        envelope = evaluate(root)
        endpoint_type = 'slope_0.15_reporting_boundary'
        slope_crossing_found = True
        slope_crossing_mu = root
        envelope_mu = root
    elif first_failure_mu is not None:
        envelope_mu = float(quick_mu[-1])
        envelope = evaluate(envelope_mu)
        endpoint_type = 'numerical_termination_before_slope_boundary'
        slope_crossing_found = False
        slope_crossing_mu = None
    else:
        envelope_mu = float(quick_mu[-1])
        envelope = evaluate(envelope_mu)
        endpoint_type = 'mu_sweep_limit_before_slope_boundary'
        slope_crossing_found = False
        slope_crossing_mu = None
    diagnostic = [(mu, met) for (mu, met) in quick if mu <= envelope_mu + 1e-08]
    diagnostic.append((envelope_mu, envelope))
    diagnostic = sorted({round(mu, 8): met for (mu, met) in diagnostic}.items())
    stroke = np.asarray([met['stroke_amp'] for (_, met) in diagnostic])
    slope = np.asarray([met['max_slope'] for (_, met) in diagnostic])
    min_dstroke = float(np.min(np.diff(stroke))) if len(stroke) > 1 else math.nan
    min_dslope = float(np.min(np.diff(slope))) if len(slope) > 1 else math.nan
    monotone_stroke = bool(np.all(np.diff(stroke) >= -2e-05))
    monotone_slope = bool(np.all(np.diff(slope) >= -2e-05))
    return {'t0_ms': t0_ms, 'r_w': r_w, 'endpoint_type': endpoint_type, 'slope_crossing_found': slope_crossing_found, 'mu_slope_crossing': slope_crossing_mu, 'envelope_mu': envelope_mu, **{f'envelope_{key}': value for (key, value) in envelope.items()}, 'first_scout_failure_mu': first_failure_mu, 'scout_cap_mu': cap, 'scout_rollover_seen': rollover_seen, 'preendpoint_sample_count': len(diagnostic), 'preendpoint_monotone_stroke': monotone_stroke, 'preendpoint_monotone_slope': monotone_slope, 'minimum_preendpoint_delta_stroke': min_dstroke, 'minimum_preendpoint_delta_slope': min_dslope, 'samples': [dict(mu=mu, **met) for (mu, met) in quick]}

def f8_classify(slope: float) -> str:
    if slope <= f8_SLOPE_LIMIT:
        return 'quantitative'
    if slope <= f8_MARGINAL_LIMIT:
        return 'marginal_diagnostic'
    return 'outside_diagnostic'

def f8_cycle_metrics(sol, cycle_samples: int=f8_CYCLE_SAMPLES) -> dict:
    crossings = np.asarray(sol.t_events[0], dtype=float)
    if crossings.size < 4:
        raise RuntimeError('fewer than three complete cycles after settling')
    (ta, tb) = (float(crossings[-2]), float(crossings[-1]))
    (tc, td) = (float(crossings[-3]), float(crossings[-2]))
    t = np.linspace(ta, tb, int(cycle_samples))
    tp = np.linspace(tc, td, int(cycle_samples))
    y = sol.sol(t)
    yp = sol.sol(tp)
    amp = 0.5 * (float(np.max(y[0])) - float(np.min(y[0])))
    amp_prev = 0.5 * (float(np.max(yp[0])) - float(np.min(yp[0])))
    q_max = float(np.max(np.abs(y[2])))
    scale = max(float(np.max(np.abs(y))), 1e-12)
    state_rel = float(np.max(np.abs(sol.sol(tb) - sol.sol(ta))) / scale)
    return {'amp': amp, 'Omega': float(2.0 * np.pi / (tb - ta)), 'q_max': q_max, 'max_slope': float(f8_SLOPE_FACTOR * q_max), 'rotation_deg': float(np.degrees(np.arctan(f8_SLOPE_FACTOR * q_max))), 'periodic_state_rel': state_rel, 'cycle_amp_rel': abs(amp - amp_prev) / max(abs(amp), 1e-12), 'state': np.asarray(sol.sol(tb), dtype=float)}

def f8_settle(P: dict, mu: float, y0: np.ndarray | None, *, warm: bool, nspan: int=f8_NSPAN, rtol: float=f8_RTOL, atol: float=f8_ATOL, max_step: float=f8_MAX_STEP, cycle_samples: int=f8_CYCLE_SAMPLES) -> dict:
    """Settle a direct-aerodynamic limit cycle from a fixed or warm start."""
    aero = n5.Aerodynamics('direct', nspan)
    state = n5.Y0.copy() if y0 is None else np.asarray(y0, dtype=float).copy()
    if np.linalg.norm(state[:4]) < 1e-11:
        state[0] += 0.02
    duration = f8_WARM_DURATION if warm else f8_BASE_DURATION
    previous_amp = None
    accepted = False
    result = None
    total_nfev = 0
    for block in range(f8_MAX_SETTLE_BLOCKS):
        sol = solve_ivp(n5.rhs, (0.0, duration), state, args=(P, float(mu), True, aero), method='DOP853', rtol=rtol, atol=atol, max_step=max_step, dense_output=True, events=n5.upward_zero)
        total_nfev += int(sol.nfev)
        if not sol.success:
            raise RuntimeError(sol.message)
        result = f8_cycle_metrics(sol, cycle_samples=cycle_samples)
        state = result['state']
        if previous_amp is not None and abs(result['amp'] - previous_amp) <= 1e-06 * max(result['amp'], 0.001) and (result['cycle_amp_rel'] <= 3e-06) and (result['periodic_state_rel'] <= 3e-05):
            accepted = True
            break
        previous_amp = result['amp']
        duration = 100.0
    if not accepted:
        raise RuntimeError('unconverged amplitude cycle after settling blocks')
    result['mu'] = float(mu)
    result['settle_blocks'] = int(block + 1)
    result['solver_nfev'] = int(total_nfev)
    return result

def f8_refine_target(P: dict, target: float, mu_guess: float, **settings) -> dict:
    """Newton refinement on the connected positive-slope branch."""
    mu = float(mu_guess)
    base = f8_settle(P, mu, None, warm=False, **settings)
    for _ in range(7):
        error = base['amp'] - target
        if abs(error) <= 1e-8:
            break
        h = max(0.01 * mu, 0.005)
        minus = f8_settle(P, mu - h, base['state'], warm=True, **settings)
        plus = f8_settle(P, mu + h, base['state'], warm=True, **settings)
        derivative = (plus['amp'] - minus['amp']) / (2.0 * h)
        if not np.isfinite(derivative) or derivative <= 1e-07:
            raise RuntimeError('target refinement left the positive-slope branch')
        step = float(np.clip(error / derivative, -0.2 * mu, 0.2 * mu))
        mu_new = mu - step
        if mu_new <= 0.0:
            raise RuntimeError('nonpositive drive proposed during target refinement')
        seed = minus['state'] if mu_new < mu else plus['state']
        mu = float(mu_new)
        base = f8_settle(P, mu, seed, warm=True, **settings)
    if abs(base['amp'] - target) > 3e-8:
        raise RuntimeError(f"target mismatch: target={target}, achieved={base['amp']}")
    return base

def f8_scout_target(P: dict, target: float, mu_guess: float, **settings) -> tuple[dict | None, dict]:
    """Follow the connected branch when local refinement from the old drive fails."""
    mu_grid = np.linspace(max(0.25, 0.18 * mu_guess), max(1.55 * mu_guess, mu_guess + 4.0), 18)
    branch = []
    state = None
    first_failure_mu = None
    first_failure_message = None
    for (index, mu) in enumerate(mu_grid):
        try:
            item = f8_settle(P, float(mu), state, warm=index > 0, **settings)
        except RuntimeError as exc:
            first_failure_mu = float(mu)
            first_failure_message = str(exc)
            break
        state = item['state']
        branch.append(item)
        if len(branch) >= 2:
            (left, right) = (branch[-2], branch[-1])
            if left['amp'] <= target <= right['amp'] and right['amp'] > left['amp']:
                (lo, hi) = (left, right)
                for _ in range(14):
                    mu_mid = 0.5 * (lo['mu'] + hi['mu'])
                    seed = lo['state'] if mu_mid - lo['mu'] <= hi['mu'] - mu_mid else hi['state']
                    mid = f8_settle(P, mu_mid, seed, warm=True, **settings)
                    if mid['amp'] < target:
                        lo = mid
                    else:
                        hi = mid
                    if abs(mid['amp'] - target) <= max(3e-05, 0.0001 * target):
                        return (mid, {'scout_max_amp': max((row['amp'] for row in branch)), 'scout_max_amp_mu': max(branch, key=lambda row: row['amp'])['mu'], 'scout_points': len(branch), 'scout_found_crossing': True})
                closest = min((lo, hi), key=lambda row: abs(row['amp'] - target))
                return (closest, {'scout_max_amp': max((row['amp'] for row in branch)), 'scout_max_amp_mu': max(branch, key=lambda row: row['amp'])['mu'], 'scout_points': len(branch), 'scout_found_crossing': True})
    if not branch:
        raise RuntimeError(first_failure_message or 'branch scout produced no solution')
    best = max(branch, key=lambda row: row['amp'])
    return (None, {'scout_max_amp': best['amp'], 'scout_max_amp_mu': best['mu'], 'scout_max_slope': best['max_slope'], 'scout_points': len(branch), 'scout_found_crossing': False, 'scout_amp_at_cap': branch[-1]['amp'], 'scout_first_failure_mu': first_failure_mu, 'scout_first_failure_message': first_failure_message})

def f9_species_t0(t0: float) -> dict:
    species = dict(n5.BEE)
    (species['a2'], species['a3']) = muscle_a2_a3(t0, 0.17, n5.BEE['wn'])
    return species

def f9_cycle_metrics(sol, cycle_samples=f9_CYCLE_SAMPLES) -> dict:
    crossings = np.asarray(sol.t_events[0], dtype=float)
    if crossings.size < 4:
        raise RuntimeError('fewer than three complete cycles')
    (ta, tb) = (float(crossings[-2]), float(crossings[-1]))
    (tc, td) = (float(crossings[-3]), float(crossings[-2]))
    t = np.linspace(ta, tb, cycle_samples)
    tp = np.linspace(tc, td, cycle_samples)
    (y, yp) = (sol.sol(t), sol.sol(tp))
    amp = 0.5 * (float(np.max(y[0])) - float(np.min(y[0])))
    amp_prev = 0.5 * (float(np.max(yp[0])) - float(np.min(yp[0])))
    omega = 2 * np.pi / (tb - ta)
    omega_prev = 2 * np.pi / (td - tc)
    qmax = float(np.max(np.abs(y[2])))
    scale = max(float(np.max(np.abs(y))), 1e-12)
    return {'amp': amp, 'Omega': float(omega), 'q_max': qmax, 'max_slope': float(f9_SLOPE_FACTOR * qmax), 'rotation_deg': float(np.degrees(np.arctan(f9_SLOPE_FACTOR * qmax))), 'periodic_state_rel': float(np.max(np.abs(sol.sol(tb) - sol.sol(ta))) / scale), 'cycle_amp_rel': abs(amp - amp_prev) / max(abs(amp), 1e-12), 'cycle_frequency_rel': abs(omega - omega_prev) / max(abs(omega), 1e-12), 'state': np.asarray(sol.sol(tb), dtype=float)}

def f9_settle(P: dict, mu: float, flexible: bool, y0=None, *, warm=False, nspan=f9_NSPAN, rtol=f9_RTOL, atol=f9_ATOL, max_step=f9_MAX_STEP, cycle_samples=f9_CYCLE_SAMPLES) -> dict:
    aero = n5.Aerodynamics('direct', nspan)
    state = n5.Y0.copy() if y0 is None else np.asarray(y0, dtype=float).copy()
    if np.linalg.norm(state[:4]) < 1e-11:
        state[0] += 0.02
    duration = 150.0 if warm else 360.0
    previous_amp = None
    accepted = False
    result = None
    nfev = 0
    for block in range(3):
        sol = solve_ivp(n5.rhs, (0.0, duration), state, args=(P, float(mu), flexible, aero), method='DOP853', rtol=rtol, atol=atol, max_step=max_step, dense_output=True, events=n5.upward_zero)
        nfev += int(sol.nfev)
        if not sol.success:
            raise RuntimeError(sol.message)
        result = f9_cycle_metrics(sol, cycle_samples)
        state = result['state']
        if previous_amp is not None and abs(result['amp'] - previous_amp) <= 1e-06 * max(result['amp'], 0.001) and (result['cycle_amp_rel'] <= 3e-06) and (result['cycle_frequency_rel'] <= 3e-06) and (result['periodic_state_rel'] <= 3e-05):
            accepted = True
            break
        previous_amp = result['amp']
        duration = 100.0
    if not accepted:
        raise RuntimeError('unconverged timing cycle after settling blocks')
    result.update(mu=float(mu), settle_blocks=block + 1, solver_nfev=nfev)
    return result

def f9_refine_target(P: dict, flexible: bool, mu_guess: float, seed=None, **settings) -> dict:
    mu = float(mu_guess)
    base = f9_settle(P, mu, flexible, seed, warm=seed is not None, **settings)
    for _ in range(8):
        error = base['amp'] - f9_TARGET
        if abs(error) <= 1e-8:
            break
        h = max(0.01 * mu, 0.005)
        minus = f9_settle(P, mu - h, flexible, base['state'], warm=True, **settings)
        plus = f9_settle(P, mu + h, flexible, base['state'], warm=True, **settings)
        derivative = (plus['amp'] - minus['amp']) / (2 * h)
        if not np.isfinite(derivative) or derivative <= 1e-07:
            raise RuntimeError('target refinement left rising branch')
        step = float(np.clip(error / derivative, -0.2 * mu, 0.2 * mu))
        mu_new = mu - step
        if mu_new <= 0:
            raise RuntimeError('nonpositive drive')
        base = f9_settle(P, mu_new, flexible, minus['state'] if mu_new < mu else plus['state'], warm=True, **settings)
        mu = mu_new
    if abs(base['amp'] - f9_TARGET) > 3e-8:
        raise RuntimeError('target mismatch')
    return base

def f9_scout_target(P: dict, flexible: bool, mu_guess: float, **settings):
    grid = np.linspace(max(0.25, 0.2 * mu_guess), max(1.8 * mu_guess, mu_guess + 6), 22)
    (branch, state) = ([], None)
    for (index, mu) in enumerate(grid):
        try:
            item = f9_settle(P, float(mu), flexible, state, warm=index > 0, **settings)
        except RuntimeError:
            break
        state = item['state']
        branch.append(item)
        if len(branch) >= 2 and branch[-2]['amp'] <= f9_TARGET <= branch[-1]['amp']:
            (lo, hi) = (branch[-2], branch[-1])
            for _ in range(15):
                mid_mu = 0.5 * (lo['mu'] + hi['mu'])
                mid = f9_settle(P, mid_mu, flexible, lo['state'], warm=True, **settings)
                if mid['amp'] < f9_TARGET:
                    lo = mid
                else:
                    hi = mid
                if abs(mid['amp'] - f9_TARGET) < 3e-05:
                    return (mid, {'scout_max_amp': max((x['amp'] for x in branch))})
            return (min((lo, hi), key=lambda x: abs(x['amp'] - f9_TARGET)), {'scout_max_amp': max((x['amp'] for x in branch))})
    return (None, {'scout_max_amp': max((x['amp'] for x in branch), default=math.nan), 'scout_max_mu': max((x['mu'] for x in branch), default=math.nan)})

def f9_matched(t0: float, r_w: float, flexible: bool, mu_guess: float, seed=None, **settings) -> tuple[dict | None, dict]:
    P = n5.dimensionless_groups(f9_species_t0(t0), r_w)
    try:
        return (f9_refine_target(P, flexible, mu_guess, seed, **settings), {})
    except RuntimeError:
        base, info = f9_scout_target(P, flexible, mu_guess, **settings)
        if base is None: return None, info
        return f9_refine_target(P, flexible, base['mu'], **settings), info

def f11_groups_with_nw(rw: float, nw: float) -> dict:
    """Use the exact historical consistent-N_W scaling."""
    p = n5.dimensionless_groups(n5.BEE, float(rw))
    ref = n5.dimensionless_groups(n5.BEE, 3.0)
    scale = ref['Nw'] / float(nw)
    p['Nw'] = float(nw)
    for key in ('eta_A', 'eta_B', 'eta_C'):
        p[key] = ref[key] * scale
    return p

def f11_amplitude(p, mu, nspan=32):
    try:
        r = n5.solve_cycle(p, float(mu), True, 'direct', nspan=nspan, t_end=260.0, rtol=2e-08, atol=2e-10, max_step=0.1, n_cycle=801)
        return float(r['stroke_amp'])
    except Exception:
        return float('nan')

def f11_tune_first_branch(p, initial):
    """Locate the first rising target crossing and diagnose non-reachability."""
    top = min(200.0, max(18.0, 2.0 * float(initial)))
    lower = max(0.25, 0.25 * float(initial))
    extension = np.geomspace(1.01 * top, 200.0, 16) if top < 198.0 else np.array([])
    grid = np.unique(np.r_[np.geomspace(lower, top, 25), np.linspace(lower, top, 21), extension])
    grid.sort()
    cache = {}
    sampled = []

    def f(mu):
        key = round(float(mu), 9)
        if key not in cache:
            cache[key] = f11_amplitude(p, float(mu))
        return cache[key] - f11_TARGET
    previous = None
    peak = (-np.inf, np.nan)
    declining = 0
    for mu in grid:
        value = f(mu)
        amp = value + f11_TARGET
        if not np.isfinite(amp):
            continue
        sampled.append((float(mu), float(amp)))
        if amp > peak[0]:
            peak = (float(amp), float(mu))
            declining = 0
        elif amp < 0.995 * peak[0]:
            declining += 1
        if previous is not None and previous[1] <= 0.0 <= value:
            root = float(brentq(f, previous[0], float(mu), xtol=0.0001, rtol=1e-05, maxiter=30))
            return (root, sampled, 'target_reached_first_rising_branch')
        previous = (float(mu), float(value))
        if peak[0] > 0.05 and amp < 0.75 * peak[0] and (declining >= 2):
            break
    if sampled:
        imax = int(np.argmax([a for (_, a) in sampled]))
        lo = sampled[max(0, imax - 1)][0]
        hi = sampled[min(len(sampled) - 1, imax + 1)][0]
        if hi > lo:
            optimum = minimize_scalar(lambda x: -(f(float(x)) + f11_TARGET), bounds=(lo, hi), method='bounded', options={'xatol': 0.002, 'maxiter': 18})
            peak_mu = float(optimum.x)
            peak_amp = float(-optimum.fun)
            sampled.append((peak_mu, peak_amp))
            sampled.sort()
            if peak_amp >= f11_TARGET:
                below = [(m, a) for (m, a) in sampled if m < peak_mu and a < f11_TARGET]
                if below:
                    lo_root = max(below, key=lambda x: x[0])[0]
                    root = float(brentq(f, lo_root, peak_mu, xtol=0.0001, rtol=1e-05, maxiter=30))
                    return (root, sampled, 'target_reached_near_rollover')
    status = 'rollover_below_target' if declining >= 2 else 'imposed_mu_scan_limit_below_target'
    return (None, sampled, status)

def loss_mu_from_frequency_structural(P: dict, omega: float, flexible: bool, gamma_theta: float) -> complex:
    """Complex drive implied at s=i*omega with thoracic structural loss."""
    s = 1j * omega
    (a2, a3) = (P['a2'], P['a3'])
    mech = s * s + 1.0 + 1j * gamma_theta
    if flexible:
        bend = s * s + 2.0 * P['zeta'] * P['r'] * s + P['r'] ** 2
        if abs(bend) < 1e-12:
            return complex(np.nan, np.nan)
        mech -= P['b_tq'] * P['g_qt'] * s ** 4 / bend
    return -mech * (s * s + a2 * s + a3) / (a3 * s)

def loss_structural_threshold(P: dict, flexible: bool, gamma_theta: float, n_scan: int) -> tuple[float, float, float]:
    """Lowest positive-drive harmonic onset for structural thorax loss."""
    if gamma_theta == 0.0 and (not flexible):
        return (0.0, 1.0, 0.0)
    omega_grid = np.linspace(0.03, 3.0, n_scan)
    imag = np.array([loss_mu_from_frequency_structural(P, w, flexible, gamma_theta).imag for w in omega_grid])
    candidates: list[tuple[float, float, float]] = []
    for i in range(n_scan - 1):
        if not np.isfinite(imag[i]) or not np.isfinite(imag[i + 1]):
            continue
        roots: list[float]
        if imag[i] == 0.0:
            roots = [float(omega_grid[i])]
        elif imag[i] * imag[i + 1] < 0.0:
            roots = [brentq(lambda w: loss_mu_from_frequency_structural(P, w, flexible, gamma_theta).imag, omega_grid[i], omega_grid[i + 1], xtol=2e-13, rtol=2e-13)]
        else:
            roots = []
        for omega in roots:
            mu = loss_mu_from_frequency_structural(P, omega, flexible, gamma_theta)
            if mu.real >= -1e-09:
                candidates.append((max(0.0, float(mu.real)), float(omega), abs(float(mu.imag))))
    if not candidates:
        raise RuntimeError('no positive structural-loss characteristic root')
    return min(candidates, key=lambda item: item[0])

def onset_jacobian(P: dict, mu: float, flexible: bool) -> np.ndarray:
    """Exact rest-state Jacobian; quadratic aerodynamic terms vanish."""
    (a2, a3) = (P['a2'], P['a3'])
    zt = P.get('zeta_theta', 0.0)
    if not flexible:
        J = np.zeros((4, 4))
        J[0, 1] = 1.0
        J[1, 0] = -1.0
        J[1, 1] = -2.0 * zt
        J[1, 2] = mu
        J[2, 3] = 1.0
        J[3, 1] = -a3
        J[3, 2] = -a3
        J[3, 3] = -a2
        return J
    (b, g, r) = (P['b_tq'], P['g_qt'], P['r'])
    zw = P['zeta']
    det = 1.0 - b * g
    J = np.zeros((6, 6))
    J[0, 1] = 1.0
    J[2, 3] = 1.0
    J[4, 5] = 1.0
    J[1, 0] = -1.0 / det
    J[1, 1] = -2.0 * zt / det
    J[1, 2] = b * r * r / det
    J[1, 3] = 2.0 * b * zw * r / det
    J[1, 4] = mu / det
    J[3, 0] = g / det
    J[3, 1] = 2.0 * g * zt / det
    J[3, 2] = -r * r / det
    J[3, 3] = -2.0 * zw * r / det
    J[3, 4] = -g * mu / det
    J[5, 1] = -a3
    J[5, 4] = -a3
    J[5, 5] = -a2
    return J

def onset_spectral_abscissa(P: dict, mu: float, flexible: bool) -> float:
    return float(np.max(np.linalg.eigvals(onset_jacobian(P, mu, flexible)).real))

def onset_crossing_eigenvalue(P: dict, mu: float, flexible: bool) -> complex:
    vals = np.linalg.eigvals(onset_jacobian(P, mu, flexible))
    positive = vals[vals.imag > 1e-07]
    if positive.size == 0:
        return complex(vals[np.argmax(vals.real)])
    return complex(positive[np.argmax(positive.real)])

def onset_eigen_threshold(P: dict, flexible: bool) -> tuple[float, float, float]:
    """Numerical Hopf threshold from the full rest-state eigenproblem."""
    alpha0 = onset_spectral_abscissa(P, 0.0, flexible)
    if abs(alpha0) <= 2e-10:
        eig = onset_crossing_eigenvalue(P, 0.0, flexible)
        return (0.0, abs(eig.imag), alpha0)
    if alpha0 > 0.0:
        raise RuntimeError(f'rest already unstable at zero drive: alpha={alpha0}')
    mus = np.r_[0.0, np.geomspace(1e-08, 200.0, 1000)]
    (previous_mu, previous_alpha) = (mus[0], alpha0)
    for mu in mus[1:]:
        alpha = onset_spectral_abscissa(P, float(mu), flexible)
        if alpha >= 0.0:
            root = brentq(lambda value: onset_spectral_abscissa(P, value, flexible), previous_mu, float(mu), xtol=2e-13, rtol=2e-13)
            eig = onset_crossing_eigenvalue(P, root, flexible)
            return (float(root), abs(eig.imag), alpha0)
        (previous_mu, previous_alpha) = (float(mu), alpha)
    raise RuntimeError('no onset below mu=200')

def limiting_cases():
    (B1, S1) = (1.8751040687, 0.7340955)
    xg = np.linspace(0, 1, 4001)
    raw = np.cosh(B1 * xg) - np.cos(B1 * xg) - S1 * (np.sinh(B1 * xg) - np.sin(B1 * xg))
    Ph = raw / (np.cosh(B1) - np.cos(B1) - S1 * (np.sinh(B1) - np.sin(B1)))
    j0 = np.trapz(Ph ** 2, xg)
    j1 = np.trapz(xg * Ph, xg)
    j2 = np.trapz(xg ** 2 * Ph, xg)
    j3 = np.trapz(xg * Ph ** 2, xg)
    (mw, L, c, m, k) = (5e-05, 0.05, 0.018, 0.0007, 2850.0)
    (Gam, ww1, CD, rho) = (0.0008, 2 * np.pi * 60.0, 3.0, 1.25)
    IO = mw * L ** 2 / 3.0
    (Ith, Mtq, Mq) = (m * Gam ** 2 + 2 * IO, 2 * mw * L * j1, 2 * mw * j0)
    (Kth, Kq) = (k * Gam ** 2, 2 * mw * j0 * ww1 ** 2)
    (CA, CB, CC) = (CD * rho * c * L ** 4 / 4.0, CD * rho * c * L ** 3 * j2, 2 * CD * rho * c * L ** 2 * j3)
    import numpy.linalg as la

    def Ein(F0, f, flexible, zeta=0.05):
        wf = 2 * np.pi * f
        ncyc = 45
        t_end = ncyc / f
        Cq = 2 * zeta * ww1 * Mq

        def rhs(t, y):
            (th, thd, q, qd) = y
            tau = Gam * F0 * np.sin(wf * t)
            if flexible:
                rth = tau - Kth * th - CA * thd * abs(thd)
                rq = -Cq * qd - (Kq - Mq * thd ** 2) * q - CB * thd * abs(thd) - CC * abs(thd) * qd
                (thdd, qdd) = la.solve(np.array([[Ith, Mtq], [Mtq, Mq]]), np.array([rth, rq]))
            else:
                thdd = (tau - Kth * th - CA * thd * abs(thd)) / Ith
                qd = qdd = 0.0
            return [thd, thdd, qd, qdd]
        tev = np.linspace(0, t_end, ncyc * 90)
        sol = solve_ivp(rhs, (0, t_end), [0, 0, 0, 0], t_eval=tev, rtol=1e-07, atol=1e-10)
        m2 = tev > 0.6 * t_end
        thd = sol.y[1][m2]
        tt = tev[m2]
        Pw = F0 * np.sin(wf * tt) * (Gam * thd)
        return np.trapz(np.abs(Pw), tt)
    ff = np.linspace(15, 35, 7)
    ratios = {F0: np.array([Ein(F0, f, True) / Ein(F0, f, False) for f in ff]) for F0 in (1.0, 1.5, 2.0)}
    fns = np.linspace(45, 175, 10)
    fwb = []
    for fn in fns:
        wn = 2 * np.pi * fn
        (a2, a3) = muscle_a2_a3(T0, KAPPA, wn)
        sp = dict(BEE)
        sp.update(fn=fn, wn=wn, a2=a2, a3=a3)
        R, _ = cycle(dimensionless_groups(sp, 6.0), 6.0, False,omega_n=wn)
        fwb.append(R['Omega'] * fn if np.isfinite(R['Omega']) else np.nan)
    fwb = np.array(fwb)
    return (ff, ratios, fns, fwb)

def hopf_threshold(P):
    """Closed-form self-excitation threshold mu_H (smallest positive root)."""
    (b, g, r, z) = (P['b_tq'], P['g_qt'], P['r'], P['zeta'])
    (a2, a3, bg) = (P['a2'], P['a3'], P['b_tq'] * P['g_qt'])
    A = lambda w: (1 - w ** 2) * (r ** 2 - w ** 2) - bg * w ** 4
    Bs = lambda w: 2 * z * r * w * (1 - w ** 2)
    C = lambda w: a3 - w ** 2
    D = lambda w: a2 * w
    F = lambda w: (r ** 2 - w ** 2) * (A(w) * C(w) - Bs(w) * D(w)) + 2 * z * r * w * (A(w) * D(w) + Bs(w) * C(w))
    ws = np.linspace(0.2, 1.4, 2000)
    fv = np.array([F(w) for w in ws])
    best = np.nan
    for i in range(len(ws) - 1):
        if fv[i] * fv[i + 1] < 0:
            wr = brentq(F, ws[i], ws[i + 1])
            mu = (A(wr) * C(wr) - Bs(wr) * D(wr)) / (2 * a3 * z * r * wr ** 2)
            if mu > 0 and (np.isnan(best) or mu < best):
                best = mu
    return best

# Corrected orchestration. Thresholds are the Phase 2–4 acceptance thresholds;
# regression tolerances are distinct from these numerical acceptance checks.
class NumericalFailure(RuntimeError):pass
class RolloverBelowTarget(NumericalFailure):
    def __init__(self,mu,amplitude,point):
        super().__init__(f'Converged first-branch local peak {amplitude:.12g} below target at gain {mu:.12g}')
        self.mu=mu;self.amplitude=amplitude;self.point=point

def json_safe(x):
    if isinstance(x,dict):return {str(k):json_safe(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [json_safe(v) for v in x]
    if isinstance(x,np.ndarray):return json_safe(x.tolist())
    if isinstance(x,np.generic):return json_safe(x.item())
    if isinstance(x,float) and not np.isfinite(x):return None
    return x

def read_json(path):return json.loads(Path(path).read_text())
def write_json(path,data):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(json_safe(data),indent=2,allow_nan=False)+'\n')
def write_csv(path,rows):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    fields=list(dict.fromkeys(k for r in rows for k,v in r.items() if not isinstance(v,(dict,list))))
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fields,extrasaction='ignore');w.writeheader()
        for r in rows:w.writerow({k:r.get(k) for k in fields})
def scalar_result(r):return {k:float(v) for k,v in r.items() if isinstance(v,(int,float,np.number)) and np.isfinite(v)}
def required(r,fields):
    for k in fields:
        if k not in r or r[k] is None or isinstance(r[k],(bool,np.bool_)) or not isinstance(r[k],(int,float,np.number)) or not np.isfinite(r[k]):raise NumericalFailure('Missing/invalid required metric: '+k)
def accepted(r):
    if 'shooting_residual' in r:
        required(r,['shooting_residual','periodic_state_rel','energy_balance_rel','max_slope'])
        return r['shooting_residual']<1e-7 and r['periodic_state_rel']<3e-7 and r['energy_balance_rel']<1e-7
    required(r,['periodic_state_rel','cycle_amp_rel','cycle_work_rel','energy_balance_rel','max_slope'])
    return r['periodic_state_rel']<=3e-7 and r['cycle_amp_rel']<=5e-6 and r['cycle_work_rel']<=3e-7 and r['energy_balance_rel']<1e-7

def enrich(r,P,mu,flexible,omega_n=None,rhs_function=None,restored=False):
    out=scalar_result(r);t=r['t'];v=r['theta_dot'];q=r['q'];qd=r['q_dot'];th=r['theta']
    chi=P['b_tq']/P['g_qt'] if flexible else 0.;b=P['b_tq'] if flexible else 0.
    E=.5*(v*v+th*th+2*b*v*qd+chi*(qd*qd+P['r']**2*q*q+P['sig_cf']*v*v*q*q))
    delta=float(E[-1]-E[0]);wt=0. if restored else r['W_trunc']
    residual=r['W_m_net']-delta-r['W_aero']-r['W_struct']-wt
    aero=Aerodynamics('direct',64)
    fun=rhs_function or (lambda ti,y:n5.rhs(ti,y,P,mu,flexible,aero))
    states=r['solution'].sol(t);dy=np.column_stack([fun(ti,y) for ti,y in zip(t,states.T)])
    pf=(b*dy[3]+r['aero_stroke_load']-v*np.abs(v)/P['Nw'])*v if flexible else np.zeros_like(v)
    out.update(mu_hat=float(mu),r_w=float(P['r']),W_trunc=wt,Delta_E=delta,energy_balance_residual=float(residual),energy_balance_rel=float(abs(residual)/max(abs(r['W_m_net']),1e-14)),phi_net=r['W_m_net']/r['W_plus'],phi_aero=r['W_aero']/r['W_plus'],Pi_flex=float(max(abs(pf))/(r['W_aero']/r['period'])),**n5.dimensional_metrics(P,r,omega_n))
    wn=BEE['wn'] if omega_n is None else omega_n
    out.update(dimensional_mu_N_s=mu*2*P['I_theta']*wn/GAMMA,fwb_Hz=out['Omega']*wn/(2*np.pi))
    out['periodic']=accepted(out);out['quantitative']=bool(out['periodic'] and out['max_slope']<=SLOPE_LIMIT)
    out['status']='computed';return out

def _cycle_integrate(P,mu,flexible=True,*,fine=False,y0=None,omega_n=None):
    opts=dict(nspan=128 if fine else 64,t_end=1200. if fine else 600.,rtol=3e-10 if fine else 3e-9,atol=3e-12 if fine else 3e-11,max_step=.03 if fine else .06,n_cycle=16001 if fine else 6001)
    for end in [opts['t_end'],1800. if not fine else 2400.]:
        opts['t_end']=end
        if y0 is None:r=n5.solve_cycle(P,mu,flexible,'direct',**opts)
        else:
            aero=Aerodynamics('direct',opts['nspan']);sol=solve_ivp(n5.rhs,(0,end),y0,args=(P,mu,flexible,aero),method='DOP853',rtol=opts['rtol'],atol=opts['atol'],max_step=opts['max_step'],dense_output=True,events=n5.upward_zero)
            if not sol.success:raise NumericalFailure(sol.message)
            ev=sol.t_events[0];ev=ev[ev>.75*end]
            if len(ev)<4:raise NumericalFailure('Insufficient complete cycles')
            r=n5.cycle_metrics(sol,ev[-2],ev[-1],P,mu,flexible,aero,n_cycle=opts['n_cycle']);prev=n5.cycle_metrics(sol,ev[-3],ev[-2],P,mu,flexible,aero,n_cycle=2001)
            r.update(cycle_amp_rel=abs(r['stroke_amp']-prev['stroke_amp'])/r['stroke_amp'],cycle_work_rel=abs(r['W_m_net']-prev['W_m_net'])/max(abs(r['W_m_net']),1e-14),solution=sol)
        s=enrich(r,P,mu,flexible,omega_n)
        if accepted(s):return s,r
    raise NumericalFailure('Unconverged cycle: '+str({k:s[k] for k in ['periodic_state_rel','cycle_amp_rel','cycle_work_rel','energy_balance_rel']}))

def cycle(P,mu,flexible=True,*,fine=False,y0=None,omega_n=None):
    """Resume only exact-input, exact-numerical-kernel checkpoints from this run.

    Dense output preceding the final complete cycle is discarded *after* all
    acceptance tests. The retained interpolants are unchanged, not resampled.
    """
    import inspect,pickle,copy,os
    import scipy
    kernel=hashlib.sha256()
    for path in [PACKAGE/'model.py',PACKAGE/'parameters.py']:
        kernel.update(path.read_bytes())
    for fn in [_cycle_integrate,enrich,accepted,required,scalar_result]:
        kernel.update(inspect.getsource(fn).encode())
    identity=dict(kernel=kernel.hexdigest(),numpy=np.__version__,scipy=scipy.__version__,P=P,mu=float(mu),flexible=bool(flexible),fine=bool(fine),y0=None if y0 is None else np.asarray(y0).tolist(),omega_n=omega_n)
    key=hashlib.sha256(json.dumps(identity,sort_keys=True).encode()).hexdigest()
    folder=OUTPUT/'checkpoints/cycles';folder.mkdir(parents=True,exist_ok=True);path=folder/(key+'.pickle')
    if path.exists():
        with path.open('rb') as f:entry=pickle.load(f)
        if entry['identity']!=identity or not accepted(entry['scalar']):raise NumericalFailure('Invalid numerical checkpoint')
        return entry['scalar'],entry['raw']
    start=time.time();s,r=_cycle_integrate(P,mu,flexible,fine=fine,y0=y0,omega_n=omega_n)
    raw=dict(r);sol=copy.copy(r['solution']);dense=sol.sol
    idx=max(0,int(np.searchsorted(dense.ts,r['t'][0]))-1)
    sol.sol=type(dense)(dense.ts[idx:],dense.interpolants[idx:]);sol.t=sol.t[-1:];sol.y=sol.y[:,-1:];raw['solution']=sol
    tmp=path.with_suffix('.'+str(os.getpid())+'.tmp')
    with tmp.open('wb') as f:pickle.dump(dict(identity=identity,scalar=s,raw=raw,seconds=time.time()-start),f,pickle.HIGHEST_PROTOCOL)
    tmp.replace(path)
    print('cycle',P['r'],float(mu),'fine',fine,'seconds',round(time.time()-start,2),flush=True)
    return s,raw

def save_cycle(name,r):
    path=OUTPUT/'cycles';path.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(path/(name+'.npz'),**{k:r[k] for k in ['t','theta','theta_dot','q','q_dot','muscle_state']},state=r['solution'].sol(r['t']))

def continuous_stroke_amplitude(r):
    """Physical half-amplitude from zero-velocity extrema of the dense cycle.

    This is the Phase4 rollover check's continuous-extrema measurement; sampled
    maxima remain sufficient away from the nearly flat local amplitude maximum.
    """
    t=r['t'];v=r['theta_dot'];sol=r['solution']
    roots=[brentq(lambda z:sol.sol(z)[1],t[i],t[i+1],xtol=1e-12) for i in np.where(v[:-1]*v[1:]<0)[0]]
    if len(roots)!=2:raise NumericalFailure('Expected two continuous stroke extrema')
    return float(np.ptp(sol.sol(np.asarray(roots))[0])/2)

def match(P,target,guess,flexible=True,*,fine=False):
    cache={};continuous=False
    def ev(mu):
        k=float(mu)
        if k not in cache:cache[k]=cycle(P,k,flexible,fine=fine)
        return cache[k]
    def amplitude(mu):
        s,r=ev(mu)
        return continuous_stroke_amplitude(r) if continuous else s['stroke_amp']
    if abs(ev(guess)[0]['stroke_amp']-target)<2e-8:return ev(guess)
    lo,hi=.97*guess,1.03*guess
    for _ in range(12):
        a,b=amplitude(lo),amplitude(hi)
        if a<=target<=b:break
        center=amplitude(guess)
        if center>max(a,b) and max(a,b)<target:
            continuous=True
            optimum=minimize_scalar(lambda u:-amplitude(float(u)),bounds=(lo,hi),method='bounded',options={'xatol':2e-5})
            if -optimum.fun<target:
                peak=ev(float(optimum.x))[0]
                raise RolloverBelowTarget(float(optimum.x),float(-optimum.fun),peak)
            hi=float(optimum.x);break
        lo=max(.05,.82*lo);hi*=1.18
    else:raise NumericalFailure('Matched target not bracketed on rising branch')
    class LocalTargetReached(Exception):pass
    def objective(u):
        error=amplitude(u)-target
        # The validated matching criterion is amplitude error1e-8. Near a
        # rollover, gain-coordinate precision alone chases integration noise.
        if continuous and abs(error)<=1e-8:raise LocalTargetReached(float(u))
        return error
    try:mu=brentq(objective,lo,hi,xtol=3e-8,rtol=2e-10)
    except LocalTargetReached as found:mu=found.args[0]
    s,r=ev(mu)
    if continuous:s=dict(s,sampled_stroke_amp=s['stroke_amp'],stroke_amp=amplitude(mu))
    if abs(s['stroke_amp']-target)>3e-8:raise NumericalFailure('Final matched-stroke error')
    return s,r

def f6_direct_cycle(t0_ms,r_w,mu,**options):
    # Scouts retain their original coarse integration; final cycles are gated.
    if options:
        r=n5.solve_cycle(f6_groups(t0_ms,r_w),mu,True,'direct',**dict(dict(nspan=64,t_end=600.,rtol=3e-9,atol=3e-11,max_step=.06,n_cycle=3001),**options));return f6_scalar(r)
    return f6_scalar(cycle(f6_groups(t0_ms,r_w),mu)[0])
def f6_quick_cycle(t0_ms,r_w,mu):return f6_direct_cycle(t0_ms,r_w,mu,nspan=32,t_end=300.,rtol=1e-8,atol=1e-10,max_step=.09,n_cycle=1201)

def mapped_jobs(label,tasks,worker,workers=1):
    # No stale caches: a task is recomputed unless a matching complete result is explicitly loaded for plotting.
    if workers==1:return [worker(x) for x in tasks]
    from concurrent.futures import ProcessPoolExecutor
    with ProcessPoolExecutor(max_workers=workers) as pool:return list(pool.map(worker,tasks))
def successful(rows):
    if any(x.get('execution_status')=='failed' for x in rows):raise NumericalFailure('Numerical job failed')
    return rows

def sensitivity(task,timing=False,fine=False):
    rw=task['r_w'];flex=task.get('tag')!='rigid';target=.2 if timing else task['target'];t0=task.get('t0',T0)
    species=dict(BEE);species['a2'],species['a3']=muscle_a2_a3(t0,KAPPA,BEE['wn']);P=dimensionless_groups(species,rw)
    try:base,_=match(P,target,task['mu_guess'],flex,fine=fine)
    except RuntimeError as exc:
        return dict(key=task['key'],status='unavailable',quantitative=False,r_w=rw,error=str(exc),**({'tag':task['tag'],'t0':t0} if timing else {'target':target}))
    estimates=[];allpoints=[]
    for h in [.01,.005]:
        pts=[]
        for sign in [-1,1]:
            if timing:
                sp=dict(BEE);sp['a2'],sp['a3']=muscle_a2_a3(t0*(1+sign*h),KAPPA,BEE['wn']);point,_=match(dimensionless_groups(sp,rw),target,base['mu_hat'],flex,fine=fine)
            else:point,_=cycle(P,base['mu_hat']*(1+sign*h),flex,fine=fine)
            pts.append(point);allpoints.append(point)
        a,b=pts;field='Omega' if timing else 'stroke_amp'
        estimates.append(dict(h=h,G_log=(math.log(b[field])-math.log(a[field]))/(math.log(1+h)-math.log(1-h)),G_diff=(b[field]-a[field])/(2*h*(t0 if timing else base['mu_hat'])),support_slope=max(base['max_slope'],a['max_slope'],b['max_slope'])))
    change=abs(estimates[0]['G_log']/estimates[1]['G_log']-1)
    row=dict(base,mu=base['mu_hat'],key=task['key'],G_log=estimates[0]['G_log'],G_diff=estimates[0]['G_diff'],support_slope=estimates[0]['support_slope'],fd_relative_change=change,fd_converged=change<.02,base=base,stencils=allpoints,estimates=estimates)
    if timing:row['amp']=base['stroke_amp']
    row.update(quantitative=bool(row['support_slope']<=.15 and change<.02 and all(p['periodic'] for p in allpoints)),**({'tag':task['tag'],'t0':t0} if timing else {'target':target}))
    return row

def table2_compute():
    P=dimensionless_groups(BEE,2.)
    return dict(lambda_t=BEE['lambda_t'],N_W=P['Nw'],a2=P['a2'],a3=P['a3'],beta=P['b_tq'],gamma=P['g_qt'],zeta=P['zeta'])
def table3_compute():
    rows=[]
    for g in [.095,.113]:
        rigid=loss_structural_threshold(dimensionless_groups(BEE,6.),False,g,16000)
        for rw in [1.3,1.6,2.,3.,5.5]:
            mu,w,res=loss_structural_threshold(dimensionless_groups(BEE,rw),True,g,16000)
            rows.append(dict(gamma_theta=g,r_w=rw,mu_H=mu,Omega_H=w,flexible_over_rigid=mu/rigid[0]))
    return rows

def configuration_signature():
    import scipy
    paths=sorted(PACKAGE.glob('*.py'))+sorted(REFERENCE.glob('*.json'))+[PACKAGE/'requirements.txt']
    signature={str(p.relative_to(PACKAGE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    signature.update(numpy=np.__version__,scipy=scipy.__version__)
    return signature

def check_configuration():
    signature=configuration_signature();dest=OUTPUT/'run_configuration.json'
    if dest.exists() and read_json(dest)!=signature:raise NumericalFailure('Configuration changed: use a clean output directory.')
    write_json(dest,signature)

def map_point(task):
    P=f11_groups_with_nw(task['r_w'],task['N_W'])
    try:s,r=match(P,.2,task['mu_guess'])
    except RuntimeError as exc:
        if isinstance(exc,RolloverBelowTarget):return dict(key=task['key'],N_W=task['N_W'],r_w=task['r_w'],status='unavailable',quantitative=False,search_status='converged_continuous_local_rollover_below_target',peak_amp=exc.amplitude,peak_mu=exc.mu,peak_point=exc.point)
        mu,sampled,status=f11_tune_first_branch(P,task['mu_guess'])
        if mu is None:
            # Refine the low-N_W rollover with converged cycles, not the coarse scout.
            if not sampled:return dict(key=task['key'],N_W=task['N_W'],r_w=task['r_w'],status='unavailable',quantitative=False,search_status=status)
            grid=sorted(sampled);i=max(range(len(grid)),key=lambda i:grid[i][1]);lo=grid[max(0,i-1)][0];hi=grid[min(len(grid)-1,i+1)][0]
            opt=minimize_scalar(lambda u:-cycle(P,float(u))[0]['stroke_amp'],bounds=(lo,hi),method='bounded',options={'xatol':2e-5})
            if -opt.fun<.2:return dict(key=task['key'],N_W=task['N_W'],r_w=task['r_w'],status='unavailable',quantitative=False,search_status='converged_rollover_below_target',peak_amp=-opt.fun)
            mu=brentq(lambda u:cycle(P,u)[0]['stroke_amp']-.2,lo,float(opt.x),xtol=3e-8,rtol=2e-10)
        s,r=match(P,.2,mu)
    return dict(s,key=task['key'],N_W=task['N_W'],target=.2)

def resolved_boundaries(number):
    """Recompute support roots/interior points; never interpolate metric values."""
    refs=read_json(REFERENCE/'production_values.json');seeds=read_json(REFERENCE/'search_seeds.json');out=[]
    if number==8:
        for b in refs['boundaries8']:
            target=b['target'];guess=b['point']['base']['mu_hat']
            def point(rw):return sensitivity(dict(key='boundary',r_w=rw,target=target,mu_guess=guess))
            lo=b['root_rw']-.02;hi=b['root_rw']+.02
            def f(rw):
                row=point(rw)
                if row['status']!='computed':raise NumericalFailure('Unavailable derivative boundary')
                return row['support_slope']-.15
            root=brentq(f,lo,hi,xtol=2e-6,rtol=1e-8);safe=root+2e-5;row=point(safe)
            if not row['quantitative']:raise NumericalFailure('Derivative interior point fails support')
            out.append(dict(target=target,root_rw=root,reported_rw=safe,point=row))
    elif number==11:
        for b in refs['boundaries11']:
            nw=b['N_W'];guess=b['point']['mu_hat'];lo=b['reported_rw']-.002;hi=b['reported_rw']+.002
            def point(rw):return map_point(dict(key='boundary',N_W=nw,r_w=rw,mu_guess=guess))
            for _ in range(12):
                mid=(lo+hi)/2;row=point(mid)
                if row.get('quantitative'):hi=mid
                else:lo=mid
            safe=hi+2e-5;row=point(safe)
            if not row.get('quantitative'):raise NumericalFailure('Map interior point fails support')
            # The manuscript reports a selected admissible boundary coordinate,
            # with the original bracket uncertainty, not the newly chosen
            # interior point. Near a rollover their loads differ measurably.
            selected=point(b['reported_rw'])
            if not selected.get('quantitative'):raise NumericalFailure('Published boundary coordinate fails support')
            lower=point(lo)
            out.append(dict(N_W=nw,reported_rw=b['reported_rw'],mechanism=b['mechanism'],lower_rw=lo,upper_rw=hi,point=selected,localized_reported_rw=safe,localized_point=row,lower_state=lower,reporting_scope='Selected approved coordinate; independently localized interval and interior check retained separately'))
    return out

def compute_figure(number,workers=1):
    check_configuration();tasks=read_json(REFERENCE/'search_seeds.json')[str(number)];rows=[]
    if number==2:
        ff,ratios,fns,fwb=limiting_cases()
        for task in tasks:
            if task['kind']=='forced':
                i=int(np.argmin(abs(ff-task['frequency'])));rows.append(dict(task,energy_ratio=ratios[task['F0']][i],role='diagnostic_prescribed_forcing_benchmark',energy_definition='integral_abs_actuator_power'))
            else:
                i=int(np.argmin(abs(fns-task['natural_frequency'])));fn=float(fns[i]);wn=2*np.pi*fn;sp=dict(BEE);sp.update(fn=fn,wn=wn);sp['a2'],sp['a3']=muscle_a2_a3(T0,KAPPA,wn)
                measured,_=cycle(dimensionless_groups(sp,6.),6.,False,omega_n=wn);rows.append(dict(measured,**task,wingbeat_frequency=fwb[i]))
    elif number==3:
        P=dimensionless_groups(BEE,6.);guess=tasks[0]['calibration_guess'];cache={}
        def ev(mu):
            if mu not in cache:cache[mu]=cycle(P,mu)
            return cache[mu]
        mu=brentq(lambda u:ev(u)[0]['W_aero']/ev(u)[0]['period']-.030,.8*guess,1.2*guess,xtol=2e-5,rtol=2e-6)
        P=dimensionless_groups(BEE,2.);s,r=cycle(P,mu);save_cycle('fig3',r);rows=[dict(s,key='fixed_drive')]
        aero=Aerodynamics('direct',64);sol=solve_ivp(n5.rhs,(0,170),Y0,args=(P,mu,True,aero),method='DOP853',rtol=3e-9,atol=3e-11,max_step=.06,dense_output=True)
        if not sol.success:raise NumericalFailure(sol.message)
        t=np.linspace(0,170,22000);np.savez_compressed(OUTPUT/'cycles/fig3_transient.npz',t=t,y=sol.sol(t))
    elif number==4:
        import figure4_analysis as h
        rows=h.compute(workers)
    elif number in [5,10]:
        for task in tasks:
            s,r=match(dimensionless_groups(BEE,task['r_w']),.2,task['mu_guess'],task['flexible'])
            if not s['quantitative']:raise NumericalFailure('Matched production figure exceeds deformation criterion')
            save_cycle(f'fig{number}_'+task['key'],r);rows.append(dict(s,key=task['key'],flexible=task['flexible']))
        if number==10:
            rigid=next(x for x in rows if x['key']=='rigid')
            for s in rows:
                for k in ['W_plus','peak_muscle_force','peak_positive_muscle_power']:s[k+'_over_rigid']=s[k]/rigid[k]
    elif number==6:
        # Panel (a) follows one upward-drive connected continuation per stiffness.
        for rw in sorted(set(t['r_w'] for t in tasks if t['panel']=='a')):
            state=Y0.copy();P=f6_groups(5.83,rw)
            for task in sorted([t for t in tasks if t['panel']=='a' and t['r_w']==rw],key=lambda x:x['mu']):
                s,r=cycle(P,task['mu'],y0=state);state=r['solution'].y[:,-1]
                if rw==7.5 and (abs(task['mu']-68.2717320127)<1e-5 or abs(task['mu']-86.3998743308)<1e-5):
                    s['stability']=floquet(P,task['mu'],r)
                    if not s['stability']['locally_stable']:raise NumericalFailure('rw=7.5 orbit unstable')
                rows.append(dict(s,key=task['key'],panel='a',mu=task['mu']))
        for task in [t for t in tasks if t['panel']=='b']:
            end=f6_boundary_task((task['t0_ms'],task['r_w'],task['root_mu']));roll='rollover' in end['endpoint_type']
            if not roll and not end['slope_crossing_found']:raise NumericalFailure('Unresolved operating-envelope endpoint')
            mu=end['envelope_mu']*(1. if roll else .9999);s,r=cycle(f6_groups(task['t0_ms'],task['r_w']),mu)
            if not s['quantitative']:raise NumericalFailure('Operating-envelope endpoint not admissible')
            if ('rollover' in task['endpoint_type'])!=roll:raise NumericalFailure('Endpoint mechanism changed')
            if roll:
                s['stability']=floquet(f6_groups(task['t0_ms'],task['r_w']),mu,r)
                if not s['stability']['locally_stable']:raise NumericalFailure('Rollover orbit unstable')
            rows.append(dict(s,mu=s['mu_hat'],key=task['key'],panel='b',t0_ms=task['t0_ms'],endpoint_type=end['endpoint_type']))
    elif number==7:
        for nw in f7_NWS:
            group=[t for t in tasks if t['N_W']==nw];P=f11_groups_with_nw(4.,nw);base,_=match(P,.22,group[0]['mu_guess']);mu=base['mu_hat']
            rigid,_=cycle(f11_groups_with_nw(6.,nw),mu,False)
            for task in group:
                s,r=cycle(f11_groups_with_nw(task['r_w'],nw),mu)
                if not s['quantitative']:raise NumericalFailure('Displayed Figure 7 case outside support')
                rows.append(dict(s,key=task['key'],N_W=nw,delta_phi_net_pp=100*(s['phi_net']-rigid['phi_net']),delta_phi_aero_pp=100*(s['phi_aero']-rigid['phi_aero'])))
    elif number in [8,9]:
        rows=[sensitivity(t,timing=number==9) for t in tasks]
        if number==8:write_json(OUTPUT/'fig8_boundaries.json',resolved_boundaries(8))
    elif number==11:
        rows=[map_point(t) for t in tasks];write_json(OUTPUT/'fig11_boundaries.json',resolved_boundaries(11))
    else:raise ValueError(number)
    write_json(OUTPUT/f'fig{number}.json',rows);write_csv(OUTPUT/f'fig{number}.csv',rows)
    return rows

def table4_compute():
    rows=[]
    def add(fig,key,q,s,**extra):rows.append(dict(figure=fig,key=key,q_max=q,max_slope=s,**extra))
    d=read_json(OUTPUT/'fig5.json')
    for r in d:add(5,r['key'],r['q_max'],r['max_slope'])
    d=[r for r in read_json(OUTPUT/'fig6.json') if r['panel']=='b']
    for roll in [False,True]:
        g=[r for r in d if ('rollover' in r['endpoint_type'])==roll]
        add(6,'rollover' if roll else 'slope',max(r['q_max'] for r in g),max(r['max_slope'] for r in g),count=len(g))
    d=read_json(OUTPUT/'fig7.json');add(7,'72 evaluated cases',max(r['q_max'] for r in d),max(r['max_slope'] for r in d))
    for b in read_json(OUTPUT/'fig8_boundaries.json'):
        d=b['point'];add(8,str(b['target']),d['support_slope']/SLOPE_FACTOR,d['support_slope'])
    d=[r for r in read_json(OUTPUT/'fig9.json') if r['status']=='computed' and r['quantitative']]
    for tag in ['rigid','rw2p5','rw2']:
        g=[r for r in d if r['tag']==tag]
        if not g:raise NumericalFailure('Missing timing family '+tag)
        r=max(g,key=lambda x:x['q_max']);add(9,tag,r['q_max'],r['max_slope'],support_q_max=max([r['q_max']]+[x['q_max'] for x in r['stencils']]),support_slope=r['support_slope'])
    d=[r for r in read_json(OUTPUT/'fig10.json') if r['key']!='rigid']
    add(10,'minimum',min(r['q_max'] for r in d),min(r['max_slope'] for r in d));add(10,'maximum',max(r['q_max'] for r in d),max(r['max_slope'] for r in d))
    for b in read_json(OUTPUT/'fig11_boundaries.json'):add(11,str(b['N_W']),b['point']['q_max'],b['point']['max_slope'])
    write_json(OUTPUT/'table4.json',rows);write_csv(OUTPUT/'table4.csv',rows);return rows

def compare_records(actual,expected):
    """Exact coordinate/status coverage plus mandatory finite metric checks."""
    records=[];lookup={r.get('key'):r for r in actual}
    if len(lookup)!=len(actual):raise NumericalFailure('Duplicate result coordinates')
    if set(lookup)!={r['key'] for r in expected}:raise NumericalFailure('Missing or unexpected result coordinates')
    coordinates=['r_w','N_W','t0','t0_ms','target','natural_frequency','frequency','F0']
    metrics=['amp','phi_net','phi_aero','Pi_flex','phase_q_minus_theta_deg','W_minus','W_struct','W_trunc','peak_muscle_power','peak_negative_muscle_power','peak_positive_power_W','mu','mu_hat','Omega','stroke_amp','q_max','max_slope','W_plus','W_m_net','W_aero','peak_muscle_force','peak_positive_muscle_power','peak_effective_force_N','mean_positive_power_W','positive_work_J','dimensional_mu_N_s','G_log','G_diff','support_slope','mu_H','eigen_mu_H','omega_H','energy_ratio','wingbeat_frequency','W_plus_over_rigid','peak_muscle_force_over_rigid','peak_positive_muscle_power_over_rigid','delta_phi_net_pp','delta_phi_aero_pp']
    for old in expected:
        new=lookup[old['key']]
        if old.get('status','computed')!=new.get('status','computed'):raise NumericalFailure('Changed availability status: '+old['key'])
        # Metadata applies to unavailable cases too; absence is not a False flag.
        if old.get('status')=='unavailable' or 'quantitative' in old:
            if type(new.get('quantitative')) is not bool:
                raise NumericalFailure('Missing/invalid quantitative classification: '+old['key'])
            if new['quantitative'] != old.get('quantitative',False):
                raise NumericalFailure('Changed admissibility: '+old['key'])
            if old.get('status')=='unavailable' and new['quantitative']:
                raise NumericalFailure('Unavailable case marked quantitative: '+old['key'])
        for label in ['panel','tag','endpoint_type','kind']:
            if label in old and new.get(label)!=old[label]:raise NumericalFailure('Wrong case label '+label+': '+old['key'])
        for k in coordinates:
            if k in old:
                required(new,[k])
                if abs(float(new[k])-float(old[k]))>1e-8:raise NumericalFailure('Wrong coordinate '+k+': '+old['key'])
        if old.get('status')=='unavailable':continue
        if new.get('kind')=='forced' and (new.get('role')!='diagnostic_prescribed_forcing_benchmark' or new.get('energy_definition')!='integral_abs_actuator_power'):raise NumericalFailure('Prescribed-forcing benchmark scope/energy convention missing')
        if new.get('stroke_amp',0)>0:
            if not accepted(new) or new.get('periodic') is False:raise NumericalFailure('Unconverged saved solution: '+old['key'])
            if new.get('quantitative') and new['max_slope']>.15:raise NumericalFailure('Unsupported solution marked quantitative')
        if 'G_log' in old:
            if len(new.get('stencils',[]))!=4 or len(new.get('estimates',[]))!=2 or not all(accepted(p) for p in new['stencils']):raise NumericalFailure('Missing/unconverged sensitivity stencil')
            support=max([new['max_slope']]+[p['max_slope'] for p in new['stencils']])
            if abs(support-new['support_slope'])>1e-12:raise NumericalFailure('Sensitivity support inconsistent with integrated stencils')
            if new.get('quantitative') and (support>.15 or new['fd_relative_change']>=.02):raise NumericalFailure('Unsupported derivative marked quantitative')
        fields=[k for k in metrics if k in old]
        if not fields:raise NumericalFailure('Empty benchmark: '+old['key'])
        required(new,fields)
        for k in fields:
            delta=abs(new[k]-old[k]);limit=2e-6+5e-4*abs(old[k])
            records.append(dict(key=old['key'],metric=k,reference=old[k],value=new[k],absolute_difference=delta,atol=2e-6,rtol=5e-4,status='PASS' if delta<=limit else 'FAIL'))
    return records

def peak_power_checks():
    rows=[]
    for n in [3,5,10]:
        for r in read_json(OUTPUT/f'fig{n}.json'):
            name='fig3' if n==3 else f'fig{n}_'+r['key'];path=OUTPUT/'cycles'/(name+'.npz')
            if not path.is_file():raise NumericalFailure('Missing cycle for peak-power test: '+name)
            with np.load(path) as d:power=r['mu_hat']*d['muscle_state']*d['theta_dot']
            val=float(max(power));required(r,['peak_positive_muscle_power'])
            rows.append(dict(key=name,positive_peak=val,absolute_peak=float(max(abs(power))),status='PASS' if abs(val-r['peak_positive_muscle_power'])<=1e-12 else 'FAIL'))
    return rows

def run_convergence(workers=1):
    resolve_production_boundary()
    records=[]
    for n in [5,10]:
        base=next(x for x in read_json(OUTPUT/f'fig{n}.json') if x['key']=='boundary');fine,_=cycle(dimensionless_groups(BEE,base['r_w']),base['mu_hat'],fine=True)
        for k in ['Omega','stroke_amp','max_slope','W_plus','W_m_net','peak_muscle_force','peak_positive_muscle_power']:
            change=abs(fine[k]-base[k])/max(abs(base[k]),1e-14);records.append(dict(key=f'{n}|{k}',relative_difference=change,status='PASS' if change<=5e-5 else 'FAIL'))
    write_json(OUTPUT/'convergence.json',records);return records

def run_sensitivity_convergence(workers=1):
    rows=[]
    seeds=read_json(REFERENCE/'search_seeds.json')
    for n in [8,9]:
        data=[r for r in read_json(OUTPUT/f'fig{n}.json') if r.get('status')=='computed' and r.get('quantitative')]
        group='target' if n==8 else 'tag'
        for value in sorted(set(r[group] for r in data)):
            base=min([r for r in data if r[group]==value],key=lambda r:r['r_w'] if n==8 else r['t0']);task=next(t for t in seeds[str(n)] if t['key']==base['key']);fine=sensitivity(task,timing=n==9,fine=True)
            change=abs(fine['G_log']/base['G_log']-1);rows.append(dict(key=f'{n}|{value}',relative_difference=change,status='PASS' if change<=5e-4 else 'FAIL'))
    write_json(OUTPUT/'sensitivity_convergence.json',rows);return rows

def verify_results():
    refs=read_json(REFERENCE/'production_values.json');checks=[];comparisons=[]
    try:
        if read_json(OUTPUT/'run_configuration.json')!=configuration_signature():raise NumericalFailure('Stale output configuration')
        checks.append(dict(key='configuration',status='PASS'))
    except (OSError,ValueError,RuntimeError) as e:checks.append(dict(key='configuration',status='FAIL',reason=str(e)))
    for n in range(2,12):
        try:
            rows=read_json(OUTPUT/f'fig{n}.json');result=compare_records(rows,refs['figures'][str(n)]);comparisons+=result
            checks.append(dict(key=f'Figure {n}',status='PASS' if result and all(x['status']=='PASS' for x in result) else 'FAIL'))
        except (OSError,ValueError,KeyError,TypeError,RuntimeError) as e:checks.append(dict(key=f'Figure {n}',status='FAIL',reason=str(e)))
    for n,group in [(8,'target'),(11,'N_W')]:
        try:
            bounds=read_json(OUTPUT/f'fig{n}_boundaries.json');expected=refs[f'boundaries{n}']
            if len(bounds)!=len(expected) or {x[group] for x in bounds}!={x[group] for x in expected}:raise NumericalFailure('Missing boundary family')
            for b in expected:
                x=next(z for z in bounds if z[group]==b[group]);required(x,['reported_rw'])
                if not accepted(x['point']):raise NumericalFailure('Boundary point lacks periodic/energy acceptance')
                if n==11:
                    required(x,['lower_rw','upper_rw','localized_reported_rw'])
                    if not x['lower_rw']<x['upper_rw']<=x['reported_rw'] or x['upper_rw']-x['lower_rw']>b['width'] or not accepted(x['localized_point']) or not x['localized_point'].get('quantitative') or x['localized_point']['max_slope']>.15:raise NumericalFailure('Independent localization/reporting mismatch')
                if n==11:
                    lower=x['lower_state']
                    if lower.get('quantitative'):raise NumericalFailure('Boundary lower state is still admissible')
                    if b['mechanism']=='connected_branch_target_search':
                        required(lower,['peak_amp'])
                        if lower.get('status')!='unavailable' or lower['peak_amp']>=.2 or not accepted(lower['peak_point']):raise NumericalFailure('Reachability boundary lacks a converged below-target peak')
                    elif lower.get('status')!='computed' or not accepted(lower) or lower['max_slope']<=.15:raise NumericalFailure('Deformation boundary lower state is not a converged diagnostic cycle')
                else:
                    p=x['point']
                    if len(p.get('stencils',[]))!=4 or not all(accepted(z) for z in p['stencils']) or p['support_slope']>.15 or p['fd_relative_change']>=.02:raise NumericalFailure('Boundary derivative support is incomplete')
                tolerance=4e-6 if n==8 else max(b['width'],3e-6)+2e-5
                if abs(x['reported_rw']-b['reported_rw'])>tolerance or not x['point'].get('quantitative'):raise NumericalFailure('Boundary localization/admissibility mismatch')
            checks.append(dict(key=f'Figure {n} boundaries',status='PASS'))
        except (OSError,ValueError,KeyError,RuntimeError) as e:checks.append(dict(key=f'Figure {n} boundaries',status='FAIL',reason=str(e)))
    for name in ['convergence','sensitivity_convergence']:
        try:
            data=read_json(OUTPUT/(name+'.json'))
            if len(data)!=(14 if name=='convergence' else 7) or len({r.get('key') for r in data})!=len(data) or any(r.get('status')!='PASS' for r in data):raise NumericalFailure('Incomplete/failed '+name)
            checks.append(dict(key=name,status='PASS'))
        except (OSError,ValueError,RuntimeError) as e:checks.append(dict(key=name,status='FAIL',reason=str(e)))
    try:
        powers=peak_power_checks();checks.append(dict(key='peak-positive power',status='PASS' if powers and all(x['status']=='PASS' for x in powers) else 'FAIL'))
    except (OSError,KeyError,ValueError,RuntimeError) as e:checks.append(dict(key='peak-positive power',status='FAIL',reason=str(e)))
    try:
        bound=read_json(OUTPUT/'production_boundary.json')
        required(bound,['root_rw','reported_rw','max_slope'])
        if abs(bound['reported_rw']-1.9814995134296145)>5e-6 or bound['max_slope']>.15:raise NumericalFailure('Production boundary mismatch')
        checks.append(dict(key='production boundary',status='PASS'))
    except (OSError,ValueError,RuntimeError) as e:checks.append(dict(key='production boundary',status='FAIL',reason=str(e)))
    for n in [2,3,4]:
        try:
            actual=read_json(OUTPUT/f'table{n}.json');expected=refs['tables'][str(n)]
            verify_table(n,actual,expected);checks.append(dict(key=f'Table {n}',status='PASS'))
        except (OSError,ValueError,KeyError,TypeError,RuntimeError) as e:checks.append(dict(key=f'Table {n}',status='FAIL',reason=str(e)))
    result=dict(status='PASS' if all(x['status']=='PASS' for x in checks) else 'FAIL',checks=checks,comparisons=comparisons)
    write_json(OUTPUT/'verification.json',result);write_csv(OUTPUT/'verification.csv',comparisons);return result

def verify_table(n,actual,expected):
    if n==2:
        required(actual,list(expected))
        if any(abs(actual[k]-expected[k])>1e-10 for k in expected):raise NumericalFailure('Table 2 mismatch')
        return
    if not actual or len(actual)!=len(expected):raise NumericalFailure('Incomplete table '+str(n))
    for a,b in zip(actual,expected):
        for k,x in b.items():
            if isinstance(x,(float,int)) and not isinstance(x,bool):
                required(a,[k])
                if abs(a[k]-x)>2e-6+5e-4*abs(x):raise NumericalFailure('Table mismatch '+str(n)+' '+k)
            elif isinstance(x,str) and a.get(k)!=x:raise NumericalFailure('Table coordinate mismatch '+k)

def verification_exit_status():
    try:d=read_json(OUTPUT/'verification.json')
    except (OSError,ValueError):return 2
    expected={f'Figure {n}' for n in range(2,12)}|{f'Table {n}' for n in [2,3,4]}|{'convergence','sensitivity_convergence','peak-positive power','production boundary','configuration','Figure 8 boundaries','Figure 11 boundaries'}
    checks=d.get('checks',[])
    return 0 if d.get('status')=='PASS' and len(checks)==len(expected) and {x.get('key') for x in checks}==expected and all(x.get('status')=='PASS' for x in checks) else 2

def enrich_report():
    d=read_json(OUTPUT/'verification.json');write_json(OUTPUT/'run_summary.json',dict(status=d['status'],scope='Corrected production suite; structural suite separately required',stage_b2='Execution evidence only; no automatic publication authorization'))

def floquet(P,mu,r):
    """Local variational flow; Phase 1/4 finite-difference RHS convention."""
    aero=Aerodynamics('direct',128);period=r['period'];initial=r['solution'].sol(r['t'][0]);fun=lambda t,y:n5.rhs(t,y,P,mu,True,aero)
    sol=solve_ivp(fun,(0,period),initial,method='DOP853',rtol=3e-11,atol=3e-13,max_step=.03,dense_output=True)
    def variational(t,z):
        y=sol.sol(t);eps=1e-6;J=np.column_stack([(fun(t,y+np.eye(6)[i]*eps)-fun(t,y-np.eye(6)[i]*eps))/(2*eps) for i in range(6)])
        return (J@z.reshape(6,6)).ravel()
    flow=solve_ivp(variational,(0,period),np.eye(6).ravel(),method='DOP853',rtol=1e-9,atol=1e-11,max_step=.03)
    if not sol.success or not flow.success:raise NumericalFailure('Floquet integration failed')
    eig=np.linalg.eigvals(flow.y[:,-1].reshape(6,6));phase=int(np.argmin(abs(eig-1)));trans=np.delete(eig,phase);error=float(abs(eig[phase]-1));largest=float(max(abs(trans)))
    return dict(neutral_error=error,largest_transverse_modulus=largest,locally_stable=bool(error<1e-4 and largest<1-1e-5))


def resolve_production_boundary():
    """Resolve theta_A=.20 deformation root independently of saved answers."""
    def point(rw):return match(dimensionless_groups(BEE,rw),.2,16.27)[0]
    root=brentq(lambda rw:point(rw)['max_slope']-.15,1.97,2.,xtol=3e-7,rtol=1e-9)
    safe=root+3e-6;s=point(safe)
    if not s['quantitative']:raise NumericalFailure('Production boundary interior margin failed')
    result=dict(root_rw=root,reported_rw=safe,max_slope=s['max_slope'],point=s)
    write_json(OUTPUT/'production_boundary.json',result);return result
