"""Portable corrected robustness protocols. Full integrations belong to Stage B2.
Kernels: Phase 2 modal_model/full_model; protocols: structural_checks/run_combined.
"""
from __future__ import annotations
import argparse
import numpy as np
from scipy.integrate import solve_ivp,simpson
from scipy.optimize import brentq
import modal_model as mm
import full_model as fm
import analysis as a
from parameters import BEE

METRICS=['Omega','stroke_amp','max_slope','W_plus','W_m_net','peak_muscle_force','peak_positive_muscle_power','Pi_flex']

def _evaluate_integrated(rw,n,mu,full=False,fine=False):
    model=(fm.FullerMultiModeModel if full else mm.MultiModeModel)(rw,n,nspan=128 if fine else 64,nx=801 if fine else 401)
    end=1200. if fine else 600.
    sol=solve_ivp(model.rhs,(0,end),model.initial_state(),args=(mu,),method='DOP853',rtol=3e-10 if fine else 3e-9,atol=3e-12 if fine else 3e-11,max_step=.0225 if fine else .045,dense_output=True,events=mm.upward_zero)
    if not sol.success:raise a.NumericalFailure(sol.message)
    ev=sol.t_events[0];ev=ev[ev>.75*end]
    if len(ev)<4:raise a.NumericalFailure('Incomplete modal cycles')
    r=mm.cycle_metrics(sol,ev[-2],ev[-1],model,mu,n_cycle=16001 if fine else 6001);prev=mm.cycle_metrics(sol,ev[-3],ev[-2],model,mu,n_cycle=2001)
    s=a.scalar_result(r);t=r['t'];v=r['theta_dot'];q=r['q'];qd=r['q_dot'];chi=model.chi[:,None];beta=model.beta[:,None];delta=chi*model.sigma
    dy=np.column_stack([model.rhs(ti,y,mu) for ti,y in zip(t,sol.sol(t).T)]);acc=dy[1];qdd=dy[3:2+2*n:2]
    E=.5*(v*v+r['theta']**2)+v*np.sum(beta*qd,axis=0)+.5*np.sum(chi*(qd*qd+model.r[:,None]**2*q*q),axis=0)+.5*v*v*np.sum(delta*q*q,axis=0)
    wt=0. if full else float(simpson(-v*acc*np.sum(delta*q*q,axis=0)-2*v*v*np.sum(delta*q*qd,axis=0),x=t))
    de=float(E[-1]-E[0]);res=r['W_m_net']-de-r['W_aero']-r['W_struct']-wt
    pf=(np.sum(beta*qdd,axis=0)+r['aero_stroke_load']-v*np.abs(v)/model.P['Nw'])*v
    s.update(mu_hat=mu,r_w=rw,n_modes=n,W_trunc=wt,energy_balance_rel=abs(res)/abs(r['W_m_net']),cycle_amp_rel=abs(r['stroke_amp']-prev['stroke_amp'])/r['stroke_amp'],cycle_work_rel=abs(r['W_m_net']-prev['W_m_net'])/abs(r['W_m_net']),peak_positive_muscle_power=float(max(mu*r['muscle_state']*v)),Pi_flex=float(max(abs(pf))/(r['W_aero']/r['period'])),phase_tip_minus_theta_deg=r['phase_tip_minus_theta_deg'],modal_energy_fraction=r['modal_energy_fraction'].tolist(),modal_rms_deformation_fraction=r['modal_rms_deformation_fraction'].tolist())
    if not a.accepted(s):raise a.NumericalFailure('Modal periodicity/energy acceptance failed')
    s['quantitative']=s['max_slope']<=.15
    return s

def evaluate(rw,n,mu,full=False,fine=False):
    """Checkpoint an independently accepted structural integration, by exact inputs."""
    import inspect,hashlib,json,time,os,scipy
    h=hashlib.sha256()
    for p in [a.PACKAGE/'modal_model.py',a.PACKAGE/'full_model.py',a.PACKAGE/'parameters.py',a.PACKAGE/'model.py']:
        h.update(p.read_bytes())
    for fn in [_evaluate_integrated,a.accepted,a.scalar_result]:h.update(inspect.getsource(fn).encode())
    identity=dict(kernel=h.hexdigest(),numpy=np.__version__,scipy=scipy.__version__,rw=float(rw),n=int(n),mu=float(mu),full=bool(full),fine=bool(fine))
    path=a.OUTPUT/'checkpoints/structural'/(hashlib.sha256(json.dumps(identity,sort_keys=True).encode()).hexdigest()+'.json')
    if path.exists():
        d=a.read_json(path)
        if d['identity']!=identity or not a.accepted(d['result']):raise a.NumericalFailure('Invalid structural checkpoint')
        return d['result']
    start=time.time();result=_evaluate_integrated(rw,n,mu,full,fine)
    tmp=path.with_suffix('.'+str(os.getpid())+'.tmp');a.write_json(tmp,dict(identity=identity,result=result,seconds=time.time()-start));tmp.replace(path)
    print('structural cycle',rw,n,mu,'full',full,'fine',fine,'seconds',round(time.time()-start,2),flush=True)
    return result

def matched(rw,n,guess,full=False):
    cache={}
    def ev(mu):
        if mu not in cache:cache[mu]=evaluate(rw,n,float(mu),full)
        return cache[mu]
    if abs(ev(guess)['stroke_amp']-.2)<2e-8:return ev(guess)
    lo,hi=.97*guess,1.03*guess
    for _ in range(8):
        if ev(lo)['stroke_amp']<.2<ev(hi)['stroke_amp']:break
        lo*=.95;hi*=1.05
    else:raise a.NumericalFailure('Structural matched target not bracketed')
    mu=brentq(lambda u:ev(u)['stroke_amp']-.2,lo,hi,xtol=3e-8,rtol=2e-10);s=ev(mu)
    if abs(s['stroke_amp']-.2)>2e-8:raise a.NumericalFailure('Structural target mismatch')
    return s

def difference(base,new):return {k:100*(new[k]/base[k]-1) for k in METRICS if base.get(k)}

def run():
    a.check_configuration();refs=a.read_json(a.REFERENCE/'production_values.json');bp=refs['figures']['10'][1];rw=bp['r_w'];mu=bp['mu_hat']
    one=matched(rw,1,mu);three=matched(rw,3,one['mu_hat']);modal=dict(protocol='matched_stroke_same_stiffness',one=one,three=three,differences_pct=difference(one,three),higher_mode_energy_pct=100*sum(three['modal_energy_fraction'][1:]),higher_mode_massweighted_squared_deformation_pct=100*sum(three['modal_rms_deformation_fraction'][1:]),refined_one=evaluate(rw,1,one['mu_hat'],fine=True),refined_three=evaluate(rw,3,three['mu_hat'],fine=True))
    same=[];combined=[]
    for value in [2.,rw]:
        prod,_=a.match(a.dimensionless_groups(BEE,value),.2,mu)
        full=evaluate(value,1,prod['mu_hat'],full=True)
        same.append(dict(protocol='same_drive_one_mode',r_w=value,truncated=prod,full=full,differences_pct=difference(prod,full),refined=evaluate(value,1,prod['mu_hat'],full=True,fine=True)))
        comp=matched(value,3,prod['mu_hat'],full=True);combined.append(dict(protocol='independently_matched_stroke_same_stiffness',r_w=value,production=prod,combined=comp,differences_pct=difference(prod,comp),refined=evaluate(value,3,comp['mu_hat'],full=True,fine=True)))
    isolated=matched(2.,1,same[0]['truncated']['mu_hat'],full=True)
    root=brentq(lambda x:matched(x,3,mu,full=True)['max_slope']-.15,1.95,rw,xtol=3e-7,rtol=1e-9);point=matched(root+5e-6,3,mu,full=True)
    if not point['quantitative']:raise a.NumericalFailure('Full-model interior boundary inadmissible')
    boundary=dict(protocol='independently_refined_boundaries',root_rw=root,reported_rw=root+5e-6,production=combined[1]['production'],combined=point,refined=evaluate(root+5e-6,3,point['mu_hat'],full=True,fine=True))
    isolated_refined=evaluate(2.,1,isolated['mu_hat'],full=True,fine=True)
    result=dict(modal=modal,same_drive=same,isolated_matched_full=isolated,isolated_matched_full_refined=isolated_refined,combined=combined,boundary=boundary)
    # Check the declared central metrics against saved corrected authority.
    checks=[]
    def compare(label,new,old,keys):
        a.required(new,keys)
        for key in keys:
            reference=old.get(key,old.get('peak_muscle_power') if key=='peak_positive_muscle_power' else None)
            if reference is None:raise a.NumericalFailure('Missing structural benchmark '+key)
            delta=abs(new[key]-reference);checks.append(dict(key=label+'|'+key,status='PASS' if delta<=2e-6+5e-4*abs(reference) else 'FAIL',value=new[key],reference=reference))
    target=refs['structural'];compare('modal one',one,target['modal']['one'],METRICS);compare('modal three',three,target['modal']['three'],METRICS)
    for i,x in enumerate(same):compare('same drive '+str(i),x['full'],target['same_drive'][i]['full'],METRICS)
    for i,x in enumerate(combined):compare('combined '+str(i),x['combined'],target['combined'][i]['combined'],METRICS)
    compare('independent boundary',point,target['boundary']['combined'],METRICS)
    compare('isolated matched full',isolated,target['matched_full']['full'],[k for k in METRICS if k in target['matched_full']['full']])
    pairs=[('modal one',one,modal['refined_one']),('modal three',three,modal['refined_three'])]
    pairs += [('same '+str(i),x['full'],x['refined']) for i,x in enumerate(same)]
    pairs += [('combined '+str(i),x['combined'],x['refined']) for i,x in enumerate(combined)]
    pairs += [('independent boundary',point,boundary['refined']),('isolated matched full',isolated,isolated_refined)]
    for label,base,fine in pairs:
        for key in METRICS:
            change=abs(fine[key]/base[key]-1)
            checks.append(dict(key=label+'|refinement|'+key,relative_change=change,status='PASS' if change<1e-4 else 'FAIL'))
        checks.append(dict(key=label+'|refined admissibility',status='PASS' if not base['quantitative'] or fine['quantitative'] else 'FAIL'))
    result.update(checks=checks,status='PASS' if checks and all(c['status']=='PASS' for c in checks) else 'FAIL');a.write_json(a.OUTPUT/'structural_verification.json',result);return result

if __name__=='__main__':
    argparse.ArgumentParser(description=__doc__).parse_args()
    raise SystemExit(0 if run()['status']=='PASS' else 2)
