"""Corrected Phase 4 onset/shooting workflow, with no private infrastructure."""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root
import model
import analysis as a
from parameters import BEE
from model import dimensionless_groups,Aerodynamics,rhs
from analysis import onset_jacobian,hopf_threshold
RWS=(1.3,2.,3.5)
def jac(P,mu):
    return onset_jacobian(P,mu,True)

def growth(P,mu):return float(np.linalg.eigvals(jac(P,mu)).real.max())

def near_onset_compute():
 rows=[]
 for r in RWS:
  P=dimensionless_groups(BEE,r);h=hopf_threshold(P);ev,V=np.linalg.eig(jac(P,h));i=np.argmax(ev.real);i=i if ev[i].imag>0 else np.argmin(abs(ev-np.conj(ev[i])))
  omega=ev[i].imag;v=V[:,i]/V[0,i];aero=Aerodynamics('direct')
  for target in (.001,.0005,.0001):
   velocity=omega*target;seed=np.real(-1j*v)*target
   def integrate(x,dense=False):
    mu=x[4];y=np.r_[0.,velocity,x[:4]]
    return solve_ivp(lambda t,y:rhs(t,y,P,mu,True,aero),(0,x[5]),y,method='DOP853',rtol=3e-11,atol=3e-13,max_step=.03,dense_output=dense)
   def residual(x):
    s=integrate(x);return (s.y[:,-1]-np.r_[0.,velocity,x[:4]])/velocity
   x0=np.r_[seed[2:],h+target,2*np.pi/omega]
   sol=root(residual,x0,tol=1e-9);x=sol.x;s=integrate(x,True);Y=s.sol(np.linspace(0,x[5],2401));amplitude=np.ptp(Y[0])/2
   # Finite-difference flow map; nontrivial Floquet multipliers exclude phase.
   y=np.r_[0.,velocity,x[:4]];columns=[]
   for j in range(6):
    dy=np.eye(6)[j]*1e-7
    ends=[]
    for sign in (1,-1):
     ss=solve_ivp(lambda t,y:rhs(t,y,P,x[4],True,aero),(0,x[5]),y+sign*dy,method='DOP853',rtol=3e-11,atol=3e-13,max_step=.03)
     ends.append(ss.y[:,-1])
    columns.append((ends[0]-ends[1])/2e-7)
   mult=np.linalg.eigvals(np.column_stack(columns));phaseidx=np.argmin(abs(mult-1));non=np.delete(mult,phaseidx)
   row=dict(r_w=r,target_linear_amp=target,stroke_amp=amplitude,mu_hat=x[4],mu_H=h,mu_minus_mu_H=x[4]-h,period=x[5],residual_scaled=float(np.max(abs(residual(x)))),root_success=sol.success,phase_multiplier_abs=float(abs(mult[phaseidx])),max_nontrivial_multiplier=float(max(abs(non))))
   rows.append(row);print(row,flush=True)
 return rows
def shooting(rw,mu,small,seedamp=None,fine=False):
 P=dimensionless_groups(BEE,rw);h=hopf_threshold(P)
 z=next(x for x in small if x['r_w']==rw and x['target_linear_amp']==.0001)
 eig,V=np.linalg.eig(jac(P,h));ii=np.argmin(abs(eig-1j*(2*np.pi/z['period'])));w=eig[ii].imag;v=V[:,ii]/V[0,ii]
 amp=seedamp or (mu-h)*z['stroke_amp']/z['mu_minus_mu_H'];seed=np.real(-1j*v)*amp
 aero=model.Aerodynamics('direct',128 if fine else 64)
 def flow(x):
  y=np.r_[0.,np.exp(x[0]),x[1:5]]
  return solve_ivp(model.rhs,(0.,x[5]),y,args=(P,mu,True,aero),method='DOP853',rtol=3e-11 if not fine else 3e-12,atol=3e-13 if not fine else 3e-14,max_step=.03 if not fine else .015,dense_output=True)
 def residual(x):return (flow(x).y[:,-1]-np.r_[0.,np.exp(x[0]),x[1:5]])/np.exp(x[0])
 x0=np.r_[np.log(w*amp),seed[2:],2*np.pi/w]
 rt=root(residual,x0,tol=1e-9);err=float(max(abs(residual(rt.x))));
 if err>=1e-7 or rt.x[-1]<=0:raise a.NumericalFailure(str((rw,mu,rt.message,err)))
 sol=flow(rt.x);r=model.cycle_metrics(sol,0.,rt.x[-1],P,mu,True,aero,n_cycle=6001 if not fine else 16001);r['solution']=sol
 r['shooting_residual']=err
 s=a.enrich(r,P,mu,True);s.update(mu_H=h,shooting_residual=err,root_success=bool(rt.success),r_w=rw,source='periodic shooting',P=P)
 return s,r

def compute(workers=1):
    tasks=a.read_json(a.REFERENCE/'search_seeds.json')['4'];small=near_onset_compute()
    if any(x['residual_scaled']>=1e-7 or x['max_nontrivial_multiplier']>=1 for x in small):raise a.NumericalFailure('Small-orbit verification failed')
    rows=[]
    for task in tasks:
        rw=task['r_w'];P=dimensionless_groups(BEE,rw);h=hopf_threshold(P)
        if task['kind']=='onset':
            e=a.onset_eigen_threshold(P,True)
            if abs(h-e[0])>1e-8:raise a.NumericalFailure('Onset/Jacobian mismatch')
            rows.append(dict(key=task['key'],kind='onset',r_w=rw,mu_H=h,eigen_mu_H=e[0],omega_H=e[1]))
        elif task['mu_hat']<=h:rows.append(dict(key=task['key'],kind='branch',r_w=rw,mu_hat=task['mu_hat'],mu_H=h,stroke_amp=0.,source='linearly stable rest'))
        else:
            s,r=shooting(rw,task['mu_hat'],small)
            st=a.floquet(P,task['mu_hat'],r)
            if not st['locally_stable']:raise a.NumericalFailure('Unstable plotted orbit')
            rows.append(dict(s,key=task['key'],kind='branch',stability=st))
    a.write_json(a.OUTPUT/'fig4_small_orbits.json',small)
    return rows
