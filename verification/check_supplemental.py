"""Verify supplemental structural fields and reject corrupted disposable outputs.
Run after the full production and structural workflows. Only outputs under this
package are read/written; actual production JSON and cycle arrays are preserved.
"""
import sys,json,copy,shutil,subprocess,os
from pathlib import Path
PACKAGE_ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(PACKAGE_ROOT))
run=PACKAGE_ROOT/'outputs/verification_extras';run.mkdir(parents=True,exist_ok=True)
import analysis as a
out=a.OUTPUT
s=a.read_json(out/'structural_verification.json');refs=a.read_json(a.REFERENCE/'production_values.json')['structural']
checks=[]
def check(key,value,reference):
 delta=abs(value-reference);checks.append(dict(key=key,value=value,reference=reference,absolute_difference=delta,status='PASS' if delta<=2e-6+5e-4*abs(reference) else 'FAIL'))
check('higher_mode_energy_pct',s['modal']['higher_mode_energy_pct'],100*sum(refs['modal']['three']['modal_energy_fraction'][1:]))
check('higher_mode_massweighted_squared_deformation_pct',s['modal']['higher_mode_massweighted_squared_deformation_pct'],100*sum(refs['modal']['three']['modal_rms_deformation_fraction'][1:]))
for name,x,y in [('modal one',s['modal']['one'],refs['modal']['one']),('modal three',s['modal']['three'],refs['modal']['three']),('boundary',s['boundary']['combined'],refs['boundary']['combined']),('isolated matched',s['isolated_matched_full'],refs['matched_full']['full'])]+[('same drive '+str(i),x['full'],refs['same_drive'][i]['full']) for i,x in enumerate(s['same_drive'])]+[('combined '+str(i),x['combined'],refs['combined'][i]['combined']) for i,x in enumerate(s['combined'])]:
 check(name+' phase',x['phase_tip_minus_theta_deg'],y.get('phase_tip_minus_theta_deg',y.get('phase_q_minus_theta_deg')))
b=s['boundary'];rb=refs['boundary']
check('independent boundary peak positive power difference (%)',100*(b['combined']['peak_positive_muscle_power']/b['production']['peak_positive_muscle_power']-1),100*(rb['combined'].get('peak_positive_muscle_power',rb['combined'].get('peak_muscle_power'))/a.read_json(a.REFERENCE/'production_values.json')['figures']['10'][1]['peak_positive_muscle_power']-1))
a.write_json(run/'supplemental_checks.json',checks)
assert len(checks)==11 and all(x['status']=='PASS' for x in checks)
if '--structural-only' in sys.argv:
 print('11 supplemental structural checks PASS');raise SystemExit(0)
# Disposable output fixture only; production files and source files are never mutated.
fixture=run/'fault_fixture';fixture.mkdir(exist_ok=True)
for p in out.glob('*.json'):shutil.copyfile(p,fixture/p.name)
(fixture/'cycles').mkdir(exist_ok=True)
for p in (out/'cycles').glob('*.npz'):os.link(p,fixture/'cycles'/p.name) if not (fixture/'cycles'/p.name).exists() else None
original={p.name:p.read_text() for p in fixture.glob('*.json')}
a.OUTPUT=fixture
baseline=a.verify_results();assert baseline['status']=='PASS' and a.verification_exit_status()==0
faults=[]
def fault(label,name,mutation):
 for k,v in original.items():(fixture/k).write_text(v)
 p=fixture/name
 if mutation is None:p.unlink()
 else:
  d=json.loads(p.read_text());mutation(d);p.write_text(json.dumps(d))
 result=a.verify_results();code=a.verification_exit_status();faults.append(dict(key=label,status='PASS' if result['status']=='FAIL' and code!=0 else 'FAIL',verification_status=result['status'],exit_status=code,failed_checks=[x for x in result['checks'] if x['status']!='PASS']))
fault('required field removed','fig10.json',lambda d:d[1].pop('W_plus'))
fault('required numerical result removed','fig10.json',None)
fault('nonfinite metric','fig10.json',lambda d:d[1].update(W_plus=float('nan')))
fault('configuration fingerprint corrupted','run_configuration.json',lambda d:d.update(numpy='invalid'))
fault('unconverged solution','fig10.json',lambda d:d[1].update(periodic_state_rel=1.0))
fault('required convergence record removed','convergence.json',lambda d:d.pop())
fault('failed positive peak comparison','fig10.json',lambda d:d[1].update(peak_positive_muscle_power=d[1]['peak_positive_muscle_power']*1.01))
fault('duplicate case identifier','fig10.json',lambda d:d.append(copy.deepcopy(d[1])))
fault('missing case identifier','fig10.json',lambda d:d.pop())
fault('unsupported deformation marked quantitative','fig10.json',lambda d:d[1].update(max_slope=.16,quantitative=True))
fault('required sensitivity stencil removed','fig8.json',lambda d:next(x for x in d if 'stencils' in x)['stencils'].pop())
fault('unavailable marked quantitative','fig9.json',lambda d:next(x for x in d if x.get('status')=='unavailable').update(quantitative=True))
fault('unavailable classification missing','fig9.json',lambda d:next(x for x in d if x.get('status')=='unavailable').pop('quantitative'))
fault('unavailable classification numeric','fig9.json',lambda d:next(x for x in d if x.get('status')=='unavailable').update(quantitative=0))
fault('unavailable coordinate changed','fig9.json',lambda d:next(x for x in d if x.get('status')=='unavailable').update(t0=99.))
a.write_json(run/'fault_checks.json',faults)
a.OUTPUT=out
# Explicit positive power failure must participate in aggregate failure.
p=next(x for x in faults if x['key']=='failed positive peak comparison')
assert any(x['key']=='peak-positive power' for x in p['failed_checks'])
print('Supplemental structural checks:',len(checks),'Fault rejection checks:',len(faults),'rejected:',sum(x['status']=='PASS' for x in faults),flush=True)
assert len(faults)==15 and all(x['status']=='PASS' for x in faults)
