"""Stage B1 regression checks: no production integration and no plot generation."""
from __future__ import annotations
import copy,json,math,unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
import numpy as np
import parameters as p
import model
import analysis as a
import modal_model as mm
import full_model as fm

class NumericTests(unittest.TestCase):
    def test_filter_independent_rates_and_peak(self):
        r3=math.log(1/p.KAPPA)/(p.T0*(1-p.KAPPA));r4=p.KAPPA*r3
        a2,a3=p.muscle_a2_a3(p.T0,p.KAPPA,p.BEE['wn'])
        self.assertAlmostEqual(a2,(r3+r4)/p.BEE['wn'],places=13);self.assertAlmostEqual(a3,r3*r4/p.BEE['wn']**2,places=13)
        self.assertAlmostEqual(r3*math.exp(-r3*p.T0),r4*math.exp(-r4*p.T0),places=12)
        self.assertAlmostEqual(a2,.7185322638000936,places=12);self.assertAlmostEqual(a3,.06411649090560806,places=12)
    def test_baseline_configuration(self):
        self.assertEqual(p.T0,.00583);self.assertEqual(p.KAPPA,.17);self.assertEqual(p.f6_T0S_MS,[3.,4.4,5.83,8.]);self.assertIn(7.5,p.f6_RW_A);self.assertEqual(p.f7_NWS,[.3,.84,1.3,3.])
    def test_coordinate_and_metric_coverage(self):
        ref=[dict(key='case',r_w=2.,stroke_amp=.2,quantitative=True)]
        with self.assertRaises(a.NumericalFailure):a.compare_records([],ref)
        with self.assertRaises(a.NumericalFailure):a.compare_records([dict(key='case',r_w=2.,quantitative=True)],ref)
        with self.assertRaises(a.NumericalFailure):a.compare_records([dict(key='case',stroke_amp=.2,quantitative=True)],ref)
        with self.assertRaises(a.NumericalFailure):a.compare_records(ref+ref,ref)
    def test_unavailable_classification_and_coordinates(self):
        ref=[dict(key='rw2p5|0.0024',status='unavailable',quantitative=False,r_w=2.5,t0=.0024,tag='rw2p5')]
        self.assertEqual(a.compare_records(copy.deepcopy(ref),ref),[])
        for changes in [dict(quantitative=True),dict(quantitative=0),dict(quantitative=None),dict(quantitative='False'),dict(r_w=2.0),dict(t0=float('nan')),dict(tag='rigid'),dict(status='computed')]:
            bad=copy.deepcopy(ref);bad[0].update(changes)
            with self.subTest(changes=changes),self.assertRaises(a.NumericalFailure):a.compare_records(bad,ref)
        for key in ['quantitative','t0','r_w','tag']:
            bad=copy.deepcopy(ref);del bad[0][key]
            with self.subTest(missing=key),self.assertRaises(a.NumericalFailure):a.compare_records(bad,ref)
    def test_figure7_complete_reference_families(self):
        rows=a.read_json(a.REFERENCE/'production_values.json')['figures']['7']
        self.assertEqual({nw:sum(x['N_W']==nw for x in rows) for nw in p.f7_NWS},{nw:18 for nw in p.f7_NWS})
        self.assertEqual(sum(x['delta_phi_net_pp']<0 for x in rows),57)
        positive=[x for x in rows if x['delta_phi_net_pp']>0]
        self.assertEqual(len(positive),15);self.assertEqual({x['N_W'] for x in positive},{.3})
        self.assertTrue(all(x['delta_phi_aero_pp']<0 for x in rows))
    def test_nonfinite_metric_rejected(self):
        with self.assertRaises(a.NumericalFailure):a.required({'a':float('nan')},['a'])
        with self.assertRaises(a.NumericalFailure):a.required({'a':True},['a'])
    def test_convergence_acceptance(self):
        r=dict(periodic_state_rel=1e-9,cycle_amp_rel=1e-8,cycle_work_rel=1e-9,energy_balance_rel=1e-10,max_slope=.149)
        self.assertTrue(a.accepted(r));r['cycle_work_rel']=1e-3;self.assertFalse(a.accepted(r))
    def test_settling_exhaustion_is_failure(self):
        fake=SimpleNamespace(nfev=1,success=True)
        metrics=dict(amp=.2,Omega=1.,state=np.ones(6),cycle_amp_rel=.01,periodic_state_rel=.01,cycle_frequency_rel=.01)
        with patch.object(a,'solve_ivp',return_value=fake),patch.object(a,'f8_cycle_metrics',return_value=metrics),patch.object(a,'f9_cycle_metrics',return_value=metrics):
            P=model.dimensionless_groups(p.BEE,2.)
            with self.assertRaises(RuntimeError):a.f8_settle(P,16.,None,warm=False)
            with self.assertRaises(RuntimeError):a.f9_settle(P,16.,True)
    def test_missing_outputs_fail_verification(self):
        with patch.object(a,'read_json',side_effect=lambda path: json.loads((a.REFERENCE/'production_values.json').read_text()) if Path(path).parent==a.REFERENCE else (_ for _ in ()).throw(FileNotFoundError('test missing output'))),patch.object(a,'write_json'),patch.object(a,'write_csv'):
            result=a.verify_results();self.assertEqual(result['status'],'FAIL');self.assertTrue(all(x['status']=='FAIL' for x in result['checks']))
    def test_incomplete_verification_report_fails_exit(self):
        for data in [{},{'status':'PASS','checks':[]},{'status':'PASS','checks':[dict(key='Figure 10',status='PASS')]}]:
            with patch.object(a,'read_json',return_value=data):self.assertEqual(a.verification_exit_status(),2)
    def test_peak_positive_is_distinct_from_absolute(self):
        P=model.dimensionless_groups(p.BEE,2.);r=dict(peak_muscle_force=1.,peak_muscle_power=3.,peak_positive_muscle_power=2.,W_plus=.5,period=2.)
        d=model.dimensional_metrics(P,r);self.assertAlmostEqual(d['peak_positive_power_W']/d['power_scale_W'],2.);self.assertAlmostEqual(d['mean_positive_power_W']/d['power_scale_W'],.25)
    def test_corrected_reference_targets(self):
        r=a.read_json(a.REFERENCE/'production_values.json');b=next(x for x in r['figures']['10'] if x['key']=='boundary')
        self.assertAlmostEqual(b['r_w'],1.9814995134296145,places=12)
        for key,val in [('W_plus_over_rigid',2.951355080586331),('peak_muscle_force_over_rigid',2.9464988520190225),('peak_positive_muscle_power_over_rigid',3.018438004629434)]:self.assertAlmostEqual(b[key],val,places=10)
        self.assertEqual(len(r['figures']['7']),72);self.assertEqual({x['N_W'] for x in r['figures']['7']},{.3,.84,1.3,3.})
        self.assertTrue(any(x['r_w']==7.5 for x in r['figures']['6']))
        for k in ['2','3','4','5','6','7','8','9','10','11']:self.assertEqual(len({x['key'] for x in r['figures'][k]}),len(r['figures'][k]))
    def test_table2_derived(self):
        ref=a.read_json(a.REFERENCE/'production_values.json')['tables']['2'];a.verify_table(2,a.table2_compute(),ref)
    def test_thoracic_loss_onset_regression(self):
        # Algebraic complex-stiffness root search only, not time integration.
        ref=a.read_json(a.REFERENCE/'production_values.json')['tables']['3']
        a.verify_table(3,a.table3_compute(),ref)
    def test_modal_full_rhs_consistency(self):
        y=np.array([.1,.2,.03,.04,.01,.02]);P=model.dimensionless_groups(p.BEE,2.);aero=model.Aerodynamics('direct',64)
        one=mm.MultiModeModel(2.,1);self.assertLess(float(max(abs(one.rhs(0,y,15)-model.rhs(0,y,P,15,True,aero)))),2e-7)
        full=fm.FullerMultiModeModel(2.,1);self.assertEqual(full.rhs(0,np.zeros(6),15).tolist(),np.zeros(6).tolist())
    def test_corrected_hopf_positive_gain(self):
        P=model.dimensionless_groups(p.BEE,2.)
        self.assertAlmostEqual(a.hopf_threshold(P),.2426050798,places=8)
        self.assertAlmostEqual(a.onset_eigen_threshold(P,True)[0],.2426050798,places=8)
    def test_bilateral_gain_conversion(self):
        P=model.dimensionless_groups(p.BEE,2.)
        self.assertAlmostEqual(2*P['I_theta'],3.9339474792280276e-11,places=20)
        gain=16.27271459983226*2*P['I_theta']*p.BEE['wn']/p.GAMMA
        self.assertAlmostEqual(gain,.009992433874690093,places=12)
    def test_table4_current_display_regressions(self):
        rows=a.read_json(a.REFERENCE/'production_values.json')['tables']['4']
        def one(fig,key):return next(x for x in rows if x['figure']==fig and x['key']==key)
        self.assertEqual(round(one(5,'rw3.0')['q_max'],4),.0366)
        self.assertEqual(round(one(5,'rw2.5')['max_slope'],4),.0772)
        self.assertEqual(one(6,'rollover')['count'],3)
        self.assertEqual(round(one(6,'rollover')['q_max'],4),.0852)
        self.assertEqual(round(one(7,'72 evaluated cases')['max_slope'],4),.1041)
        self.assertEqual(round(one(9,'rw2p5')['q_max'],4),.1053)
        self.assertEqual(round(one(9,'rw2p5')['support_q_max'],4),.1088)
        self.assertEqual(round(one(11,'0.3')['q_max'],4),.0980)
    def test_restored_one_mode_terms(self):
        P=model.dimensionless_groups(p.BEE,2.);y=np.array([.1,.2,.03,.04,.01,.02]);th,v,q,qd,f,fd=y
        ath,aq=model.Aerodynamics('direct',64).scalar(v,qd,P,True);delta=P['b_tq']/P['g_qt']*P['sig_cf'];mass=1+delta*q*q
        rt=15*f-th-ath-2*delta*q*qd*v;rq=-2*P['zeta']*P['r']*qd-(P['r']**2-P['sig_cf']*v*v)*q-aq;det=mass-P['b_tq']*P['g_qt']
        reference=np.array([v,(rt-P['b_tq']*rq)/det,qd,(mass*rq-P['g_qt']*rt)/det,fd,-P['a2']*fd-P['a3']*f-P['a3']*v])
        self.assertLess(float(max(abs(fm.FullerMultiModeModel(2.,1).rhs(0,y,15)-reference))),2e-7)
    def test_historical_data_never_current_targets(self):
        refs=a.read_json(a.REFERENCE/'production_values.json');self.assertEqual(refs['schema_version'],2);self.assertEqual(refs['baseline']['t0_s'],p.T0)
        self.assertNotIn('fig4_final',refs);self.assertNotIn('fig10',refs)

if __name__=='__main__':unittest.main(verbosity=2)
