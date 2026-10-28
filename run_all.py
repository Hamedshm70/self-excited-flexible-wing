"""Corrected package interface. Full integrations are computationally expensive."""
from __future__ import annotations
import argparse,os,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.dont_write_bytecode=True
os.environ['MPLCONFIGDIR']=str(HERE/'outputs/.mplconfig')
os.environ.setdefault('OPENBLAS_NUM_THREADS','1');os.environ.setdefault('OMP_NUM_THREADS','1')

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--figure',nargs='+',type=int,choices=range(2,12));p.add_argument('--table',nargs='+',type=int,choices=[2,3,4]);p.add_argument('--workers',type=int,default=1)
    p.add_argument('--verify',action='store_true',help='Strictly verify all existing production results; writes verification evidence only')
    p.add_argument('--plot-only',action='store_true',help='Render only previously computed output, never reference targets')
    p.add_argument('--tests',action='store_true',help='Lightweight tests; no integrations or plots')
    p.add_argument('--structural',action='store_true',help='Expensive separate structural robustness suite')
    args=p.parse_args()
    if args.workers<1:p.error('workers must be positive')
    if args.tests:
        import unittest,tests_numeric
        result=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromModule(tests_numeric));return 0 if result.wasSuccessful() else 2
    import analysis as a
    if args.verify:a.verify_results();return a.verification_exit_status()
    if args.structural:
        import structural_verification as s
        return 0 if s.run()['status']=='PASS' else 2
    if args.plot_only and not args.figure:p.error('--plot-only requires --figure')
    figs=args.figure if args.figure is not None else ([] if args.table else list(range(2,12)))
    for n in figs:
        if not args.plot_only:a.compute_figure(n,args.workers)
        import figures
        figures.render(n)
    if args.plot_only:return 0
    tables=args.table if args.table is not None else ([2,3,4] if args.figure is None else [])
    for n in tables:
        a.check_configuration()
        data=a.table2_compute() if n==2 else a.table3_compute() if n==3 else a.table4_compute()
        a.write_json(a.OUTPUT/f'table{n}.json',data)
    if args.figure is None and args.table is None:
        a.run_convergence(args.workers);a.run_sensitivity_convergence(args.workers);a.verify_results();a.enrich_report();return a.verification_exit_status()
    return 0
if __name__=='__main__':raise SystemExit(main())
