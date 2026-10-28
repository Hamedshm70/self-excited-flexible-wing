# Final public reproducibility candidate — verification report

Version 0.1.0-rc2; completed 2026-10-27. **READY AFTER MINOR CORRECTIONS:** scientific preparation and scoped numerical reproduction passed; software licensing, authorship/citation metadata and publication approval remain author decisions. No online publication occurred.

**Manuscript results: 13 tested, 13 PASS, 0 PARTIAL, 0 FAIL, 0 NOT TESTED**, within the documented supported/diagnostic scope. This does not mean every sampled parameter case converged: excluded numerical failures are explicitly listed below. Figure 1 is a schematic outside the computational generator.

## Reproduction results

| Manuscript result | Reproduction status | Numerical discrepancy | Notes |
|---|---|---|---|
| Figure 2 | PASS | max absolute difference 1.98111e-07; 1.498e-05 of allowed limit | 31 case IDs; 10 supported flags, 0 diagnostic flags, 0 unavailable. |
| Figure 3 | PASS | max absolute difference 2.3629e-07; 5.882e-05 of allowed limit | 1 case IDs; 1 supported flags, 0 diagnostic flags, 0 unavailable. |
| Figure 4 | PASS | max absolute difference 5.99675e-09; 9.963e-09 of allowed limit | 189 case IDs; 73 supported flags, 0 diagnostic flags, 0 unavailable. |
| Figure 5 | PASS | max absolute difference 4.76995e-08; 5.291e-05 of allowed limit | 4 case IDs; 4 supported flags, 0 diagnostic flags, 0 unavailable. |
| Figure 6 | PASS | max absolute difference 0.000352984; 0.0521 of allowed limit | 240 case IDs; 171 supported flags, 69 diagnostic flags, 0 unavailable. |
| Figure 7 | PASS | max absolute difference 1.10147e-06; 0.003458 of allowed limit | 72 case IDs; 72 supported flags, 0 diagnostic flags, 0 unavailable. |
| Figure 8 | PASS | max absolute difference 4.04636e-07; 0.003187 of allowed limit | 162 case IDs; 126 supported flags, 33 diagnostic flags, 3 unavailable. |
| Figure 9 | PASS | max absolute difference 0.0301186; 0.1162 of allowed limit | 44 case IDs; 33 supported flags, 2 diagnostic flags, 9 unavailable. |
| Figure 10 | PASS | max absolute difference 4.76995e-08; 0.0001192 of allowed limit | 7 case IDs; 7 supported flags, 0 diagnostic flags, 0 unavailable. |
| Figure 11 | PASS | max absolute difference 1.37489e-06; 0.001531 of allowed limit | 160 case IDs; 153 supported flags, 5 diagnostic flags, 2 unavailable. |
| Table 2 | PASS | derived constants within 1e-10 of targets | All reported rounded values checked. |
| Table 3 | PASS | within 2e-6 + 5e-4 times target magnitude | All ten onset-ratio entries match manuscript rounding. |
| Table 4 | PASS | within 2e-6 + 5e-4 times target magnitude | All 25 rows checked; rounded manuscript values checked separately. Figure 7 now uses all 72 cases. |

Counts use explicit quantitative flags; onset/rest records and the Figure 2 limiting-model benchmark are not assigned coupled-model support by inference. All 21 prescribed-forcing records retain their diagnostic benchmark role. All **910** expected figure-case identifiers are present. **17,848 scalar comparisons and 20 strict verification groups passed**. The 88 separate manuscript-value checks use half of the final reported rounding unit as absolute tolerance; approximate boundary values retain their displayed precision. Independent classification/coordinate inspection passed all 910 records, including 14 unavailable records.

The largest fraction of an allowed scalar tolerance is 0.116154, for `rw2p5|0.00290415656198` / `G_diff`: computed -518.564452599, target -518.594571172, absolute discrepancy 0.030118574. Numerical agreement is approximate, not an assertion of exact identity.

## Scientific implementation and Figure 7

The dSA dynamics use alpha2=r3+r4 and alpha3=r3*r4, kappa=0.17 and baseline t_o=5.83 ms. Timing-sensitivity cases, including 4.4 ms as a study value, are retained. Production uses direct complete-local-normal-velocity aerodynamic projection. Bilateral dimensional work, force and power use twice the half-system inertia. Positive work, net work, mean positive power, absolute peak and peak positive power remain distinct. Governing equations, baseline physical parameters, solvers and acceptance tolerances were not changed.

Figure 7 independently reproduced 18 cases per family at N_W=0.3,0.84,1.3,3.0: **57 negative and 15 positive net-work-fraction changes; all 15 positive cases in N_W=0.3; all 72 aerodynamic-work-fraction changes negative**. Gains are selected once per family at r_w=4 and then held fixed. Fractions were independently recomputed from their work definitions against freshly calculated rigid references. Maximum q=0.0756185439985, slope=0.104089342508, matching Table 4 values 0.0756 and 0.1041. The added family also recalibrated successfully from a gain guess displaced by 20%.

## Central Figure 10 point

| Quantity | Freshly calculated |
|---|---:|
| r_w | 1.98149951343 |
| mu_hat | 16.2727145998 |
| stroke_amp | 0.199999999996 |
| W_plus_over_rigid | 2.95135508843 |
| peak_force_over_rigid | 2.94649886524 |
| peak_positive_power_over_rigid | 3.01843803181 |
| peak_effective_force_mN | 67.9048496579 |
| mean_positive_power_mW | 0.154398787383 |
| positive_work_uJ | 1.50883369289 |
| frequency_Hz | 102.329891035 |

Positive work and peak positive power were independently reaggregated from the freshly integrated cycle arrays. The displayed dimensional values also match the final manuscript rounding.

## Numerical acceptance and robustness

Time-integrated production cycles pass state closure <=3e-7, amplitude drift <=5e-6, work drift <=3e-7 and relative energy imbalance <1e-7. Figure 4 shooting branches instead use shooting residual <1e-7, state closure <3e-7 and energy imbalance <1e-7; nine near-onset checks and 78 saved local Floquet classifications from Figures 4/6 passed. Quantitative support additionally requires max slope <=0.15. Diagnostic branches and unavailable cases are recorded separately from reference agreement. Sensitivity support covers base and perturbed solutions and its existing finite-difference criterion. The 14 cycle-refinement metric comparisons at the Figure 5/10 boundary points pass the 5e-5 relative criterion; seven selected sensitivity refinements pass 5e-4. Scalar agreement uses absolute tolerance 2e-6 plus relative tolerance 5e-4. Table 2 uses 1e-10. Boundary checks retain the distributed width/localization limits: Figure 8 reported coordinate 4e-6; Figure 11 max(reference width,3e-6)+2e-5. Structural fine/coarse checks retain 1e-4.

Structural verification passed **136 checks plus 11 supplemental phase/higher-mode checks**. Protocols remain distinct:

| Protocol | Status | Interpretation |
|---|---|---|
| One versus three modes | PASS | Independently matched stroke at equal stiffness. |
| Restored structural terms | PASS | Checked separately from modal truncation. |
| Isolated matched-stroke comparison | PASS | One-mode restored-term comparison at r_w=2. |
| Same drive | PASS | Fixed production drive at two stiffnesses. |
| Combined equal stiffness | PASS | Three modes and restored terms, independently matched stroke at two stiffnesses. |
| Independently refined boundaries | PASS | Different model boundaries localized separately. |

Combined-model root r_w=1.96974051485; admissible interior 1.96974551485. Higher-mode energy share 0.0334574421439%; mass-weighted squared-deformation share 0.00108829120072%. These are protocol-specific results, not universal model-error estimates.

## Verifier and graphical corrections

All 20 lightweight tests passed, including the new unavailable-classification regression and complete Figure 7 family test. The public `verification/check_supplemental.py` helper passed all 11 supplemental structural checks and all 15 disposable corrupted-output fixtures caused explicit aggregate verification failure: missing fields/results, nonfinite values, stale configuration, unconverged cycles, missing refinements/stencils, duplicate/missing cases, unsupported deformation, and unavailable-case classification/coordinate errors. Peak-positive-power discrepancies contribute to aggregate failure. Actual regenerated classifications remain correct; no scientific value was altered to satisfy the verifier.

All newly generated Figures 2–11 were visually inspected. Figure 4 restores the stable-rest shading, analytical-onset markers and the full onset range (the old upper axis limit clipped the curve). Figure 5 restores the red dotted/shaded boundary loop. Figure 7 includes all four families. Figure 9 restores open admissible boundary-adjacent symbols and the sensitivity hat. Labels, units, branches, boundaries, gaps and legends were checked. Colors, markers, dash patterns and layout retain some stylistic differences from the manuscript artwork; numerical coordinates and scientific meaning are preserved. The protected final PDFs for Figures 4/5/7/9 were rendered and inspected directly. Final manuscript artwork is unchanged.

## Environment, independence and runtime

Tested environment: macOS-26.5-arm64-64bit, Python 3.9.6 (default, May 22 2026, 11:13:45) , architecture arm64. Direct dependencies: NumPy 2.0.2, SciPy 1.13.1, Matplotlib 3.9.4. All installed distribution versions are recorded in verification_evidence/summary.json. Only this platform/interpreter combination was tested; the full dependency-compatible Python range and Windows/Linux were not verified.

A new virtual environment was created and dependencies installed only with requirements.txt, using the README procedure. The first network-restricted installation attempt failed; the allowed dependency-index retry succeeded. No old environment, production output or checkpoint was copied. Runtime imports were checked. macOS sandbox restrictions denied original research-tree contents, network access and writes outside the preparation candidate. Negative read/network probes passed. Original documentation, manuscript and hashes were read separately by the auditor; they were not numerical-runtime inputs.

Production used the documented individual `run_all.py --figure N` commands in parallel, then `--table 2 3 4`, unchanged cycle/sensitivity refinement functions, `--structural` and `--verify`. Every production calculation was cold initially. Only compatible accepted-cycle checkpoints created during this campaign were reused. Independent sensitivity families and the exact distributed boundary-loop bodies were prepared in parallel; their saved dictionaries exactly equal the documented workflow results. The serial default orchestration was not separately timed as one uninterrupted command.

Measured preparation/validation wall time to report completion: 72.5 minutes, including installation, orchestration and overlapping checks; not summed CPU time. See exact stage runtimes below. Serial/default runtime may be longer. A second source copy, initially empty of outputs/checkpoints, independently ran all 20 tests and Figures 7/10 in 369.5 seconds, with primary outputs and sources denied. All 79 cases passed target comparisons and had zero scalar difference from the primary run. This used the same fresh dependency environment and OS, so it is not a second-platform test.

| Stage | Measured seconds |
|---|---:|
| figure5 | 15.3 |
| figure3 | 29.0 |
| figure10 | 13.9 |
| figure2 | 80.4 |
| figure4 | 154.2 |
| figure7 | 344.9 |
| structural | 772.4 |
| figure9 | 1509.8 |
| figure6 | 2107.0 |
| figure8 | 2245.9 |
| figure11 | 3606.3 |
| tables | 2.1 |
| refinement | 4.2 |
| strict | 1.4 |
| central | 0.8 |
| faults_and_supplemental | 13.9 |
| Early refinement (overlapped; includes waiting for Figure 8) | 1720.0 |

## Explicit limitations and unavailable cases

Figure 2(a) retains its diagnostic prescribed-forcing/absolute-input-work qualification. All nine fresh 45-versus-90-cycle duration checks passed; 6 checked flexible cases exceed the coupled-model slope criterion, maximum 0.3266101. This is numerical verification of a limiting benchmark, not quantitative coupled-model or independent experimental validation. Historical comparison markers are internal model samples.

Fourteen sampled cases remain explicitly unavailable, matching the authoritative availability targets. Five are numerical failures (three unconverged Figure 8 cycles and two Figure 9 integration failures); they are not accepted solutions or proof of physical non-reachability. The other cases report local first-branch peaks below target, not global absence of other branches. Supported-domain curves are reproduced; failed/excluded cases are not silently discarded.

| Figure | Case | Recorded reason |
|---|---|---|
| 8 | `0.5\|3` | Unconverged cycle: {'periodic_state_rel': 0.11318697968681342, 'cycle_amp_rel': 0.06657000834358293, 'cycle_work_rel': 0.031749215133796554, 'energy_balance_rel': 4.4784470202588167e-10} |
| 8 | `1.0\|5.2` | Unconverged cycle: {'periodic_state_rel': 0.49981300699744224, 'cycle_amp_rel': 0.008916699118275808, 'cycle_work_rel': 0.0009751097449775053, 'energy_balance_rel': 1.964833895356109e-10} |
| 8 | `1.0\|5.25` | Unconverged cycle: {'periodic_state_rel': 0.5336157846165965, 'cycle_amp_rel': 0.000555700422036377, 'cycle_work_rel': 0.002167087252587497, 'energy_balance_rel': 2.633711612587975e-10} |
| 9 | `rw2p5\|0.0024` | Converged first-branch local peak 0.166965247911 below target at gain 14.1108813722 |
| 9 | `rw2p5\|0.0027064244137` | Converged first-branch local peak 0.193325865498 below target at gain 16.2743616423 |
| 9 | `rw2\|0.0024` | direct integration failed: Required step size is less than spacing between numbers. |
| 9 | `rw2\|0.0027064244137` | direct integration failed: Required step size is less than spacing between numbers. |
| 9 | `rw2\|0.00305197212794` | Converged first-branch local peak 0.133491664986 below target at gain 13.4511676252 |
| 9 | `rw2\|0.00344163828207` | Converged first-branch local peak 0.153269736446 below target at gain 15.8918989747 |
| 9 | `rw2\|0.00388105577904` | Converged first-branch local peak 0.173373139998 below target at gain 18.9861545335 |
| 9 | `rw2\|0.00437657671304` | Converged first-branch local peak 0.193671357663 below target at gain 22.8899555647 |
| 9 | `rw2\|0.0044` | Converged first-branch local peak 0.194575295699 below target at gain 23.0844471929 |
| 11 | `0.3\|3.147640625` | converged_continuous_local_rollover_below_target |
| 11 | `0.4\|2.63693451118` | converged_continuous_local_rollover_below_target |

Reduced structural representation, quasi-steady aerodynamics, local branch-search scope and a single tested platform limit interpretation and portability claims. No biological predictive accuracy is established by these numerical checks.

## Integrity, distribution and remaining author decisions

All 16 candidate scientific files retain their pre-test SHA-256 hashes. All 188 monitored original files retain their initial hashes, including the prior 31-file distribution, C3 reports, manuscript, reviewer response, final figure PDFs and selected authoritative validation records. No original production source, manuscript, reviewer response or historical C3 evidence was modified. An exhaustive snapshot of all filesystem/provider metadata was not taken; generated/relocated files naturally have new timestamps.

The clean release directory and ZIP contain the ten primary Python sources and supplemental verification helper, pinned requirements, five reference JSON files, README, current provenance/report, draft citation metadata, version, compact verification evidence and SHA256SUMS. Environments, checkpoints, calculated outputs, private paths, temporary files, internal logs and old packaging reports are excluded. Static portability/privacy scans and runtime dependency checks passed. Exact-source and reference changes are restricted to the necessary corrections described in PROVENANCE.md; the six other Python files and requirements remain identical to the prior candidate.

Technical suitability: the source package can be used for a GitHub repository or archival deposition, with documented numerical scope and portable relative paths. **Before public release, the author must approve a software license and software-authorship/citation metadata.** Both remain pending by explicit author instruction. No LICENSE was assigned; CITATION.cff is clearly marked draft. The title follows the current main.tex, as confirmed by the author. No journal acceptance, publication date, repository URL or DOI is invented. Cite a public repository/deposition in the manuscript only after it exists and its final metadata is verified. Nothing was pushed, uploaded, published or assigned a DOI.
