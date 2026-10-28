# A Self-Excited Wing–Thorax Model of Asynchronous Insect Flight: Emergent Wingbeat and the Cost of Wing Flexibility

Python reproducibility package for the manuscript by Seyed Hamed Mirtalebi and Rozhin Hajian (unpublished manuscript, Bioinspiration & Biomimetics resubmission). Publication status and DOI are not asserted. This candidate computes Figures 2–11 and Tables 2–4 using a self-excited wing–thorax model with a reduced structural mode and quasi-steady, direct local-velocity aerodynamic drag. It solves emergent periodic wingbeats, computes work, actuator demand and sensitivities, and locates deformation-supported operating boundaries. Figure 1 is a schematic and is not generated.

## Installation

Python 3.9–3.12 is the dependency-compatible range; the verification report identifies the actually tested version/platform. From this directory:

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
source .venv/bin/activate
```

Windows: use `.venv\Scripts\python.exe` and `.venv\Scripts\Activate.ps1`. Dependencies are pinned in requirements.txt. Installation needs the package index; computation needs no network, MATLAB, LaTeX, private input data, or original research directory.

## Reproduce and verify

```sh
python -B run_all.py --tests
python -B run_all.py
python -B run_all.py --structural
python -B run_all.py --verify
python -B verification/check_supplemental.py
```

Run from the package directory. The default workflow independently integrates Figures 2–11, generates Tables 2–4 and figure PDFs/PNGs, performs cycle/sensitivity refinement and strict verification. The separate structural workflow tests one/three-mode convergence, restored terms, isolated matched stroke, same drive, combined equal stiffness and independently refined boundaries. Both production and structural verification must pass. The supplemental command checks phase/higher-mode fields and 15 disposable corrupted-output fixtures, preserving actual production results. Lightweight tests alone do not establish reproduction.

Results and verification reports are written to `outputs/`; figures to `outputs/figures/`. To calculate selected results or render already computed data:

```sh
python -B run_all.py --figure 7 10
python -B run_all.py --figure 4 5 7 9 --plot-only
python -B run_all.py --table 2 3 4
```

Table 4 requires the associated figure/boundary outputs. Plot-only and verify modes do not simulate missing cases. Partial runs intentionally fail complete-suite verification. Allow several hours and several GB of working storage; measured campaign runtimes are reported in REPRODUCIBILITY_REPORT.md. `--workers` does not parallelize every continuation. Independent figure commands can run in separate processes; preserve each command's serial continuation order.

Exact-input, numerical-kernel and dependency fingerprints gate locally generated accepted-cycle checkpoints. To start cold, move/remove your own `outputs/` directory first. Do not combine outputs made with different source/reference versions. Load only your own or trusted checkpoints: Python pickle checkpoints are not safe against malicious files. No checkpoints or precomputed simulation outputs are distributed as calculated results.

## Scientific conventions and limitations

Finite-amplitude production cases use a purely elastic thorax; Table 3 separately examines harmonic thoracic structural loss. The dSA filter uses alpha2=r3+r4 and alpha3=r3*r4, selected kappa=0.17 and t_o=5.83 ms. The 3, 4.4, 5.83 and 8 ms timing studies are retained; 4.4 ms is not the baseline. Bilateral dimensional work/force/power use twice the half-system inertia. Mean positive power, peak positive power and absolute peak power are distinct.

Figure 7 includes all four N_W families (0.3, 0.84, 1.3, 3.0), 18 stiffnesses each. Each family's drive is calibrated at r_w=4 to stroke half-amplitude 0.22 and then held fixed, including its rigid reference. The target JSONs are comparisons and search seeds, never substitutes for independently integrated solutions. Figure 7/Table 4 use all 72 cases.

Quantitative coupled-model results require maximum wing slope <=0.15 plus periodic and energy acceptance. Diagnostic continuations are separated; unavailable cases retain explicit reasons and are never treated as quantitative. Some numerical failures are part of the documented study and do not prove physical non-reachability. Local branch searches do not prove global uniqueness. Reduced structural and quasi-steady aerodynamic assumptions restrict physical interpretation; numerical consistency is not independent experimental validation.

Figure 2(a) is a diagnostic prescribed-forcing limiting-model benchmark using absolute input work, separate from the positive-work definition of the self-excited model. Six of nine checked benchmark flexible cases exceed the coupled-model deformation criterion. Figure 2 historical markers are internal historical model samples, not independent experiments or literature validation. Structural protocol differences are scoped checks, not universal model-error bounds.

Production acceptance: state closure <=3e-7, amplitude drift <=5e-6, work drift <=3e-7, relative energy imbalance <1e-7. Scalar comparison tolerance is |computed-reference| <=2e-6+5e-4|reference|. Figure 4 shooting branches use residual <1e-7 plus state/energy checks and local stability checks. Selected cycle refinement <=5e-5 and sensitivity refinement <=5e-4. Additional boundary/structural tolerances and actual outcomes are documented in REPRODUCIBILITY_REPORT.md.

## Integrity, provenance and citation

Verify the distributed files with `shasum -a 256 -c SHA256SUMS` (GNU: `sha256sum -c SHA256SUMS`). Hashes detect alteration, not publisher authenticity. See PROVENANCE.md for source history and changes, and REPRODUCIBILITY_REPORT.md for verification scope and limitations.

Associated manuscript: Seyed Hamed Mirtalebi and Rozhin Hajian, “A Self-Excited Wing–Thorax Model of Asynchronous Insect Flight: Emergent Wingbeat and the Cost of Wing Flexibility,” unpublished manuscript. Manuscript DOI: **not assigned/verified**. Software citation: use this package title and version, approved software authors, repository URL and archival DOI only after the author supplies/approves them. CITATION.cff remains explicit draft metadata while software authorship is unresolved. No DOI is invented.

No software license is assigned. Public reuse permission, software authorship/citation and publication authorization remain author decisions. Nothing has been uploaded or published.
