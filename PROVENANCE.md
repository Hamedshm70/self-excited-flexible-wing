# Provenance and scoped changes

Version 0.1.0-rc2 derives from the existing 31-file Stage C2 candidate, independently tested in Stage C3. Original packages and C3 reports remain preserved separately as historical evidence. The final manuscript is the authority for reported results; corrected validated production kernels determine the calculations.

The model, modal/full-model kernels, Figure 4 analysis, structural workflow, command interface and pinned requirements are unchanged. The corrected dSA sum/product rates and 5.83 ms / 0.17 baseline were already present. No equation, physical constant, numerical integration tolerance or acceptance tolerance was changed.

Necessary changes: extend Figure 7 from 54 to 72 cases with the validated N_W=0.3 family; extend Table 4 maxima to all four families; validate typed quantitative classifications and coordinates before skipping unavailable-case metrics; add targeted regression coverage; restore Figure 4 stable-rest shading, onset markers and full onset range; restore Figure 5 red dotted/shaded boundary loop; add Figure 9 open admissible endpoint symbols and the sensitivity hat; update documentation and manuscript title. Plot changes use computed data only and do not change numerical coordinates.

The 18 added comparison targets and calibration seeds were copied from the authoritative corrected Phase 4 Figure 7 records, without changing their numerical values. Existing 54 Figure 7 targets and all other figure targets are preserved. Only the Figure 7 Table 4 reference row changes to the validated four-family maximum. References are targets, not calculation outputs. Historical source arrays and reference/provenance.json retain their original assembly history and do not declare current validation status.

SHA256SUMS lists every public distribution file except itself, sorted by relative POSIX path, as SHA-256 followed by two spaces and filename. Locally generated outputs, environments and private preparation evidence are excluded from the ZIP. No private execution log or absolute research path is needed to reproduce the calculations.

Licensing options for author review include MIT/BSD-3-Clause (permissive), Apache-2.0 (permissive with explicit patent provisions), or GPL-3.0 (copyleft). This is a choice list, not a license assignment. Approved software authorship, repository/deposition information and DOI must come from the author/repository. No repository, upload or public identifier has been created.

The manuscript title follows the current main.tex, confirmed by the author during preparation. Licensing and software authorship remain pending by the author’s explicit instruction.

The supplemental verification helper makes the existing additional structural comparisons and disposable corruption tests reproducible from public outputs alone; it adds no model equations or simulation behavior.
