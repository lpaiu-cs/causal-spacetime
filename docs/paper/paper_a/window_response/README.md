# Window-response study companion

Study date: 2026-09-21. Local branch: `codex/observation-window-response`.
This is the new focused study, not a replacement for the earlier frozen Paper A
submission. Start with [the Korean result report](research_report_ko.md) or
[the English manuscript source](article.tex). The final PDF is
`output/pdf/article.pdf` at the repository root.

## Result and population

The primary result is `PREDICTIONS_SUPPORTED`: all eight false-positive gates
and all eight prediction-accuracy gates passed. There are 51,200 independent
fixed-N posets, 3,200 per arm/condition, at N=300 and 1200. Amplitudes are 1/8
for R+/R-, sqrt(2) for R0, and 1/64 for D. The diamond amplitude was a
prespecified fallback, not a located 50%-power transition.

All inference is conditional on the specified plane-wave geometry, acquisition
window and known N. It is not unknown-window identification, general curvature
recovery, first Weyl detection, or evidence that the diamond's full order is
uninformative. Read the manuscript's sampling/numerical coverage distinction.

## Files and chronological provenance

All paths below are relative to `docs/paper/paper_a/` unless otherwise shown.
Producer SHAs refer to real commits reachable on this branch, not synthetic
review commits. Nothing in the frozen history was amended, rebased or squashed.

| Artifact | Producer / role |
|---|---|
| `supplementary_research_plan_2026-09-21.md`, `novelty_audit_2026-09-21.md` | Accepted scope recorded at `951420b`; historical pre-execution text retained |
| `window_response_theory.md` | First box derivation `75aecbb`; diamond reduction added at `e21e9de` after grid but before pilot |
| `window_response_g1_results.json` | `f3d4427056f327465fc45758aed281d262cb2850` |
| `window_response_protocol.md` | Pre-pilot freeze at `2cb6a1b` |
| `window_response_grid_results.json` | `34de74d`; 24 numerical-integral conditions, not holdout posets |
| `window_response_pilot_results.json` | `e21e9de`; 3,840 posets; raw statistic dictionaries included |
| `window_response_confirmation_freeze.json` | Computed at `a8997b2`, committed at `f6fe4ac` before confirmation |
| `window_response_confirmation_raw.npz` | Produced at `f6fe4ac`; shape `(8,2,3200,3)` |
| `window_response_confirmation_results.json` | `f6fe4ac`; raw/freeze hashes, all 16 intervals, formal verdict and cost |
| `window_response/analysis.json` | `498e78e`; independent recomputation, descriptive secondary comparisons and timing |
| `window_response/derive_exact.py` | Post-confirmation symbolic reconstruction of already frozen coefficients; reads no poset data |
| `window_response/manifest.json` | Raw and LF-normalized hashes of the delivered study files; not a pre-data registration |

NPZ axes are: curved design cell in freeze order; arm (flat,curved); replicate;
statistics `(f,chain3,escalation_count)`. Complete point arrays and relation
matrices are not retained. They are deterministically regenerable from the
recorded samplers and independent seed namespace
`[20260921,6401,region_index,N_index,arm,replicate]`.
The pilot uses namespace 6301 and quadrature 6201/6202. No seed was reused
between those populations. Record all software versions when regenerating.

The baseline statistical identity is standard U-statistic theory, not a new
result. The new calculations are the explicit window responses and the
diamond conditional-degree integrals. The normal approximation is empirical
at the confirmed conditions, not an exact finite-sample law.

## Software and inspection

Campaign environment: CPython 3.14.0, NumPy 2.5.3, SciPy 1.18.1 on Windows.
Plots: Matplotlib 3.11.2. Tests: pytest 9.1.1. Lint: Ruff 0.16.8.
Symbolic reconstruction was checked separately on CPython 3.11 with SymPy
1.14.0. The numerical campaign does not import SymPy.
The manuscript uses ordinary LaTeX packages and was built with Tectonic 0.17.0.

From the repository root, with the experiment dependencies installed:

```powershell
python -m pip install -r docs/paper/paper_a/window_response/requirements.txt
python -m pytest tests/test_window_response.py -q
python experiments/positive_control/window_response_report.py --audit-only
python docs/paper/paper_a/window_response/derive_exact.py
tectonic -X compile docs/paper/paper_a/window_response/article.tex --outdir output/pdf
```

The test and audit commands passed: 18 tests and an independent artifact audit.
`--audit-only` does not generate posets, overwrite analyses, or rerun timing.
The symbolic command requires SymPy and verifies the exact expressions by
polynomial/beta integrals, independently of the sampled statistics.
Full report regeneration without `--audit-only` refreshes plots and the
host-dependent timing benchmark; preserve the original `analysis.json` first.

This campaign reused `experiments/positive_control/p14_plane_wave.py` for the
scalar near-boundary fallback. New implementation is limited to
`window_response.py`, `window_response_g1.py`, `window_response_run.py`,
`window_response_confirm.py`, and `window_response_report.py` in that directory.
No earlier stage constants, outcomes or digest-locked inputs were changed.
During development, SciPy was temporarily added to the experiment extra in
`pyproject.toml`. Final packaging found that file was digest-locked by the older
oracle campaigns. It was restored byte-for-byte to the parent revision;
the new dependency is now declared in this study's requirements and CI install
step. Original new-study producer commits and their installed versions remain
preserved. No previous oracle campaign ran under the temporary change.

## Exact rerun from the local Git bundle

The companion bundle `outputs/window_response/window_response_2026-09-21.bundle`
contains this branch and its reachable history. It can be cloned without a
remote account. It is a local handoff artifact, not a public archive deposit.

```powershell
git -c core.autocrlf=true clone -b codex/observation-window-response window_response_2026-09-21.bundle window-response-rerun
Set-Location window-response-rerun
python -m venv .venv
.venv/Scripts/python.exe -m pip install numpy==2.5.3 scipy==1.18.1 matplotlib==3.11.2 pytest==9.1.1
.venv/Scripts/python.exe experiments/positive_control/window_response_report.py --audit-only
```

Use Python 3.14 to reproduce the recorded environment. On other platforms,
adjust the venv executable path. NumPy/SciPy version changes can change random
generation or rounding, so they produce a new reproduction record rather than
silently replacing the original.

For a full confirmation rerun, work **only in that disposable clone**. Preserve
the original results outside the clone, remove the two confirmation outputs
from the cloned checkout, and commit that deliberate change before running
`window_response_confirm.py confirm`. The runner requires a clean committed
state, an ancestral freeze, matching development hashes and absent output
files. It retains the original frozen forecasts and seed layout. Do not run it
on the original research checkout or change seeds after inspecting results.
Record the resulting new producer commit; compare statistics, not ZIP container
timestamps. The observed original run cost was 4,165 CPU seconds / 551 wall
seconds with eight workers; other hardware may differ.

JSON hashes in the original campaign were computed over Windows CRLF bytes.
Git stores their equivalent LF blobs. The audit accepts only reversible LF/CRLF
conversion while checking both recorded digests and pre-data content. The
`core.autocrlf=true` clone above restores the original input byte convention
for the unmodified historical runner. NPZ hashes are exact binary hashes.
The manifest records both text forms explicitly; content changes remain errors.

The original frozen stage commands, in order, were:

```text
python experiments/positive_control/window_response_g1.py
python experiments/positive_control/window_response_run.py grid
python experiments/positive_control/window_response_run.py pilot
python experiments/positive_control/window_response_confirm.py freeze
python experiments/positive_control/window_response_confirm.py confirm
```

They are a provenance recipe, not an instruction to overwrite existing outputs.
Regenerating development stages requires their own clean output-free commits,
recorded versions, original protocols and explicit comparison to the originals.

## Review observations

The audit initially caught Git LF versus Windows CRLF, then a NumPy boolean
serialization issue in the report generator. Both were corrected in separate
commits. Neither changed simulation data, forecasts, thresholds, primary scores,
or the experimental runner. All primary intervals were independently recomputed
from saved raw statistics, and a regression check rejects a changed JSON value.

All ten final PDF pages were rendered and visually inspected; the final build
has no missing references or overfull boxes. Eighteen targeted checks, exact
symbolic reconstruction, and the artifact audit passed. This is automated
verification, not independent human peer review.

The bibliography identifies primary sources and their inspection limits. The
earlier novelty audit is bounded; no exhaustive priority conclusion is implied.
AI assistance, absence of external funding and the declared absence of competing
interests are disclosed in the manuscript. A new journal submission and a public
archive deposit have not been performed.
