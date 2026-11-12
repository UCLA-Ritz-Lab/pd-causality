# Phase 0 — toolkit capabilities (CausalFM, TabPFN)

Date: 2026-11-11. Machine: 4 CPU cores, 30 GB RAM, **no GPU**. Env: `pd-qtl` (Python 3.13).

## What was installed

Added to `pd-qtl` via pip (CPU wheels):

| Package | Version | Note |
|---|---|---|
| torch | 2.7.1+cpu | |
| tabpfn | 2.0.9 | **Pinned.** CausalFM's bundled fork (`src/tabpfn`) imports `tabpfn.*` internals from this exact release |
| scikit-learn | 1.9.0 → **1.6.1** | Downgraded; tabpfn 2.0.9 requires <1.7. No PD-QTLs code uses sklearn |
| networkx, einops, tensorboard, huggingface-hub | latest | CausalFM data generators / training |

Rollback snapshot (pre-install): `conda list --explicit` and `pip freeze` saved in the session scratchpad
(`pd-qtl_before_phase0.{explicit,pip}.txt`). Reinstall with `pip install scikit-learn==1.9.0` if needed.

CausalFM-toolkit is a git clone at `tools/CausalFM-toolkit` (commit `bb74ef7`, 2026-01-14); git-ignored.
It is not pip-installable; import it by adding that folder to `sys.path` (the `causalfm` package imports `src.tabpfn...`).
TabPFN v2 regressor weights cached in `~/.cache/tabpfn/` (HuggingFace download failed; fallback URL worked).

## CausalFM: what's actually available

| Setting | Checkpoint released | Treatment A | Instrument Z | Status |
|---|---|---|---|---|
| Back-door CATE | `checkpoints_standard/best_model.pth` | binary | — | not tested (not our setting) |
| Front-door | `checkpoints_FD/best_model.pth` | binary | — | not tested |
| IV, binary Z | `checkpoints_IV_binary/best_model.pth` | **binary** | binary | **works** (see below) |
| IV, continuous Z | `checkpoints_IV_conti/best_model.pth` | **binary** | continuous | **broken: near-constant output** |

Model: per-feature transformer (TabPFN-v2-style), 2.8M params, in-context learning (no fitting).
Inputs: X (n × p covariates), Z (n × 1), A (n × 1), Y (n × 1), X_test. Output: CATE at X_test plus a
3-component GMM (`gmm_pi/mu/sigma`), which is the only uncertainty output. It is a predictive distribution over CATE, not a calibrated CI for the ATE.
Normalisation: standardise X and Y on training rows (`causalfm.data.normalize_data`). CATE comes out on the standardised-Y scale.

### Toolkit bugs found
- `causalfm.models.IVModel.estimate_cate` calls `self.model.estimate_cate`, which doesn't exist (the method is
  `estimate_cate_iv`), and passes args as (x, z, a, y). The real order is **`estimate_cate_iv(X, A, Y, Z, X_test, return_gmm=True)`**.
  Call `IVModel(...).model.estimate_cate_iv` directly. `IVModel.forward` has the same arg-order problem.
- The IV generators only produce **binary treatments** (`DATA_IV/gen_iv_data_*.py`: `A ~ Bernoulli(sigmoid(f(X,Z,U)))`).

### Smoke-test results (`scripts/py/phase0_causalfm_iv_smoke.py`, plus ad-hoc probes)

Bundled test sets (10 per setting, n = 1024, p = 10, 80/20 in-context split):

| Checkpoint | mean PEHE | mean \|ATE error\| | naive diff-in-means \|ATE error\| | pred-vs-true ATE corr / slope |
|---|---|---|---|---|
| IV binary-Z | 0.42 | **0.10** | 0.26 | 0.79 / 0.78 |
| IV continuous-Z | 0.58 | 0.24 | 0.18 | −0.30 / 0.00 |

- Continuous-Z checkpoint returns CATE ≈ −0.01 for every input (sd of predicted ATE across datasets 0.002
  vs true 0.33; same on synthetic data with true effect 0, 0.5, 1.0). It looks like a collapsed or mismatched checkpoint. **Don't use it.**
- Binary-Z checkpoint tracks the truth but shrinks toward 0 (synthetic, strong Z: true 0 / 0.5 / 1.0 → est 0.12 / 0.36 / 0.93).
- Shapes: ran without error for p = 3, 10, 30 and n up to 3000. It was trained on p = 10 (x0–x9), so other p are out of the training distribution.
- CPU runtime: ~6 s at n = 1024, ~14 s at n = 1000 with p = 30, ~43 s at n = 3000. Fine for applied work. For Phase 2 this is
  ~1–2 h per 1000 replicates per setting at n ≈ 1000.

## TabPFN v2 (official, 2.0.9)

- `TabPFNRegressor(device="cpu")`: Friedman1 test, 750 train / 250 test, 10 features: R² 0.96, 95% predictive-interval
  coverage 0.944, **65 s on CPU** (n_estimators = 4). Limits per docs: ≤10k rows, ≤500 features.
- Usable as a nonparametric first/second-stage learner (e.g. inside DML-IV or a control-function IV). Quantile output
  (`output_type="quantiles"`) gives predictive intervals, not CIs for a causal effect.
- CPU speed is the constraint: inside cross-fitting × simulation replicates it will be slow. Budget compute, or use a GPU.

## Implications for this project

1. **There's no released CausalFM model for our setting** (continuous exposure = methylation, instrument = 0/1/2 dosage).
   The only working IV checkpoint needs binary Z and binary A. Using it means dichotomising both the SNP (carrier vs not)
   and methylation (e.g. above/below median). That throws away information and changes the estimand to an effect of a
   dichotomised exposure, which violates exclusion if the rest of the methylation distribution matters.
2. **Decision (2026-11-11): (c) is the primary PFN arm and (a) is a labelled comparator** to quantify the cost of dichotomising. (b) isn't pursued. Options considered:
   - (a) Use CausalFM-IV-binary as-is with dichotomised Z/A, as a clearly-labelled comparator only.
   - (b) Retrain an IV model with a continuous-treatment prior (generators and `train_iv/*` scripts exist, but A would need
     changing to continuous). Training on CPU is likely impractical, so this needs a GPU and is in scope only if agreed.
   - (c) Use TabPFN as a component learner (DML-IV / control function) for the PFN arm instead of CausalFM.
3. Uncertainty from either PFN is a predictive distribution, not a frequentist CI. Phase 2 coverage must be checked empirically.
