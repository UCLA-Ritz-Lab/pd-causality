"""Phase 0 smoke test: run released CausalFM IV checkpoints on the toolkit's bundled IV test sets.

Run from repo root:  python scripts/py/phase0_causalfm_iv_smoke.py
Reports PEHE and ATE error per dataset vs. a naive difference-in-means baseline, plus CPU runtime.
"""
import glob
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch

TOOLKIT = Path(__file__).resolve().parents[2] / "tools" / "CausalFM-toolkit"
sys.path.insert(0, str(TOOLKIT))

from causalfm.data import normalize_data, normalize_ite  # noqa: E402
from causalfm.evaluation import compute_ate_error, compute_pehe  # noqa: E402
from causalfm.models import IVModel  # noqa: E402

torch.manual_seed(0)
rows = []
for setting in ["conti", "binary"]:
    model = IVModel.from_pretrained(str(TOOLKIT / f"checkpoints/checkpoints_IV_{setting}/best_model.pth"), device="cpu")
    model.eval_mode()
    files = sorted(glob.glob(str(TOOLKIT / f"DATA_IV/iv_{setting}_TEST/*.csv")))
    for f in files:
        df = pd.read_csv(f)
        x_cols = [c for c in df.columns if c.startswith("x")]
        X_norm, Y_norm, _, y_scaler = normalize_data(df[x_cols].values, df["outcome"].values, df["y0"].values, df["y1"].values)
        ite, _ = normalize_ite(df["y0"].values, df["y1"].values, y_scaler)
        X = torch.FloatTensor(X_norm)
        Z = torch.FloatTensor(df["z"].values).unsqueeze(1)
        A = torch.FloatTensor(df["treatment"].values).unsqueeze(1)
        Y = torch.FloatTensor(Y_norm).unsqueeze(1)
        n = int(0.8 * len(X))
        t0 = time.perf_counter()
        # IVModel.estimate_cate in the toolkit calls a non-existent method and passes args as (x, z, a, y);
        # the underlying model's method is estimate_cate_iv(X, A, Y, Z, X_test). Call it directly.
        with torch.no_grad():
            res = model.model.estimate_cate_iv(X[:n], A[:n], Y[:n], Z[:n], X[n:], return_gmm=True)
        dt = time.perf_counter() - t0
        cate = res["cate"].cpu().numpy().ravel()
        a, y = df["treatment"].values[:n], Y_norm[:n].ravel()
        dim = y[a == 1].mean() - y[a == 0].mean()  # naive (confounded) baseline
        rows.append(dict(
            setting=setting, dataset=Path(f).stem, n=len(df), n_x=len(x_cols),
            pehe=compute_pehe(cate, ite[n:]), ate_err=compute_ate_error(cate, ite[n:]),
            naive_ate_err=abs(dim - ite[n:].mean()), true_ate=ite[n:].mean(), sec=dt,
            has_uncertainty="gmm_sigma" in res,
        ))

out = pd.DataFrame(rows)
pd.set_option("display.width", 200)
print(out.round(3).to_string(index=False))
print("\nMeans by setting:")
print(out.groupby("setting")[["pehe", "ate_err", "naive_ate_err", "sec"]].mean().round(3))
