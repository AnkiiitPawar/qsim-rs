"""Plot docs/grover.csv (written by `cargo run --release --example grover_curve`).

usage: python scripts/plot_grover.py [csv] [output_png]
"""

import csv
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

src = Path(sys.argv[1] if len(sys.argv) > 1 else "docs/grover.csv")
out = Path(sys.argv[2] if len(sys.argv) > 2 else "docs/grover.png")

curve, dist = [], []
with open(src) as fh:
    for row in csv.DictReader(fh):
        if row["kind"] == "curve":
            curve.append((int(row["x"]), float(row["simulated"]), float(row["theory"])))
        else:
            dist.append((int(row["x"]), float(row["simulated"])))

r, sim, theory = zip(*curve)
best = max(curve, key=lambda c: c[1])

fig, ax = plt.subplots(1, 2, figsize=(13, 4.4))
ax[0].plot(r, theory, "-", color="0.6", lw=2, label=r"theory $\sin^2((2r+1)\theta)$")
ax[0].plot(r, sim, "o", ms=4, label="simulated")
ax[0].axvline(best[0], ls="--", color="k", lw=0.8)
ax[0].set(xlabel="Grover iterations r", ylabel="P(marked item)",
          title=f"Grover search over N = 256: peak {best[1]:.4f} at r = {best[0]}")
ax[0].legend(loc="upper center", bbox_to_anchor=(0.68, 0.99), fontsize=9)

k, p = zip(*dist)
ax[1].bar(k, p, width=1.0)
ax[1].set(xlabel="basis state", ylabel="probability", yscale="log", ylim=(1e-8, 2.0),
          title=f"Distribution after r = {best[0]}: one spike, 255 states near zero")

fig.tight_layout()
out.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(out, dpi=110)
print(f"wrote {out}")
