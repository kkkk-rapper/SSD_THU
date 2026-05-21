#!/usr/bin/env python3
"""Plot the 122 keV energy spectrum from hits.csv (per-event deposited energy)."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---- per-event total deposited energy ------------------------------------
d = np.genfromtxt("hits.csv", delimiter=",", names=True)
ev = d["event"].astype(np.int64)
uev, inv = np.unique(ev, return_inverse=True)
etot = np.bincount(inv, weights=d["edep_keV"])      # keV, one entry per event
N = len(etot)
print(f"{len(ev)} hits -> {N} events")

E0 = 122.0
compton_edge = 2 * E0**2 / (511.0 + 2 * E0)          # 39.4 keV
escape = E0 - 9.886                                  # Ge Kalpha escape ~112 keV
peak = np.sum((etot > 121.0) & (etot < 123.0))
print(f"photopeak {peak} ({100*peak/N:.1f}%)  Compton edge {compton_edge:.1f} keV")

bins = np.arange(0, 126.0, 0.25)
fig, (axL, axR) = plt.subplots(1, 2, figsize=(15, 6))
fig.suptitle(f"sim_project — 122 keV gamma energy spectrum in HPGe "
             f"({N:,} events with energy deposit)", fontsize=13, fontweight="bold")

def features(ax):
    ax.axvline(E0, color="red", ls="--", lw=1, alpha=.7)
    ax.axvline(compton_edge, color="green", ls=":", lw=1, alpha=.7)
    ax.axvline(escape, color="orange", ls=":", lw=1, alpha=.7)

# ---- (A) raw Geant4 deposited-energy spectrum ----------------------------
axL.hist(etot, bins=bins, histtype="step", color="navy", lw=1.2)
axL.set_yscale("log")
features(axL)
axL.annotate("122 keV\nfull-energy peak", xy=(E0, N*0.5), xytext=(88, N*0.35),
             color="red", fontsize=9, ha="center",
             arrowprops=dict(arrowstyle="->", color="red"))
axL.annotate(f"Ge X-ray escape\n~{escape:.0f} keV", xy=(escape, 300),
             xytext=(75, 4000), color="orange", fontsize=9, ha="center",
             arrowprops=dict(arrowstyle="->", color="orange"))
axL.annotate(f"Compton edge\n{compton_edge:.0f} keV", xy=(compton_edge, 200),
             xytext=(52, 3000), color="green", fontsize=9, ha="center",
             arrowprops=dict(arrowstyle="->", color="green"))
axL.text(20, 60, "Compton\ncontinuum", color="gray", fontsize=9, ha="center")
axL.set_xlabel("Deposited energy (keV)"); axL.set_ylabel("Counts / 0.25 keV")
axL.set_title("Raw Geant4 spectrum (true deposited energy)")
axL.set_xlim(0, 126); axL.grid(alpha=.25, which="both")

# ---- (B) with a typical HPGe energy resolution (1 keV FWHM) --------------
fwhm = 1.0
sigma = fwhm / 2.355
rng = np.random.default_rng(1)
meas = etot + rng.normal(0.0, sigma, N)
axR.hist(meas, bins=bins, histtype="stepfilled", color="purple", alpha=.8)
axR.set_yscale("log")
features(axR)
axR.set_xlabel("Energy (keV)"); axR.set_ylabel("Counts / 0.25 keV")
axR.set_title(f"With detector resolution applied ({fwhm:.1f} keV FWHM)")
axR.set_xlim(0, 126); axR.grid(alpha=.25, which="both")
axR.text(0.97, 0.95,
         f"photopeak: {100*peak/N:.0f}% of events\n"
         f"escape peak: ~{escape:.0f} keV\n"
         f"Compton edge: {compton_edge:.0f} keV",
         transform=axR.transAxes, ha="right", va="top", fontsize=9,
         bbox=dict(boxstyle="round", fc="white", ec="gray", alpha=.9))

fig.tight_layout(rect=(0, 0, 1, 0.95))
fig.savefig("spectrum.png", dpi=120)
print("wrote spectrum.png")
