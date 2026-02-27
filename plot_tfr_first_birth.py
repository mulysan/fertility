"""
Scatter plot: TFR vs Mean Age at Childbearing (proxy for avg age at first birth) by country.

From the age-specific fertility rates (ASFR) we derive:
  TFR  = 5 * sum(ASFR_i) / 1000          (each group spans 5 years, rates are per 1 000)
  MAC  = sum(ASFR_i * midpoint_i) / sum(ASFR_i)   (mean age at childbearing)

MAC is the standard demographic proxy for mean age at first birth when
parity-specific data are unavailable.
"""

import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from scipy import stats

# ── 1. Load data ──────────────────────────────────────────────────────────────
df = pd.read_csv("fertility-rate-by-age-group.csv")

# Age-group midpoints (5-year groups: 10-14, 15-19, …, 50-54)
midpoints = np.array([12, 17, 22, 27, 32, 37, 42, 47, 52])
cols      = ["f12", "f17", "f22", "f27", "f32", "f37", "f42", "f47", "f52"]

# ── 2. Keep only ISO-3 countries, most recent year ────────────────────────────
iso3 = re.compile(r"^[A-Z]{3}$")
df = df[df["Code"].apply(lambda c: bool(iso3.match(str(c))))]

# Most recent year with data for each country
latest = df.sort_values("Year").groupby("Entity").last().reset_index()

# ── 3. Compute TFR and mean age at childbearing (MAC) ────────────────────────
asfr = latest[cols].values                    # births per 1 000 women per year

total_rate = asfr.sum(axis=1)                 # sum of ASFRs
tfr   = 5 * total_rate / 1000                 # TFR (children per woman)
mac   = (asfr * midpoints).sum(axis=1) / total_rate  # MAC (years)

latest = latest.assign(TFR=tfr, MAC=mac)

# Drop rows with missing / degenerate values
latest = latest.dropna(subset=["TFR", "MAC"])
latest = latest[latest["TFR"] > 0]

# ── 4. Broad regional colouring (UN geoscheme, hand-coded by ISO-3) ──────────
# Assign each country to a broad region for colour
africa = {
    "DZA","AGO","BEN","BWA","BFA","BDI","CPV","CMR","CAF","TCD","COM","COD",
    "COG","CIV","DJI","EGY","GNQ","ERI","SWZ","ETH","GAB","GMB","GHA","GIN",
    "GNB","KEN","LSO","LBR","LBY","MDG","MWI","MLI","MRT","MUS","MAR","MOZ",
    "NAM","NER","NGA","RWA","STP","SEN","SLE","SOM","ZAF","SSD","SDN","TZA",
    "TGO","TUN","UGA","ZMB","ZWE","REU","MYT","SHN",
}
europe = {
    "ALB","AND","AUT","BLR","BEL","BIH","BGR","HRV","CYP","CZE","DNK","EST",
    "FIN","FRA","DEU","GRC","HUN","ISL","IRL","ITA","XKX","LVA","LIE","LTU",
    "LUX","MLT","MDA","MCO","MNE","NLD","MKD","NOR","POL","PRT","ROU","RUS",
    "SMR","SRB","SVK","SVN","ESP","SWE","CHE","UKR","GBR","VAT",
}
americas = {
    "ATG","ARG","BHS","BRB","BLZ","BOL","BRA","CAN","CHL","COL","CRI","CUB",
    "DMA","DOM","ECU","SLV","GRD","GTM","GUY","HTI","HND","JAM","MEX","NIC",
    "PAN","PRY","PER","KNA","LCA","VCT","SUR","TTO","USA","URY","VEN",
}
oceania = {
    "AUS","FJI","KIR","MHL","FSM","NRU","NZL","PLW","PNG","WSM","SLB","TON",
    "TUV","VUT",
}

def get_region(code):
    if code in africa:   return "Africa"
    if code in europe:   return "Europe"
    if code in americas: return "Americas"
    if code in oceania:  return "Oceania"
    return "Asia & Middle East"

latest["Region"] = latest["Code"].apply(get_region)

palette = {
    "Africa":           "#e06c1a",
    "Europe":           "#3a7eca",
    "Americas":         "#2ca02c",
    "Asia & Middle East": "#9467bd",
    "Oceania":          "#d62728",
}

# ── 5. OLS regression ─────────────────────────────────────────────────────────
slope, intercept, r, p, se = stats.linregress(latest["MAC"], latest["TFR"])
x_line = np.linspace(latest["MAC"].min() - 0.5, latest["MAC"].max() + 0.5, 200)
y_line = slope * x_line + intercept

# ── 6. Plot ───────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(11, 8))

for region, grp in latest.groupby("Region"):
    ax.scatter(grp["MAC"], grp["TFR"],
               color=palette[region], label=region,
               s=55, alpha=0.80, edgecolors="white", linewidths=0.4, zorder=3)

ax.plot(x_line, y_line, color="#333333", lw=1.8, ls="--", zorder=4,
        label=f"OLS fit  (r = {r:.2f}, p < 0.001)" if p < 0.001
              else f"OLS fit  (r = {r:.2f}, p = {p:.3f})")

# Label a selection of notable / extreme countries
label_codes = {
    "NER","TCD","MLI","SOM","AGO",          # high TFR
    "KOR","JPN","HKG","SGP","ITA","ESP",    # lowest TFR
    "SWE","NOR","FIN","ISL",                # high MAC
    "USA","GBR","FRA","BRA","CHN","IND",    # large / interesting
    "NGA","ETH","COD",                      # populous Africa
}
offsets = {  # (dx, dy) in data units — tweak to reduce overlap
    "NER": (-0.6, 0.15), "TCD": (0.15, -0.35), "MLI": (0.15, 0.1),
    "SOM": (-1.5, 0.1),  "AGO": (0.15, 0.1),
    "KOR": (-1.4, 0.05), "JPN": (0.15, -0.3),
    "ITA": (-1.3, 0.1),  "ESP": (0.15, 0.1),
    "SWE": (0.15, 0.1),  "NOR": (-1.5, 0.1), "ISL": (-1.6, -0.3),
    "FIN": (-1.4, -0.3),
    "USA": (0.15, 0.1),  "GBR": (-1.4, 0.1), "FRA": (0.15, -0.3),
    "BRA": (0.15, 0.1),  "CHN": (-1.4, -0.3), "IND": (0.15, 0.1),
    "NGA": (0.15, 0.1),  "ETH": (-1.6, -0.3), "COD": (0.15, 0.1),
    "HKG": (0.15, 0.1),  "SGP": (0.15, -0.3),
}
for _, row in latest.iterrows():
    if row["Code"] not in label_codes:
        continue
    dx, dy = offsets.get(row["Code"], (0.15, 0.1))
    ax.annotate(
        row["Code"],
        xy=(row["MAC"], row["TFR"]),
        xytext=(row["MAC"] + dx, row["TFR"] + dy),
        fontsize=7.5, color="#222222",
        arrowprops=dict(arrowstyle="-", color="#aaaaaa", lw=0.7),
        path_effects=[pe.withStroke(linewidth=2, foreground="white")],
        zorder=5,
    )

ax.set_xlabel("Mean Age at Childbearing (years)\n"
              "[proxy for average age at first birth]", fontsize=12)
ax.set_ylabel("Total Fertility Rate (children per woman)", fontsize=12)
ax.set_title("TFR vs Average Age at First Birth by Country\n"
             f"(most recent available year per country, n = {len(latest)})",
             fontsize=13, fontweight="bold")

ax.legend(framealpha=0.9, fontsize=9, loc="upper right")
ax.grid(True, alpha=0.25, linestyle=":")
ax.tick_params(labelsize=10)

plt.tight_layout()
plt.savefig("tfr_vs_first_birth_age.png", dpi=150, bbox_inches="tight")
print(f"Saved tfr_vs_first_birth_age.png  (n={len(latest)}, r={r:.3f})")
