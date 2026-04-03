"""
Generate all paper figures for the CBPI (Cohort-Based Postponement Index) paper.
Methodology mirrors postponement_analysis.R but runs in Python.
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from matplotlib.lines import Line2D
import warnings
warnings.filterwarnings("ignore")

DATA_PATH = "."
plt.rcParams.update({
    "font.family": "serif",
    "font.size": 10,
    "axes.titlesize": 11,
    "axes.labelsize": 10,
    "legend.fontsize": 9,
    "figure.dpi": 150,
})

# ── 1. LOAD DATA ──────────────────────────────────────────────────────────────
print("Loading data...")
asfr_raw = pd.read_csv(f"{DATA_PATH}/fertility-rate-by-age-group.csv")
asfr_raw.columns = [c.lower() for c in asfr_raw.columns]
asfr_raw = asfr_raw.rename(columns={"entity": "country"})

mac_raw = pd.read_csv(f"{DATA_PATH}/period-average-age-of-mothers.csv")
mac_raw.columns = ["country", "code", "year", "mean_age"]

# Keep only 3-letter ISO country codes (drop regional aggregates)
asfr_raw = asfr_raw[asfr_raw["code"].notna() & (asfr_raw["code"].str.len() == 3)].copy()

# ── 2. RESHAPE AND CONVERT ────────────────────────────────────────────────────
age_cols = [c for c in asfr_raw.columns if c.startswith("f") and c[1:].isdigit()]
asfr_long = asfr_raw.melt(
    id_vars=["country", "code", "year"],
    value_vars=age_cols,
    var_name="age_str",
    value_name="rate_per_1000"
)
asfr_long["age"] = asfr_long["age_str"].str[1:].astype(int)
asfr_long["tfr_contrib"] = asfr_long["rate_per_1000"] * 5 / 1000

# Focus on ages 22–42 (5-year groups: 20-24 … 40-44)
asfr_long = asfr_long[asfr_long["age"].between(22, 42)].copy()

# ── 3. COHORT ALIGNMENT ───────────────────────────────────────────────────────
asfr_long["birth_year"] = asfr_long["year"] - asfr_long["age"]
asfr_long = asfr_long[asfr_long["birth_year"] % 5 == 1].copy()

# ── 4. COHORT TOTAL FERTILITY F(c) ────────────────────────────────────────────
cohort_totals = (
    asfr_long.groupby(["country", "code", "birth_year"])
    .agg(F_cohort=("tfr_contrib", "sum"), n_ages=("tfr_contrib", "count"))
    .reset_index()
)
asfr_long = asfr_long.merge(cohort_totals, on=["country", "code", "birth_year"])

# ── 5. REFERENCE AGE DISTRIBUTION FROM 1951 BIRTH COHORT ─────────────────────
ref = asfr_long[asfr_long["birth_year"] == 1951].copy()
ref["S0"] = ref["tfr_contrib"] / ref["F_cohort"]
ref = ref[["country", "code", "age", "S0"]]

# Only keep countries with complete 1951 reference (all 5 age groups)
ref_complete = ref.groupby(["country", "code"])["age"].count()
ref_complete = ref_complete[ref_complete >= 5].reset_index()[["country", "code"]]

asfr_long = asfr_long.merge(ref_complete, on=["country", "code"])
asfr_long = asfr_long.merge(ref, on=["country", "code", "age"], how="left")

# ── 6. COUNTERFACTUAL BIRTHS ──────────────────────────────────────────────────
asfr_long["cf_births"] = asfr_long["F_cohort"] * asfr_long["S0"]

# ── 7. COLLAPSE TO PERIOD YEAR ────────────────────────────────────────────────
# Mark cohorts that have ALL 5 age-group observations (complete cohort fertility).
# A cohort born in year c is observed at ages 22,27,32,37,42 in years c+22…c+42.
# We need all 5 ages present in the data before we trust F_cohort.
cohort_age_counts = (
    asfr_long.groupby(["country", "code", "birth_year"])["age"]
    .count()
    .reset_index(name="n_ages_cohort")
)
complete_cohorts = cohort_age_counts[cohort_age_counts["n_ages_cohort"] >= 5][
    ["country", "code", "birth_year"]
]
# Keep only rows belonging to complete cohorts
asfr_long_complete = asfr_long.merge(
    complete_cohorts, on=["country", "code", "birth_year"]
)

# Collapse to period year — a year is valid only if ALL 5 age groups in that
# year come from complete cohorts (n_complete_ages == 5)
tfr_ts = (
    asfr_long_complete.groupby(["country", "code", "year"])
    .agg(
        tfr_cf=("cf_births", "sum"),
        n_complete_ages=("tfr_contrib", "count"),
    )
    .reset_index()
)
tfr_ts = tfr_ts[tfr_ts["n_complete_ages"] >= 5].copy()

# Merge actual TFR (from the raw long data — not filtered by cohort completeness)
actual_tfr = (
    asfr_long.groupby(["country", "code", "year"])["tfr_contrib"]
    .sum()
    .reset_index(name="tfr_actual")
)
tfr_ts = tfr_ts.merge(actual_tfr, on=["country", "code", "year"])
tfr_ts["timing_effect"] = tfr_ts["tfr_actual"] - tfr_ts["tfr_cf"]

print(f"Countries with data: {tfr_ts['country'].nunique()}")

# ── 8. COUNTRY-LEVEL SUMMARY (early 1970s vs 2015–21) ────────────────────────
# Use early 1970s (start of postponement) vs. late 1990s–early 2000s
# (latest window with fully complete cohort fertility for all age groups).
# Data through 2023 means the youngest cohort at age 42 in 2003 was born 1961;
# that cohort completed all 5 ages by 2003, so 2003 is roughly our safe endpoint.
early = (
    tfr_ts[tfr_ts["year"].between(1973, 1978)]
    .groupby(["country", "code"])
    .agg(tfr_early=("tfr_actual", "mean"),
         ctfr_early=("tfr_cf", "mean"),
         te_early=("timing_effect", "mean"))
    .reset_index()
)
recent = (
    tfr_ts[tfr_ts["year"].between(1998, 2008)]
    .groupby(["country", "code"])
    .agg(tfr_recent=("tfr_actual", "mean"),
         ctfr_recent=("tfr_cf", "mean"),
         te_recent=("timing_effect", "mean"))
    .reset_index()
)
summary = early.merge(recent, on=["country", "code"])
summary["delta_tfr"]      = summary["tfr_early"]  - summary["tfr_recent"]
summary["delta_ctfr"]     = summary["ctfr_early"] - summary["ctfr_recent"]
summary["timing_contrib"] = summary["te_early"]   - summary["te_recent"]
# Require delta_tfr > 0.3 to avoid noisy CBPI from tiny denominators
summary["cbpi"] = np.where(
    summary["delta_tfr"] > 0.30,
    summary["timing_contrib"] / summary["delta_tfr"],
    np.nan
)

# Region lookup
region_map = {
    "CHN":"East Asia & Pacific","JPN":"East Asia & Pacific","KOR":"East Asia & Pacific",
    "SGP":"East Asia & Pacific","AUS":"East Asia & Pacific","NZL":"East Asia & Pacific",
    "VNM":"East Asia & Pacific","THA":"East Asia & Pacific","IDN":"East Asia & Pacific",
    "PHL":"East Asia & Pacific","MYS":"East Asia & Pacific","MNG":"East Asia & Pacific",
    "KHM":"East Asia & Pacific","MMR":"East Asia & Pacific","TWN":"East Asia & Pacific",
    "DEU":"Europe","FRA":"Europe","GBR":"Europe","ITA":"Europe","ESP":"Europe",
    "PRT":"Europe","BEL":"Europe","NLD":"Europe","AUT":"Europe","CHE":"Europe",
    "SWE":"Europe","NOR":"Europe","DNK":"Europe","FIN":"Europe","IRL":"Europe",
    "GRC":"Europe","POL":"Europe","CZE":"Europe","HUN":"Europe","SVK":"Europe",
    "SVN":"Europe","HRV":"Europe","BGR":"Europe","ROU":"Europe","LTU":"Europe",
    "LVA":"Europe","EST":"Europe","BLR":"Europe","UKR":"Europe","MDA":"Europe",
    "SRB":"Europe","BIH":"Europe","MKD":"Europe","ALB":"Europe","LUX":"Europe",
    "ISL":"Europe","RUS":"Europe","MNE":"Europe",
    "BRA":"Latin America","MEX":"Latin America","COL":"Latin America","ARG":"Latin America",
    "PER":"Latin America","VEN":"Latin America","CHL":"Latin America","ECU":"Latin America",
    "BOL":"Latin America","PRY":"Latin America","URY":"Latin America","PAN":"Latin America",
    "CRI":"Latin America","GTM":"Latin America","HND":"Latin America","SLV":"Latin America",
    "NIC":"Latin America","DOM":"Latin America","CUB":"Latin America","HTI":"Latin America",
    "JAM":"Latin America","TTO":"Latin America","GUY":"Latin America","SUR":"Latin America",
    "USA":"North America","CAN":"North America",
    "IND":"South Asia","PAK":"South Asia","BGD":"South Asia","NPL":"South Asia",
    "LKA":"South Asia","AFG":"South Asia","MDV":"South Asia","BTN":"South Asia",
    "TUR":"Middle East & N. Africa","IRN":"Middle East & N. Africa",
    "EGY":"Middle East & N. Africa","DZA":"Middle East & N. Africa",
    "MAR":"Middle East & N. Africa","TUN":"Middle East & N. Africa",
    "SAU":"Middle East & N. Africa","IRQ":"Middle East & N. Africa",
    "SYR":"Middle East & N. Africa","JOR":"Middle East & N. Africa",
    "LBN":"Middle East & N. Africa","ISR":"Middle East & N. Africa",
    "YEM":"Middle East & N. Africa","OMN":"Middle East & N. Africa",
    "ARE":"Middle East & N. Africa","KWT":"Middle East & N. Africa",
    "QAT":"Middle East & N. Africa","BHR":"Middle East & N. Africa",
    "LBY":"Middle East & N. Africa",
}
# More complete region lookup using UN M.49 classification
# Add extra codes for territories and countries not in the core lookup
extra_regions = {
    # Eastern Europe / Central Asia
    "GEO":"Europe","ARM":"Europe","AZE":"Europe",
    "KAZ":"Central Asia","KGZ":"Central Asia","TKM":"Central Asia",
    "TJK":"Central Asia","UZB":"Central Asia",
    # Pacific territories and small island states
    "TON":"East Asia & Pacific","WSM":"East Asia & Pacific",
    "VUT":"East Asia & Pacific","SLB":"East Asia & Pacific",
    "KIR":"East Asia & Pacific","FSM":"East Asia & Pacific",
    "MHL":"East Asia & Pacific","PLW":"East Asia & Pacific",
    "NRU":"East Asia & Pacific","TUV":"East Asia & Pacific",
    "COK":"East Asia & Pacific","NIU":"East Asia & Pacific",
    "TKL":"East Asia & Pacific","WLF":"East Asia & Pacific",
    # European micro-states
    "SMR":"Europe","AND":"Europe","MCO":"Europe","LIE":"Europe",
    "GGY":"Europe","JEY":"Europe","IMN":"Europe","GIB":"Europe",
    "FRO":"Europe","ALA":"Europe","SJM":"Europe",
    # Caribbean
    "CUW":"Latin America","ABW":"Latin America","SXM":"Latin America",
    "BLM":"Latin America","MAF":"Latin America","AIA":"Latin America",
    "MSR":"Latin America","VGB":"Latin America","VIR":"Latin America",
    "PRI":"Latin America","BMU":"North America",
    "CYM":"Latin America","TCA":"Latin America","BHS":"Latin America",
    "BRB":"Latin America","ATG":"Latin America","DMA":"Latin America",
    "GRD":"Latin America","KNA":"Latin America","LCA":"Latin America",
    "VCT":"Latin America",
    # East Asia additions
    "HKG":"East Asia & Pacific","MAC":"East Asia & Pacific",
    "PRK":"East Asia & Pacific","TLS":"East Asia & Pacific",
    "BRN":"East Asia & Pacific",
    # Middle East additions
    "PSE":"Middle East & N. Africa","GEQ":"Sub-Saharan Africa",
    # South Asia additions
    "SRI":"South Asia",
    # Africa (ensure coverage)
    "ERI":"Sub-Saharan Africa","SOM":"Sub-Saharan Africa",
    "SSD":"Sub-Saharan Africa","CAF":"Sub-Saharan Africa",
}
region_map.update(extra_regions)

def assign_region(code):
    if pd.isna(code):
        return "Other"
    r = region_map.get(code)
    if r:
        return r
    # Heuristic: use code ranges for unlisted African countries
    african_prefix = ["NG","GH","SN","CI","CM","ZA","KE","ET","TZ",
                      "UG","MZ","ZW","ZM","AO","CD","BJ","BF","ML",
                      "NE","TD","SO","BI","RW","MG","MW","SL","LR",
                      "GN","GM","GW","CV","ST","SC","MU","KM","DJ",
                      "ER","CF","CG","GA","GQ","SZ","LS","BW","NA",
                      "TG","BN","SS","SD"]
    if any(code.startswith(p) for p in african_prefix):
        return "Sub-Saharan Africa"
    return "Other"

summary["region"] = summary["code"].apply(assign_region)
tfr_ts["region"] = tfr_ts["code"].apply(assign_region)

REGION_COLORS = {
    "East Asia & Pacific":    "#1f78b4",
    "Europe":                 "#e31a1c",
    "Latin America":          "#33a02c",
    "North America":          "#ff7f00",
    "South Asia":             "#6a3d9a",
    "Middle East & N. Africa":"#b15928",
    "Sub-Saharan Africa":     "#888888",
    "Central Asia":           "#a65628",
    "Other":                  "#cccccc",
}

summary.to_csv(f"{DATA_PATH}/postponement_summary.csv", index=False)
print("Summary CSV saved.")

# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 1: Actual vs. Counterfactual TFR — illustrative countries
# ─────────────────────────────────────────────────────────────────────────────
print("Generating Figure 1...")
illustrative = ["Chile", "Spain", "South Korea", "Sweden",
                "Israel", "United States", "Japan", "Brazil"]

fig, axes = plt.subplots(2, 4, figsize=(14, 7))
axes = axes.flatten()

for i, ctry in enumerate(illustrative):
    ax = axes[i]
    d = tfr_ts[tfr_ts["country"] == ctry].sort_values("year")
    if d.empty:
        ax.set_visible(False)
        continue
    ax.plot(d["year"], d["tfr_actual"], color="#1f4e79", lw=1.8, label="Actual TFR")
    ax.plot(d["year"], d["tfr_cf"],     color="#c00000", lw=1.8, ls="--",
            label="Counterfactual TFR\n(1951 timing)")
    ax.fill_between(d["year"], d["tfr_actual"], d["tfr_cf"],
                    where=d["tfr_actual"] < d["tfr_cf"],
                    alpha=0.12, color="#c00000", label="_")
    ax.set_title(ctry, fontweight="bold")
    ax.set_xlabel("Year")
    ax.set_ylabel("TFR")
    ax.set_xlim(1970, 2010)
    ax.grid(axis="y", alpha=0.3)
    ax.tick_params(axis="x", rotation=30)

# Shared legend
handles = [
    Line2D([0],[0], color="#1f4e79", lw=2, label="Actual TFR"),
    Line2D([0],[0], color="#c00000", lw=2, ls="--", label="Counterfactual TFR (1951 timing)"),
]
fig.legend(handles=handles, loc="lower center", ncol=2, fontsize=10,
           bbox_to_anchor=(0.5, -0.01))
fig.suptitle("Actual and Counterfactual TFR by Year\n"
             "(shaded area = postponement depresses period TFR)",
             fontsize=12, fontweight="bold", y=1.01)
plt.tight_layout()
fig.savefig(f"{DATA_PATH}/fig1_actual_vs_cf.pdf", bbox_inches="tight")
fig.savefig(f"{DATA_PATH}/fig1_actual_vs_cf.png", bbox_inches="tight", dpi=200)
plt.close()
print("  fig1 done.")

# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 2: Cohort age profiles (births by birth year and age group)
# ─────────────────────────────────────────────────────────────────────────────
print("Generating Figure 2...")
cohort_countries = ["Israel", "Chile", "South Korea", "Spain"]
age_colors = {22:"#1f78b4", 27:"#33a02c", 32:"#e31a1c", 37:"#ff7f00", 42:"#6a3d9a"}
age_labels  = {22:"Age 20–24", 27:"Age 25–29", 32:"Age 30–34",
               37:"Age 35–39", 42:"Age 40–44"}

fig, axes = plt.subplots(2, 2, figsize=(12, 8))
axes = axes.flatten()

for i, ctry in enumerate(cohort_countries):
    ax = axes[i]
    d = asfr_long[
        (asfr_long["country"] == ctry) &
        (asfr_long["age"].between(22, 42)) &
        (asfr_long["birth_year"] >= 1910)
    ].sort_values("birth_year")
    for age_val in [22, 27, 32, 37, 42]:
        sub = d[d["age"] == age_val]
        ax.plot(sub["birth_year"], sub["tfr_contrib"],
                color=age_colors[age_val], lw=1.5,
                label=age_labels[age_val])
    ax.set_title(ctry, fontweight="bold")
    ax.set_xlabel("Birth year of mother")
    ax.set_ylabel("TFR contribution (births per woman)")
    ax.grid(axis="y", alpha=0.3)

handles2 = [Line2D([0],[0], color=age_colors[a], lw=2, label=age_labels[a])
            for a in [22, 27, 32, 37, 42]]
fig.legend(handles=handles2, loc="lower center", ncol=5,
           bbox_to_anchor=(0.5, -0.02))
fig.suptitle("Births by Birth Year and Age Group\n"
             "(postponement = younger ages fall, older ages rise)",
             fontsize=12, fontweight="bold")
plt.tight_layout()
fig.savefig(f"{DATA_PATH}/fig2_cohort_age_profiles.pdf", bbox_inches="tight")
fig.savefig(f"{DATA_PATH}/fig2_cohort_age_profiles.png", bbox_inches="tight", dpi=200)
plt.close()
print("  fig2 done.")

# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 3: Cross-country scatter — CBPI vs total TFR decline
# ─────────────────────────────────────────────────────────────────────────────
print("Generating Figure 3...")
scatter = summary[
    summary["cbpi"].notna() &
    summary["cbpi"].between(-0.5, 2.0) &
    (summary["delta_tfr"] > 0.3) &
    (summary["region"] != "Other")
].copy()

# Labels for notable countries
label_countries = {
    "South Korea","Japan","Spain","Italy","Greece","Sweden","France",
    "Israel","United States","Chile","Brazil","India","Nigeria",
    "Germany","South Africa","China","Russia","Canada","Portugal",
    "Poland","Hungary","Iran","Turkey","Mexico","Argentina",
}

fig, ax = plt.subplots(figsize=(11, 7))
for region, grp in scatter.groupby("region"):
    ax.scatter(grp["delta_tfr"], grp["cbpi"],
               color=REGION_COLORS.get(region, "#888888"),
               s=45, alpha=0.82, label=region, zorder=3)

# OLS trend line
from numpy.polynomial.polynomial import polyfit as npfit
x = scatter["delta_tfr"].values
y = scatter["cbpi"].values
mask = np.isfinite(x) & np.isfinite(y)
coef = np.polyfit(x[mask], y[mask], 1)
xline = np.linspace(x[mask].min(), x[mask].max(), 200)
ax.plot(xline, np.polyval(coef, xline), color="black",
        ls="--", lw=1.2, zorder=2, label="OLS fit")

# Annotate selected countries
for _, row in scatter.iterrows():
    if row["country"] in label_countries:
        ax.annotate(row["country"],
                    xy=(row["delta_tfr"], row["cbpi"]),
                    xytext=(4, 3), textcoords="offset points",
                    fontsize=7.5, color="#333333")

ax.axhline(0,   color="grey", lw=0.8, ls=":")
ax.axhline(0.5, color="grey", lw=0.8, ls=":")
ax.axhline(1.0, color="grey", lw=0.8, ls=":")
ax.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1, decimals=0))
ax.set_xlabel("Total TFR Decline (early 1970s to 2015–21)", fontsize=11)
ax.set_ylabel("Share of Decline Due to Postponement (CBPI)", fontsize=11)
ax.set_title("Cohort-Based Postponement Index Across Countries", fontsize=13,
             fontweight="bold")
ax.legend(loc="upper right", fontsize=8, framealpha=0.9)
ax.grid(alpha=0.25)
plt.tight_layout()
fig.savefig(f"{DATA_PATH}/fig3_postponement_scatter.pdf", bbox_inches="tight")
fig.savefig(f"{DATA_PATH}/fig3_postponement_scatter.png", bbox_inches="tight", dpi=200)
plt.close()
print("  fig3 done.")

# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 4: Regional average timing effect over time
# ─────────────────────────────────────────────────────────────────────────────
print("Generating Figure 4...")
regions_to_plot = ["East Asia & Pacific", "Europe", "Latin America",
                   "North America", "Sub-Saharan Africa"]

region_ts = (
    tfr_ts[tfr_ts["region"].isin(regions_to_plot)]
    .groupby(["region", "year"])["timing_effect"]
    .mean()
    .reset_index()
)

fig, ax = plt.subplots(figsize=(10, 6))
ax.axhline(0, color="black", lw=0.9)
for region in regions_to_plot:
    d = region_ts[region_ts["region"] == region].sort_values("year")
    ax.plot(d["year"], d["timing_effect"],
            color=REGION_COLORS[region], lw=2.0, label=region)

ax.set_xlabel("Year", fontsize=11)
ax.set_ylabel("Timing Effect on TFR (actual − counterfactual)", fontsize=11)
ax.set_title("Regional Average Timing (Postponement) Effect on Period TFR\n"
             "Negative = postponement depresses period TFR below quantum",
             fontsize=12, fontweight="bold")
ax.legend(fontsize=9)
ax.grid(alpha=0.25)
ax.set_xlim(1970, 2010)
plt.tight_layout()
fig.savefig(f"{DATA_PATH}/fig4_timing_effect_regions.pdf", bbox_inches="tight")
fig.savefig(f"{DATA_PATH}/fig4_timing_effect_regions.png", bbox_inches="tight", dpi=200)
plt.close()
print("  fig4 done.")

# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 5: Bar chart — top and bottom 15 by CBPI
# ─────────────────────────────────────────────────────────────────────────────
print("Generating Figure 5...")
bar_df = summary[
    summary["cbpi"].notna() &
    summary["cbpi"].between(-0.5, 2.0) &
    (summary["delta_tfr"] > 0.3) &
    (summary["region"] != "Other")
].copy().sort_values("cbpi")

top15    = bar_df.nlargest(15, "cbpi")
bottom15 = bar_df.nsmallest(15, "cbpi")
bar_data = pd.concat([bottom15, top15]).drop_duplicates(subset="country")
bar_data = bar_data.sort_values("cbpi")
bar_data["color"] = bar_data["cbpi"].apply(
    lambda v: "#1f4e79" if v >= bar_data["cbpi"].nlargest(15).min() else "#c00000"
)

fig, ax = plt.subplots(figsize=(9, 9))
bars = ax.barh(bar_data["country"], bar_data["cbpi"],
               color=bar_data["color"], alpha=0.85)
ax.axvline(0, color="black", lw=0.8)
ax.xaxis.set_major_formatter(mtick.PercentFormatter(xmax=1, decimals=0))
ax.set_xlabel("Share of TFR Decline Attributed to Postponement (CBPI)", fontsize=10)
ax.set_title("Countries by Cohort-Based Postponement Index\n"
             "Top and bottom 15 (minimum TFR decline = 0.1)",
             fontsize=12, fontweight="bold")
legend_handles = [
    plt.Rectangle((0,0),1,1, color="#1f4e79", alpha=0.85, label="High postponement share"),
    plt.Rectangle((0,0),1,1, color="#c00000", alpha=0.85, label="Low postponement share"),
]
ax.legend(handles=legend_handles, fontsize=9)
ax.grid(axis="x", alpha=0.25)
plt.tight_layout()
fig.savefig(f"{DATA_PATH}/fig5_postponement_bar.pdf", bbox_inches="tight")
fig.savefig(f"{DATA_PATH}/fig5_postponement_bar.png", bbox_inches="tight", dpi=200)
plt.close()
print("  fig5 done.")

# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 6: Israel deep-dive
# ─────────────────────────────────────────────────────────────────────────────
print("Generating Figure 6...")
isr = tfr_ts[tfr_ts["country"] == "Israel"].sort_values("year")

fig, axes = plt.subplots(1, 2, figsize=(13, 5))

# Panel A: actual vs counterfactual
ax = axes[0]
ax.plot(isr["year"], isr["tfr_actual"], color="#1f4e79", lw=2, label="Actual TFR")
ax.plot(isr["year"], isr["tfr_cf"],     color="#c00000", lw=2, ls="--",
        label="Counterfactual TFR (1951 timing)")
ax.fill_between(isr["year"], isr["tfr_actual"], isr["tfr_cf"],
                where=isr["tfr_actual"] < isr["tfr_cf"],
                alpha=0.15, color="#c00000")
ax.set_title("Actual vs. Counterfactual TFR", fontweight="bold")
ax.set_xlabel("Year"); ax.set_ylabel("TFR")
ax.legend(); ax.grid(alpha=0.25)
ax.set_xlim(1970, 2010)

# Panel B: timing effect alone
ax = axes[1]
ax.axhline(0, color="black", lw=0.8)
ax.fill_between(isr["year"], isr["timing_effect"], 0,
                where=isr["timing_effect"] < 0, alpha=0.4,
                color="#c00000", label="Postponement depresses TFR")
ax.fill_between(isr["year"], isr["timing_effect"], 0,
                where=isr["timing_effect"] >= 0, alpha=0.4,
                color="#33a02c", label="Early timing boosts TFR")
ax.plot(isr["year"], isr["timing_effect"], color="#333333", lw=1.5)
ax.set_title("Timing Effect on TFR\n(Actual − Counterfactual)", fontweight="bold")
ax.set_xlabel("Year"); ax.set_ylabel("TFR points")
ax.legend(fontsize=8); ax.grid(alpha=0.25)
ax.set_xlim(1970, 2010)

fig.suptitle("Israel: Decomposing TFR into Quantum and Timing Components",
             fontsize=13, fontweight="bold")
plt.tight_layout()
fig.savefig(f"{DATA_PATH}/fig6_israel_decomposition.pdf", bbox_inches="tight")
fig.savefig(f"{DATA_PATH}/fig6_israel_decomposition.png", bbox_inches="tight", dpi=200)
plt.close()
print("  fig6 done.")

# ─────────────────────────────────────────────────────────────────────────────
# PRINT TOP/BOTTOM TABLE
# ─────────────────────────────────────────────────────────────────────────────
valid = summary[
    summary["cbpi"].notna() &
    summary["cbpi"].between(-0.5, 2.0) &
    (summary["delta_tfr"] > 0.3) &
    (summary["region"] != "Other")
].copy()

print("\n--- Top 15 countries by CBPI (highest postponement share) ---")
print(
    valid.nlargest(15, "cbpi")
    [["country","region","tfr_early","tfr_recent","delta_tfr","ctfr_recent","cbpi"]]
    .to_string(index=False, float_format="%.2f")
)

print("\n--- Bottom 15 countries by CBPI (lowest postponement share) ---")
print(
    valid.nsmallest(15, "cbpi")
    [["country","region","tfr_early","tfr_recent","delta_tfr","ctfr_recent","cbpi"]]
    .to_string(index=False, float_format="%.2f")
)

print("\nAll figures saved to:", DATA_PATH)
