"""
================================================================
Day 01 — 16S rRNA Alpha Diversity Analysis (REAL DATA)
Author  : Subhadip Jana
Dataset : peerj32 — LGG Probiotic vs Placebo intervention
          44 samples × 130 real gut taxa
          Source: microbiome R package (Lahti et al.)

Study Design:
  • Group: LGG probiotic (16 samples) vs Placebo (28 samples)
  • Time:  Before (1) vs After (2) probiotic intervention
  • Gender: Female (30) vs Male (14)

Research Question:
  Does LGG probiotic supplementation alter gut microbiome
  alpha diversity compared to placebo?
================================================================
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu, kruskal
from itertools import combinations
import warnings
warnings.filterwarnings("ignore")

# ─────────────────────────────────────────────────────────────
# SECTION 1: LOAD DATA
# ─────────────────────────────────────────────────────────────

print("🔬 Loading peerj32 dataset...")
otu_raw = pd.read_csv("data/otu_table.csv",  index_col=0)
meta    = pd.read_csv("data/metadata.csv",   index_col=0)

# Transpose OTU: rows = samples, cols = taxa
otu_df  = otu_raw.T
taxa    = otu_df.columns.tolist()

# Merge with metadata
df = otu_df.copy()
df["group"]  = meta["group"]
df["time"]   = meta["time"].astype(str)
df["gender"] = meta["gender"]
df["subject"]= meta["subject"]
df["group_time"] = df["group"] + "_T" + df["time"]

print(f"✅ {len(df)} samples × {len(taxa)} taxa")
print(f"   Groups    : {df['group'].value_counts().to_dict()}")
print(f"   Time pts  : {df['time'].value_counts().to_dict()}")
print(f"   Gender    : {df['gender'].value_counts().to_dict()}")

count_df = df[taxa].astype(float)

# ─────────────────────────────────────────────────────────────
# SECTION 2: ALPHA DIVERSITY FUNCTIONS
# ─────────────────────────────────────────────────────────────

def shannon(c):
    c = c[c > 0]; p = c / c.sum()
    return -np.sum(p * np.log(p))

def simpson(c):
    c = c[c > 0]; n = c.sum()
    return 1 - np.sum(c*(c-1))/(n*(n-1)) if n > 1 else 0

def observed_otus(c):
    return int(np.sum(c > 0))

def chao1(c):
    f1 = np.sum(c == 1); f2 = np.sum(c == 2)
    s  = observed_otus(c)
    return s + (f1**2)/(2*f2) if f2 > 0 else s + f1*(f1-1)/2

def evenness(c):
    h = shannon(c); s = observed_otus(c)
    return h / np.log(s) if s > 1 else 0

# ─────────────────────────────────────────────────────────────
# SECTION 3: COMPUTE METRICS
# ─────────────────────────────────────────────────────────────

print("\n📊 Computing alpha diversity metrics...")
records = []
for sid, row in count_df.iterrows():
    c = row.values
    records.append({
        "SampleID"        : sid,
        "Group"           : df.loc[sid, "group"],
        "Time"            : df.loc[sid, "time"],
        "Gender"          : df.loc[sid, "gender"],
        "Group_Time"      : df.loc[sid, "group_time"],
        "Shannon"         : round(shannon(c), 4),
        "Simpson"         : round(simpson(c), 4),
        "Observed_OTUs"   : observed_otus(c),
        "Chao1"           : round(chao1(c), 2),
        "Pielou_Evenness" : round(evenness(c), 4),
        "Total_Reads"     : round(c.sum(), 2),
    })

div_df = pd.DataFrame(records)
div_df.to_csv("outputs/alpha_diversity_results.csv", index=False)
print("✅ Saved → outputs/alpha_diversity_results.csv")

# ─────────────────────────────────────────────────────────────
# SECTION 4: STATISTICAL TESTS
# ─────────────────────────────────────────────────────────────

METRICS = ["Shannon", "Simpson", "Observed_OTUs", "Chao1", "Pielou_Evenness"]

print("\n" + "="*60)
print("STATISTICAL RESULTS")
print("="*60)

# Test 1: LGG vs Placebo
print("\n[1] LGG vs Placebo (Mann-Whitney U)")
stat_group = {}
for m in METRICS:
    d1 = div_df[div_df["Group"]=="LGG"][m].values
    d2 = div_df[div_df["Group"]=="Placebo"][m].values
    _, p = mannwhitneyu(d1, d2, alternative="two-sided")
    stat_group[m] = p
    sig = "***" if p<0.001 else "**" if p<0.01 else "*" if p<0.05 else "ns"
    print(f"  {m:20s}: p={p:.4f} {sig}")

# Test 2: Time 1 vs Time 2
print("\n[2] Time 1 vs Time 2 (Mann-Whitney U)")
stat_time = {}
for m in METRICS:
    d1 = div_df[div_df["Time"]=="1"][m].values
    d2 = div_df[div_df["Time"]=="2"][m].values
    _, p = mannwhitneyu(d1, d2, alternative="two-sided")
    stat_time[m] = p
    sig = "***" if p<0.001 else "**" if p<0.01 else "*" if p<0.05 else "ns"
    print(f"  {m:20s}: p={p:.4f} {sig}")

# Test 3: 4-group comparison (LGG_T1, LGG_T2, Placebo_T1, Placebo_T2)
print("\n[3] 4-group Kruskal-Wallis")
stat_4g = {}
for m in METRICS:
    groups_data = [div_df[div_df["Group_Time"]==g][m].values
                   for g in ["LGG_T1","LGG_T2","Placebo_T1","Placebo_T2"]]
    _, p = kruskal(*groups_data)
    stat_4g[m] = p
    sig = "***" if p<0.001 else "**" if p<0.01 else "*" if p<0.05 else "ns"
    print(f"  {m:20s}: p={p:.4f} {sig}")

# ─────────────────────────────────────────────────────────────
# SECTION 5: DASHBOARD
# ─────────────────────────────────────────────────────────────

print("\n🎨 Generating dashboard...")

PAL_GROUP    = {"LGG": "#E74C3C", "Placebo": "#3498DB"}
PAL_TIME     = {"1": "#F39C12",   "2": "#27AE60"}
PAL_4G       = {"LGG_T1": "#E74C3C", "LGG_T2": "#C0392B",
                "Placebo_T1": "#3498DB", "Placebo_T2": "#1A5276"}
MLABELS = {
    "Shannon"         : "Shannon Entropy (H')",
    "Simpson"         : "Simpson Index (1-D)",
    "Observed_OTUs"   : "Observed OTUs",
    "Chao1"           : "Chao1 Richness",
    "Pielou_Evenness" : "Pielou's Evenness",
}
ORDER_2G  = ["LGG", "Placebo"]
ORDER_4G  = ["LGG_T1", "LGG_T2", "Placebo_T1", "Placebo_T2"]

fig = plt.figure(figsize=(22, 18))
fig.suptitle(
    "16S rRNA Alpha Diversity — REAL DATA\n"
    "LGG Probiotic vs Placebo Intervention | peerj32 dataset\n"
    "44 samples × 130 gut taxa (Lahti et al.)",
    fontsize=15, fontweight="bold", y=0.99
)

# ── Row 1: Violin plots — LGG vs Placebo (5 metrics) ──
for i, m in enumerate(METRICS):
    ax = fig.add_subplot(4, 5, i+1)
    sns.violinplot(data=div_df, x="Group", y=m, palette=PAL_GROUP,
                   inner=None, alpha=0.45, order=ORDER_2G, ax=ax)
    sns.stripplot(data=div_df, x="Group", y=m, palette=PAL_GROUP,
                  size=5, jitter=True, order=ORDER_2G, alpha=0.8, ax=ax)
    p   = stat_group[m]
    sig = "***" if p<0.001 else "**" if p<0.01 else "*" if p<0.05 else "ns"
    ax.set_title(f"{MLABELS[m]}\nLGG vs Placebo: {sig}",
                 fontsize=8, fontweight="bold")
    ax.set_xlabel(""); ax.tick_params(labelsize=8)

# ── Row 2: 4-group boxplots (LGG/Placebo × Time1/Time2) ──
for i, m in enumerate(METRICS):
    ax = fig.add_subplot(4, 5, i+6)
    sns.boxplot(data=div_df, x="Group_Time", y=m,
                palette=PAL_4G, order=ORDER_4G,
                width=0.5, linewidth=1.2, ax=ax,
                flierprops={"marker":"o","markersize":3})
    p   = stat_4g[m]
    sig = "***" if p<0.001 else "**" if p<0.01 else "*" if p<0.05 else "ns"
    ax.set_title(f"{MLABELS[m]}\nKW: {sig}", fontsize=8, fontweight="bold")
    ax.set_xticklabels(["LGG\nT1","LGG\nT2","Plac\nT1","Plac\nT2"], fontsize=7)
    ax.set_xlabel("")

# ── Row 3 left: Shannon detailed — LGG vs Placebo ──
ax3 = fig.add_subplot(4, 2, 5)
sns.boxplot(data=div_df, x="Group", y="Shannon", palette=PAL_GROUP,
            order=ORDER_2G, width=0.4, linewidth=1.5, ax=ax3,
            flierprops={"marker":"o","markersize":4})
sns.stripplot(data=div_df, x="Group", y="Shannon", palette=PAL_GROUP,
              size=5, jitter=True, order=ORDER_2G, alpha=0.6, ax=ax3)
for i, g in enumerate(ORDER_2G):
    mv = div_df[div_df["Group"]==g]["Shannon"].mean()
    ax3.hlines(mv, i-0.2, i+0.2, colors="black", lw=2, linestyles="--")
    ax3.text(i, mv+0.01, f"μ={mv:.3f}", ha="center",
             fontsize=10, fontweight="bold")
ax3.set_title("Shannon Entropy — LGG vs Placebo", fontweight="bold", fontsize=11)
ax3.set_ylabel("Shannon Entropy (H')"); ax3.set_xlabel("")

# ── Row 3 right: Shannon over time ──
ax4 = fig.add_subplot(4, 2, 6)
sns.boxplot(data=div_df, x="Time", y="Shannon", hue="Group",
            palette=PAL_GROUP, width=0.5, ax=ax4,
            flierprops={"marker":"o","markersize":3})
ax4.set_title("Shannon Entropy: Before vs After\n(by Group)",
              fontweight="bold", fontsize=11)
ax4.set_xlabel("Time Point (1=Before, 2=After)")
ax4.set_ylabel("Shannon Entropy (H')")
ax4.legend(title="Group", fontsize=9)

# ── Row 4 left: Radar chart ──
ax5 = fig.add_subplot(4, 2, 7, polar=True)
cats   = ["Shannon","Simpson","Richness","Chao1","Evenness"]
N      = len(cats)
angles = [n/float(N)*2*np.pi for n in range(N)]
angles += angles[:1]
for g, color in PAL_GROUP.items():
    vals = div_df[div_df["Group"]==g][METRICS].mean().values
    mn   = div_df[METRICS].min().values
    mx   = div_df[METRICS].max().values
    norm = (vals - mn) / (mx - mn + 1e-9)
    norm = np.append(norm, norm[0])
    ax5.plot(angles, norm, "o-", lw=2, color=color, label=g)
    ax5.fill(angles, norm, alpha=0.1, color=color)
ax5.set_xticks(angles[:-1])
ax5.set_xticklabels(cats, fontsize=9)
ax5.set_title("Diversity Profile\n(Normalized)", fontweight="bold",
              fontsize=11, pad=20)
ax5.legend(loc="upper right", bbox_to_anchor=(1.35, 1.1), fontsize=9)

# ── Row 4 right: Summary stats table ──
ax6 = fig.add_subplot(4, 2, 8)
ax6.axis("off")
tdata, rlabels = [], []
for g in ORDER_2G:
    sub = div_df[div_df["Group"]==g]
    tdata.append([f"{sub[m].mean():.3f}±{sub[m].std():.3f}" for m in METRICS])
    rlabels.append(g)
# Add p-values row
tdata.append([
    f"p={stat_group[m]:.3f} {'*' if stat_group[m]<0.05 else 'ns'}"
    for m in METRICS])
rlabels.append("p-value")
tbl = ax6.table(cellText=tdata, rowLabels=rlabels,
                colLabels=["Shannon","Simpson","Obs.OTUs","Chao1","Evenness"],
                cellLoc="center", loc="center")
tbl.auto_set_font_size(False); tbl.set_fontsize(8); tbl.scale(1.2, 2.5)
for i, g in enumerate(ORDER_2G):
    tbl[(i+1,-1)].set_facecolor(PAL_GROUP[g])
    tbl[(i+1,-1)].set_text_props(color="white", fontweight="bold")
tbl[(3,-1)].set_facecolor("#BDC3C7")
ax6.set_title("Mean ± SD + Statistics", fontweight="bold", fontsize=11, pad=30)

plt.tight_layout(rect=[0,0,1,0.96])
plt.savefig("outputs/alpha_diversity_dashboard.png", dpi=150, bbox_inches="tight")
plt.close()
print("✅ Dashboard saved → outputs/alpha_diversity_dashboard.png")

# ─────────────────────────────────────────────────────────────
# FINAL SUMMARY
# ─────────────────────────────────────────────────────────────
print("\n" + "="*60)
print("FINAL SUMMARY — LGG vs Placebo")
print("="*60)
print(div_df.groupby("Group")[METRICS].mean().round(3).to_string())
print("\n" + "="*60)
print("FINAL SUMMARY — Time 1 vs Time 2")
print("="*60)
print(div_df.groupby("Time")[METRICS].mean().round(3).to_string())
print("\n✅ All outputs saved!")
