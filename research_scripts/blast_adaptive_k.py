# =====================================================================
# TASK 2 - Adaptive word-size (K) selection for BLASTN  (wide K range)
# =====================================================================
# Setup (Ubuntu terminal, once):
#   sudo apt install ncbi-blast+
#   pip install biopython pandas numpy matplotlib seaborn
#
# SOURCE (established part)
#   Ma B., Tromp J., Li M. (2002) "PatternHunter: faster and more sensitive
#   homology search", Bioinformatics 18(3):440-445.
#   Keich U., Li M., Ma B., Tromp J. (2004) "On spaced seeds for similarity
#   search", Discrete Applied Mathematics 138(3):253-263.
#
#   Model: a homologous region is a 0/1 string of length N (1 = match), each
#   position matching independently with probability p. A BLASTN hit with
#   word size K is a run of K consecutive 1s. Sensitivity of K:
#
#       S(K, p, N) = Pr[ at least one run of K consecutive matches ]
#
#   Computed exactly by dynamic programming over f_i(j) = Pr[no hit yet after
#   i positions AND the current run of matches has length j], 0 <= j < K:
#       f_{i+1}(0)   = (1 - p) * sum_j f_i(j)       (mismatch resets the run)
#       f_{i+1}(j+1) = p * f_i(j)    for j < K-1    (match extends the run)
#       S = 1 - sum_j f_N(j)                         (K-run reached -> hit)
#
#   Cost side (PatternHunter): random seed hits shrink ~4x for every +1 in K:
#       E_rand = 2 * (m-K+1) * (n-K+1) / 4^K
#
# WHAT IS MY ADAPTATION (not a published BLASTN formula)
#   * p is estimated by a cheap pilot BLASTN run (K = K_MIN) on a few queries.
#   * Rule: K* = the LARGEST K in [K_MIN, K_MAX] whose predicted
#     S(K, p_hat, N) >= TARGET.  Larger K -> fewer random hits -> faster.
#   * Pilot time is counted in the adaptive method's total time.
#
# WHAT IS NEW IN THIS VERSION
#   * K is searched over a WIDE range (K_MIN..K_MAX), not just 7..15.
#   * 5 fixed K values + random K + adaptive K are compared.
#   * 5 identity scenarios -> adaptive K picks 5 different, sensible K values.
#   * New figures: sensitivity "landscape" S(K, p) with the target contour,
#     method x scenario heatmaps (red box = missed target), speed-up heatmap,
#     log-time trade-off per scenario, and a predicted-vs-measured model check.

import os, shutil, subprocess, time, random, itertools
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Rectangle
from Bio import Entrez, SeqIO

try:
    from IPython.display import display
except ImportError:                                   # plain terminal fallback
    def display(x):
        print(x.to_string())

# ---------------------------- settings -------------------------------
Entrez.email = "your_email@example.com"
WORKDIR = "blast_work"
SUBJECT_ACC = "NC_000913.3"                  # real E. coli K-12 MG1655 genome

K_MIN, K_MAX = 8, 28                         # searchable word-size range (change freely)
K_VALUES = list(range(K_MIN, K_MAX + 1))
FIXED_KS = [8, 11, 14, 18, 24]               # the 5 fixed-K baselines (must lie in K range)
PILOT_K = K_MIN                              # most sensitive K, used only for the pilot

QLEN = 100                                   # query length = region length N
N_Q = 40                                     # queries per scenario
PILOT_N = 5                                  # queries used for the pilot estimate of p
TARGET = 0.90                                # required predicted sensitivity
N_REPEATS = 3
N_RANDOM = 500                               # random-K draws
MIN_OVERLAP = 40
SEED = 7

# 5 scenarios: real K-12 fragments mutated to known identity levels
SCENARIOS = {"98% identity": 0.98, "93% identity": 0.93, "88% identity": 0.88,
             "83% identity": 0.83, "78% identity": 0.78}
BLAST_FIXED = dict(task="blastn", reward=2, penalty=-3, gapopen=5, gapextend=2,
                   evalue=1e-3, dust="no", soft_masking="false", num_threads=1)

assert shutil.which("blastn") and shutil.which("makeblastdb"), "Install ncbi-blast+ first"
assert all(k in K_VALUES for k in FIXED_KS), "FIXED_KS must lie inside K_MIN..K_MAX"
os.makedirs(WORKDIR, exist_ok=True)

# ------------------------ the adaptive method ------------------------
def hit_probability(K, p, N):
    """S(K,p,N): probability of >=1 run of K matches in N positions (DP above)."""
    f = [0.0] * K; f[0] = 1.0                 # f[j] = Pr(no hit yet, current run = j)
    for _ in range(N):
        new = [0.0] * K
        for j, pr in enumerate(f):
            new[0] += pr * (1 - p)            # mismatch -> run resets to 0
            if j + 1 < K:
                new[j + 1] += pr * p          # match -> run grows
            # else: run reaches K -> a hit, removed from the "no hit" mass
        f = new
    return 1.0 - sum(f)

def expected_random_hits(K, m, n):
    """Expected random word hits for one query (both strands)."""
    return 2 * (m - K + 1) * (n - K + 1) / 4 ** K

def select_K(p_hat, N, target=TARGET):
    """Largest K in [K_MIN, K_MAX] that still meets the predicted-sensitivity target."""
    S = {k: hit_probability(k, p_hat, N) for k in K_VALUES}
    ok = [k for k in K_VALUES if S[k] >= target]
    return (max(ok) if ok else K_MIN), S

def sanity_check():
    """Brute-force check of the DP on a tiny case, so you can verify it by hand."""
    K, p, N = 3, 0.8, 6
    total = 0.0
    for bits in itertools.product([0, 1], repeat=N):
        prob = np.prod([p if b else 1 - p for b in bits])
        if "1" * K in "".join(map(str, bits)):
            total += prob
    print(f"sanity check K={K}, p={p}, N={N}: enumeration={total:.6f}  DP={hit_probability(K, p, N):.6f}")

# ---------------------------- helpers --------------------------------
def fetch_genome(acc):
    path = f"{WORKDIR}/{acc}.fasta"
    if not os.path.exists(path):
        h = Entrez.efetch(db="nuccore", id=acc, rettype="fasta", retmode="text")
        open(path, "w").write(h.read()); h.close(); time.sleep(0.4)
    rec = SeqIO.read(path, "fasta")
    return str(rec.seq).upper(), rec.description

def mutate(seq, identity, rng):
    s = list(seq)
    for pos in rng.sample(range(len(s)), round((1 - identity) * len(s))):
        s[pos] = rng.choice([b for b in "ACGT" if b != s[pos]])
    return "".join(s)

def make_queries(genome, identity, rng):
    q = {}
    for i in range(N_Q):
        while True:
            start = rng.randrange(0, len(genome) - QLEN)
            frag = genome[start:start + QLEN]
            if set(frag) <= set("ACGT"):
                break
        q[f"q{i}"] = (mutate(frag, identity, rng), start)
    return q

def write_fasta(path, seqs):
    with open(path, "w") as f:
        for n, s in seqs.items():
            f.write(f">{n}\n{s}\n")

def run_blastn(query_fa, db, k):
    cmd = ["blastn", "-query", query_fa, "-db", db, "-word_size", str(k),
           "-outfmt", "6 qseqid sstart send pident length bitscore"]
    for key, val in BLAST_FIXED.items():
        cmd += [f"-{key}", str(val)]
    t0 = time.perf_counter()
    out = subprocess.run(cmd, capture_output=True, text=True, check=True).stdout
    return time.perf_counter() - t0, out

def recovered(out, queries):
    found = set()
    for line in out.strip().splitlines():
        q, s0, s1, *_ = line.split("\t")
        lo, hi = sorted((int(s0), int(s1)))
        t0 = queries[q][1] + 1; t1 = t0 + QLEN - 1
        if min(hi, t1) - max(lo, t0) + 1 >= MIN_OVERLAP:
            found.add(q)
    return found

def estimate_p(out):
    """Pilot estimate: length-weighted identity over each query's best HSP."""
    best = {}
    for line in out.strip().splitlines():
        q, _, _, pid, ln, bits = line.split("\t")
        if q not in best or float(bits) > best[q][2]:
            best[q] = (float(pid), int(ln), float(bits))
    if not best:
        return None
    return sum(p * l for p, l, _ in best.values()) / sum(l for _, l, _ in best.values()) / 100

# ------------------------- data + database ---------------------------
sanity_check()
genome, desc = fetch_genome(SUBJECT_ACC)
print(f"[REAL] subject: {desc[:70]} ({len(genome):,} bp)")
sfa, db = f"{WORKDIR}/t2_subject.fa", f"{WORKDIR}/t2_db"
write_fasta(sfa, {"subject": genome})
if not os.path.exists(db + ".nsq"):
    subprocess.run(["makeblastdb", "-in", sfa, "-dbtype", "nucl", "-out", db],
                   check=True, capture_output=True)

# ------------------------- run everything ----------------------------
fixed_rows, adaptive_rows, predicted = [], [], {}
for si, (sc, identity) in enumerate(SCENARIOS.items()):
    rng = random.Random(SEED + si)
    queries = make_queries(genome, identity, rng)
    qfa = f"{WORKDIR}/t2_q{si}.fa"; write_fasta(qfa, {q: v[0] for q, v in queries.items()})
    pilot_ids = rng.sample(list(queries), PILOT_N)
    pfa = f"{WORKDIR}/t2_pilot{si}.fa"; write_fasta(pfa, {q: queries[q][0] for q in pilot_ids})
    run_blastn(qfa, db, 11)                                    # warm-up (not recorded)
    print(f"\n[{sc}] running fixed K = {K_MIN}..{K_MAX} and adaptive ...")

    # (a) fixed K = every K in range (also the pool that "random K" draws from)
    for rep in range(N_REPEATS):
        ks = K_VALUES[:]; rng.shuffle(ks)
        for k in ks:
            secs, out = run_blastn(qfa, db, k)
            fixed_rows.append(dict(scenario=sc, K=k, rep=rep, ms=secs / N_Q * 1000,
                                   sens=len(recovered(out, queries)) / N_Q))

    # (b) adaptive: pilot -> p_hat -> K* -> final run (pilot time included)
    for rep in range(N_REPEATS):
        t0 = time.perf_counter()
        _, pout = run_blastn(pfa, db, PILOT_K)
        p_hat = estimate_p(pout)
        if p_hat is None:                                      # pilot found nothing -> be cautious
            print("   warning: pilot found no hits, using p_hat = 0.70")
            p_hat = 0.70
        p_hat = min(p_hat, 0.999)
        K_sel, S = select_K(p_hat, QLEN)
        _, out = run_blastn(qfa, db, K_sel)
        total = time.perf_counter() - t0
        adaptive_rows.append(dict(scenario=sc, rep=rep, p_hat=p_hat, K_sel=K_sel,
                                  ms=total / N_Q * 1000, sens=len(recovered(out, queries)) / N_Q))
    predicted[sc] = (p_hat, K_sel, S)

    tab = pd.DataFrame({"K": K_VALUES,
                        "predicted S(K)": [S[k] for k in K_VALUES],
                        "E_rand per query": [expected_random_hits(k, QLEN, len(genome)) for k in K_VALUES]})
    tab["meets target"] = tab["predicted S(K)"] >= TARGET
    print(f"=== {sc}: true identity {identity:.2f}, pilot p_hat = {p_hat:.3f} -> selected K = {K_sel} ===")
    display(tab.round(4))

# ------------------------- summaries ---------------------------------
fixed = (pd.DataFrame(fixed_rows).groupby(["scenario", "K"])
           .agg(sens=("sens", "mean"), ms=("ms", "mean"), ms_sd=("ms", "std"))
           .reset_index().fillna(0))
ad = pd.DataFrame(adaptive_rows).groupby("scenario").agg(
    K_sel=("K_sel", "first"), p_hat=("p_hat", "first"),
    sens=("sens", "mean"), ms=("ms", "mean"), ms_sd=("ms", "std")).fillna(0)

rng = random.Random(SEED)
draws = [rng.choice(K_VALUES) for _ in range(N_RANDOM)]       # random K from the whole range
METHODS = [f"Fixed K={k}" for k in FIXED_KS] + ["Random K", "Adaptive K"]
res, oracle, base_ms = {}, {}, {}
for sc in SCENARIOS:
    f = fixed[fixed.scenario == sc].set_index("K")
    res[sc] = {}
    for k in FIXED_KS:                                         # (sens, sens_sd, ms, ms_sd)
        res[sc][f"Fixed K={k}"] = (f.loc[k, "sens"], 0.0, f.loc[k, "ms"], f.loc[k, "ms_sd"])
    rs, rt = f.loc[draws, "sens"], f.loc[draws, "ms"]
    res[sc]["Random K"] = (rs.mean(), rs.std(), rt.mean(), rt.std())
    res[sc]["Adaptive K"] = (ad.loc[sc, "sens"], 0.0, ad.loc[sc, "ms"], ad.loc[sc, "ms_sd"])
    good = [k for k in K_VALUES if f.loc[k, "sens"] >= TARGET]
    oracle[sc] = max(good) if good else K_MIN                  # best K seen in hindsight
    base_ms[sc] = f.loc[K_MIN, "ms"]                           # baseline = most sensitive K
rand_k_mean, rand_k_sd = np.mean(draws), np.std(draws)

# scoreboard
rows = []
for m in METHODS:
    sens = [res[sc][m][0] for sc in SCENARIOS]
    ms = [res[sc][m][2] for sc in SCENARIOS]
    sp = [base_ms[sc] / res[sc][m][2] for sc in SCENARIOS]
    rows.append(dict(method=m, mean_sens=np.mean(sens), min_sens=np.min(sens),
                     scenarios_met_target=f"{sum(s >= TARGET for s in sens)}/{len(SCENARIOS)}",
                     _met=sum(s >= TARGET for s in sens),
                     mean_ms=np.mean(ms), mean_speedup_vs_K_MIN=np.mean(sp)))
board = (pd.DataFrame(rows).sort_values(["_met", "mean_speedup_vs_K_MIN"], ascending=False)
           .drop(columns="_met").reset_index(drop=True))
print(f"\nSCOREBOARD (target sensitivity = {TARGET}; speed-up relative to fixed K={K_MIN})")
display(board.round(3))
print("\nSelected K per scenario (adaptive vs hindsight-best fixed K):")
for sc in SCENARIOS:
    print(f"  {sc:14s} adaptive K = {predicted[sc][1]:2d} | best fixed K in hindsight = {oracle[sc]:2d}")

# ------------------------- colours / saving --------------------------
blues = plt.cm.Blues(np.linspace(0.35, 0.9, len(FIXED_KS)))
MCOL = {m: blues[i] for i, m in enumerate(METHODS[:len(FIXED_KS)])}
MCOL["Random K"], MCOL["Adaptive K"] = "#fd8d3c", "#31a354"
SC_COLS = dict(zip(SCENARIOS, ["#1f77b4", "#2ca02c", "#ff7f0e", "#d62728", "#9467bd"]))

def save(fig, name):
    fig.savefig(name, dpi=150, bbox_inches="tight")

# ------------------------- FIGURE 1: how K is selected ---------------
fig1, ax = plt.subplots(1, 2, figsize=(15, 5.2))
for sc, (p_hat, K_sel, S) in predicted.items():
    c = SC_COLS[sc]
    ax[0].plot(K_VALUES, [S[k] for k in K_VALUES], "-o", ms=4, color=c, label=f"{sc}  (p_hat={p_hat:.2f})")
    ax[0].plot(K_sel, S[K_sel], "*", ms=19, mec="k", color=c, zorder=5)
ax[0].axhline(TARGET, ls="--", c="gray", label=f"target S = {TARGET}")
ax[0].set(xlabel="Word size K", ylabel="Predicted sensitivity S(K, p_hat, N)",
          title="(a) Selection rule (stars = selected K)")
ax[0].legend(fontsize=8); ax[0].grid(alpha=0.3)

x, w = np.arange(len(SCENARIOS)), 0.26
ax[1].bar(x - w, [predicted[s][1] for s in SCENARIOS], w, color=MCOL["Adaptive K"], label="Adaptive K")
ax[1].bar(x,     [oracle[s] for s in SCENARIOS], w, color="#756bb1", label="Best fixed K in hindsight")
ax[1].bar(x + w, [rand_k_mean] * len(SCENARIOS), w, yerr=rand_k_sd, capsize=3,
          color=MCOL["Random K"], label="Random K (mean ± SD)")
ax[1].set_xticks(x); ax[1].set_xticklabels(SCENARIOS, fontsize=8)
ax[1].set(ylabel="Word size K", title="(b) Selected K per scenario"); ax[1].legend(fontsize=8)
ax[1].grid(axis="y", alpha=0.3)
fig1.tight_layout(); save(fig1, "t2_fig1_selection.png")

# ------------------------- FIGURE 2: sensitivity landscape -----------
P_GRID = np.round(np.linspace(0.60, 1.00, 41), 3)
S_GRID = np.array([[hit_probability(k, p, QLEN) for k in K_VALUES] for p in P_GRID])
fig2, a = plt.subplots(figsize=(12, 6.2))
mesh = a.pcolormesh(K_VALUES, P_GRID, S_GRID, shading="nearest", cmap="YlGnBu", vmin=0, vmax=1)
cs = a.contour(K_VALUES, P_GRID, S_GRID, levels=[TARGET], colors="red", linewidths=2.8)
a.clabel(cs, fmt={TARGET: f"S = {TARGET}"}, fontsize=10)
for sc, (p_hat, K_sel, _) in predicted.items():
    a.plot(K_sel, p_hat, "*", ms=20, mec="k", color=SC_COLS[sc], zorder=6)
    a.annotate(f"{sc}\nK*={K_sel}", (K_sel, p_hat), textcoords="offset points", xytext=(8, -4),
               fontsize=8, fontweight="bold", color="white",
               bbox=dict(boxstyle="round,pad=0.2", fc="black", alpha=0.45, lw=0))
a.set(xlabel="Word size K", ylabel="Per-base identity p",
      title=f"Sensitivity landscape S(K, p, N={QLEN}) - red line = target; stars = adaptive choices")
fig2.colorbar(mesh, ax=a, label="Predicted sensitivity")
fig2.tight_layout(); save(fig2, "t2_fig2_landscape.png")

# ------------------------- FIGURE 3: method x scenario heatmaps ------
sens_mat = pd.DataFrame({sc: {m: res[sc][m][0] for m in METHODS} for sc in SCENARIOS}).loc[METHODS]
time_mat = pd.DataFrame({sc: {m: res[sc][m][2] for m in METHODS} for sc in SCENARIOS}).loc[METHODS]
speed_mat = pd.DataFrame({sc: {m: base_ms[sc] / res[sc][m][2] for m in METHODS} for sc in SCENARIOS}).loc[METHODS]

sens_annot = sens_mat.astype(object).copy()
for m in METHODS:
    for sc in SCENARIOS:
        txt = f"{sens_mat.loc[m, sc]:.2f}"
        if m == "Adaptive K":
            txt += f"\nK={predicted[sc][1]}"
        sens_annot.loc[m, sc] = txt
speed_annot = speed_mat.astype(object).copy()
for m in METHODS:
    for sc in SCENARIOS:
        speed_annot.loc[m, sc] = f"{speed_mat.loc[m, sc]:.1f}x"

def heat(a, mat, title, cmap, annot, vmin=None, vmax=None, fmt=""):
    sns.heatmap(mat, ax=a, cmap=cmap, annot=annot, fmt=fmt, vmin=vmin, vmax=vmax,
                annot_kws={"size": 10}, linewidths=1.0, linecolor="white",
                cbar_kws={"shrink": 0.7})
    a.set_title(title, fontweight="bold", fontsize=12)
    a.set_xlabel(""); a.set_ylabel("")
    a.set_yticklabels(mat.index, rotation=0, fontsize=10)
    a.set_xticklabels(mat.columns, rotation=20, ha="right", fontsize=9)

fig3, ax = plt.subplots(1, 3, figsize=(24, 7))
heat(ax[0], sens_mat, f"Sensitivity (red box = below target {TARGET})", "YlGnBu", sens_annot, vmin=0.4, vmax=1.0)
for r, m in enumerate(METHODS):
    for c, sc in enumerate(SCENARIOS):
        if sens_mat.loc[m, sc] < TARGET:
            ax[0].add_patch(Rectangle((c, r), 1, 1, fill=False, ec="red", lw=3))
heat(ax[1], time_mat, "Mean time per query (ms)", "rocket_r", True, fmt=".2f")
heat(ax[2], speed_mat, f"Speed-up vs fixed K={K_MIN}", "viridis", speed_annot)
fig3.suptitle("Method x scenario comparison", fontsize=15, fontweight="bold")
fig3.tight_layout(); save(fig3, "t2_fig3_method_heatmaps.png")

# ------------------------- FIGURE 4: trade-off + model check ---------
fig4, axes = plt.subplots(2, 3, figsize=(21, 11)); axes = axes.ravel()
for a, sc in zip(axes[:5], SCENARIOS):
    f = fixed[fixed.scenario == sc].sort_values("K")
    a.axhspan(TARGET, 1.05, color="#31a354", alpha=0.10)
    a.plot(f.ms, f.sens, "-", color="#9ecae1", lw=1.5, zorder=1)
    sc_pts = a.scatter(f.ms, f.sens, c=f.K, cmap="viridis", s=45, zorder=3, edgecolor="k", lw=0.4)
    for _, r in f.iterrows():
        if int(r.K) in FIXED_KS:
            a.annotate(f"K={int(r.K)}", (r.ms, r.sens), textcoords="offset points", xytext=(5, 6), fontsize=8)
    rs, rss, rt, rts = res[sc]["Random K"]
    a.errorbar(rt, rs, xerr=rts, yerr=rss, fmt="X", ms=12, color=MCOL["Random K"], mec="k",
               capsize=3, zorder=4, label="Random K (mean ± SD)")
    as_, _, at, ats = res[sc]["Adaptive K"]
    a.errorbar(at, as_, xerr=ats, fmt="*", ms=22, color=MCOL["Adaptive K"], mec="k",
               capsize=3, zorder=5, label=f"Adaptive (K={predicted[sc][1]})")
    a.axhline(TARGET, ls="--", c="gray")
    a.set_xscale("log")
    a.set(xlabel="Time per query (ms, log scale)", ylabel="Sensitivity", title=sc)
    a.set_ylim(max(0, f.sens.min() - 0.05), 1.05)
    a.grid(alpha=0.3, which="both"); a.legend(fontsize=8, loc="lower right")
    fig4.colorbar(sc_pts, ax=a, label="K", pad=0.01)

# 6th panel: model check, predicted (from known identity) vs measured sensitivity
a = axes[5]
for sc, identity in SCENARIOS.items():
    f = fixed[fixed.scenario == sc].sort_values("K")
    pred = [hit_probability(k, identity, QLEN) for k in f.K]
    a.scatter(pred, f.sens, s=36, color=SC_COLS[sc], alpha=0.85, edgecolor="k", lw=0.3, label=sc)
a.plot([0, 1], [0, 1], "k--", lw=1)
a.set(xlabel="Predicted S(K, p, N) from the DP model", ylabel="Measured sensitivity (BLASTN)",
      title="Model check: predicted vs measured (all K, all scenarios)")
a.set_xlim(-0.02, 1.02); a.set_ylim(-0.02, 1.02); a.grid(alpha=0.3); a.legend(fontsize=8)
fig4.suptitle("Sensitivity-time trade-off per scenario (green band = meets target)",
              fontsize=15, fontweight="bold")
fig4.tight_layout(); save(fig4, "t2_fig4_tradeoff.png")

print("\nSaved: t2_fig1_selection.png, t2_fig2_landscape.png, "
      "t2_fig3_method_heatmaps.png, t2_fig4_tradeoff.png")
plt.show()
