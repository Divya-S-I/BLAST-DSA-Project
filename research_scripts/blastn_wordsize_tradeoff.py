# =====================================================================
# TASK 1 - BLASTN word-size (K) trade-off: sensitivity vs execution time
# =====================================================================
# Setup (Ubuntu terminal, once):
#   sudo apt install ncbi-blast+
#   pip install biopython pandas numpy matplotlib seaborn
#
# DESIGN
#  1. One subject genome per row (real NCBI genome, or a random synthetic
#     genome for the two SYNTHETIC controls).
#  2. Queries = 120-bp fragments cut from the subject at KNOWN positions,
#     then given an exact number of substitutions -> exact identity.
#     (ground truth is therefore known, no guessing)
#  3. Same queries + same BLASTN settings for every K; only -word_size changes.
#  4. Sensitivity = fraction of queries with >=1 HSP overlapping the TRUE source
#     locus by >= MIN_OVERLAP bp.  Time = wall-clock of the blastn process
#     / number of queries, repeated N_REPEATS times (K order shuffled).
#  5. On-target alignment score = mean BLAST bitscore of the best HSP that
#     overlaps the true locus, averaged over the recovered queries.

import os, shutil, subprocess, time, random, math
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import Entrez, SeqIO

try:
    from IPython.display import display
except ImportError:                                   # plain terminal fallback
    def display(x):
        print(x.to_string())

# ---------------------------- settings -------------------------------
Entrez.email = "your_email@example.com"      # NCBI asks for an e-mail
WORKDIR = "blast_work"                       # temp FASTA/BLAST db files (no CSV)
K_VALUES = list(range(7, 16))                # word sizes 7..15
QUERY_LEN = 120                              # bp per query
N_QUERIES = 175                              # per dataset (divisible by 7 levels)
MIN_OVERLAP = 40                             # bp overlap with true locus = "recovered"
N_REPEATS = 3                                # timing repeats
SEED = 42
SYNTH_GENOME_LEN = 4_000_000
SHOW_LINE_PLOTS = False                      # False = only the block (heatmap) figure is shown; True = also show line graphs

# every BLASTN parameter except word_size is fixed here
BLAST_FIXED = dict(task="blastn", reward=2, penalty=-3, gapopen=5, gapextend=2,
                   evalue=1e-3, dust="no", soft_masking="false", num_threads=1)

# Identity levels of the planted homologs (controlled assumption, not measured):
INTRA_SPECIES = [0.98, 0.95, 0.92, 0.89, 0.86, 0.83, 0.80]   # same species / close relatives
INTER_GENUS   = [0.92, 0.89, 0.86, 0.83, 0.80, 0.77, 0.74]   # more distant relatives

# (row label exactly as in your heatmap, NCBI accession, kind, identity levels)
# The printed record description lets you verify each accession.
DATASETS = [
    ("K-12 self (control)",          "NC_000913.3", "real",      [1.00]),
    ("E. coli DH10B",                "NC_010473.1", "real",      INTRA_SPECIES),
    ("E. coli BL21(DE3)",            "NC_012892.2", "real",      INTRA_SPECIES),
    ("E. coli ATCC 8739",            "NC_010468.1", "real",      INTRA_SPECIES),
    ("E. coli O157:H7 Sakai",        "NC_002695.2", "real",      INTRA_SPECIES),
    ("E. coli CFT073",               "NC_004431.1", "real",      INTRA_SPECIES),
    ("Shigella flexneri 301",        "NC_004337.2", "real",      INTRA_SPECIES),
    ("Salmonella LT2",               "NC_003197.2", "real",      INTER_GENUS),
    ("Klebsiella HS11286",           "NC_016845.1", "real",      INTER_GENUS),
    ("SYNTHETIC 90% id (control)",   None,          "synthetic", [0.90]),
    ("SYNTHETIC 80% id (control)",   None,          "synthetic", [0.80]),
]
LABELS = [d[0] for d in DATASETS]

assert shutil.which("blastn") and shutil.which("makeblastdb"), "Install ncbi-blast+ first"
os.makedirs(WORKDIR, exist_ok=True)

# ---------------------------- helpers --------------------------------
def fetch_genome(acc):
    """Download a genome FASTA once (cached in WORKDIR), return (sequence, description)."""
    path = f"{WORKDIR}/{acc}.fasta"
    if not os.path.exists(path):
        h = Entrez.efetch(db="nuccore", id=acc, rettype="fasta", retmode="text")
        open(path, "w").write(h.read()); h.close(); time.sleep(0.4)
    rec = SeqIO.read(path, "fasta")
    return str(rec.seq).upper(), rec.description

def random_genome(length, rng):
    """Synthetic subject: i.i.d. bases, 50% GC."""
    return "".join(rng.choices("ACGT", k=length))

def mutate(seq, identity, rng):
    """Substitute exactly round((1-identity)*len) positions -> exact identity."""
    s = list(seq)
    for pos in rng.sample(range(len(s)), round((1 - identity) * len(s))):
        s[pos] = rng.choice([b for b in "ACGT" if b != s[pos]])
    return "".join(s)

def make_queries(genome, levels, rng):
    """Cut fragments at known positions and mutate them. Returns {qid: (seq, start0)}."""
    per_level = N_QUERIES // len(levels)
    queries = {}
    for lvl in levels:
        for i in range(per_level):
            while True:                                   # avoid N / ambiguous bases
                start = rng.randrange(0, len(genome) - QUERY_LEN)
                frag = genome[start:start + QUERY_LEN]
                if set(frag) <= set("ACGT"):
                    break
            queries[f"q{len(queries)}_id{round(lvl*100)}"] = (mutate(frag, lvl, rng), start)
    return queries

def write_fasta(path, seqs):
    with open(path, "w") as f:
        for name, s in seqs.items():
            f.write(f">{name}\n{s}\n")

def make_db(fasta, prefix):
    if not os.path.exists(prefix + ".nsq"):
        subprocess.run(["makeblastdb", "-in", fasta, "-dbtype", "nucl", "-out", prefix],
                       check=True, capture_output=True)

def run_blastn(query_fa, db, k):
    """One real BLASTN run. Returns (seconds, tabular output)."""
    cmd = ["blastn", "-query", query_fa, "-db", db, "-word_size", str(k),
           "-outfmt", "6 qseqid sstart send bitscore"]
    for key, val in BLAST_FIXED.items():
        cmd += [f"-{key}", str(val)]
    t0 = time.perf_counter()
    out = subprocess.run(cmd, capture_output=True, text=True, check=True).stdout
    return time.perf_counter() - t0, out

def recovered_queries(out, queries):
    """Ground truth check: HSP must overlap the planted source locus.
    Returns {qid: best on-target bitscore}."""
    found = {}
    for line in out.strip().splitlines():
        q, s0, s1, bits = line.split("\t")
        lo, hi = sorted((int(s0), int(s1)))
        t0 = queries[q][1] + 1                      # true locus, 1-based
        t1 = t0 + QUERY_LEN - 1
        if min(hi, t1) - max(lo, t0) + 1 >= MIN_OVERLAP:
            found[q] = max(found.get(q, 0.0), float(bits))
    return found

def wilson(h, n, z=1.96):
    """95% Wilson confidence interval for a proportion h/n."""
    p = h / n; d = 1 + z*z/n
    c = (p + z*z/(2*n)) / d
    w = z * math.sqrt(p*(1-p)/n + z*z/(4*n*n)) / d
    return c - w, c + w

# ------------------------- run the experiment ------------------------
rows = []
for i, (label, acc, kind, levels) in enumerate(DATASETS):
    rng = random.Random(SEED + i)
    if acc:
        genome, desc = fetch_genome(acc)
        print(f"[REAL ] {label:28s} <- {acc}: {desc[:70]} ({len(genome):,} bp)")
    else:
        genome = random_genome(SYNTH_GENOME_LEN, rng)
        print(f"[SYNTH] {label:28s} <- random genome ({len(genome):,} bp)")

    queries = make_queries(genome, levels, rng)
    qfa, sfa, db = f"{WORKDIR}/q_{i}.fa", f"{WORKDIR}/s_{i}.fa", f"{WORKDIR}/db_{i}"
    write_fasta(qfa, {q: v[0] for q, v in queries.items()})
    write_fasta(sfa, {"subject": genome})
    make_db(sfa, db)
    run_blastn(qfa, db, 11)                          # warm-up (file cache), not recorded

    for rep in range(N_REPEATS):
        ks = K_VALUES[:]; rng.shuffle(ks)            # shuffle K order to avoid drift bias
        for k in ks:
            secs, out = run_blastn(qfa, db, k)
            rec = recovered_queries(out, queries)
            rows.append(dict(order=i, dataset=label, kind=kind, K=k, repeat=rep,
                             time_s=secs, n_queries=len(queries),
                             recovered=len(rec),
                             mean_score=np.mean(list(rec.values())) if rec else np.nan))
    print("        done.")

df = pd.DataFrame(rows)
df["sensitivity"] = df.recovered / df.n_queries
df["ms_per_query"] = df.time_s / df.n_queries * 1000
df["s_per_query"] = df.time_s / df.n_queries

# ------------------------- summaries ---------------------------------
per = (df.groupby(["order", "dataset", "kind", "K"])
         .agg(sensitivity=("sensitivity", "mean"), ms=("ms_per_query", "mean"),
              ms_sd=("ms_per_query", "std"), hits=("recovered", "mean"),
              sec=("s_per_query", "mean"), score=("mean_score", "mean"),
              n=("n_queries", "first")).reset_index())

# pooled curves: real genomes vs synthetic controls (kept separate on purpose)
GROUPS = {"Real genomes (9 rows)": "real", "Synthetic controls (2 rows)": "synthetic"}
COLORS = {"Real genomes (9 rows)": "#1f77b4", "Synthetic controls (2 rows)": "#d62728"}
pooled = []
for g, kind in GROUPS.items():
    sub = df[df.kind == kind]
    for k in K_VALUES:
        s = sub[sub.K == k]
        rep_means = s.groupby("repeat").ms_per_query.mean()    # group mean per repeat
        s0 = s[s.repeat == 0]                                  # BLAST is deterministic
        h, n = s0.recovered.sum(), s0.n_queries.sum()
        lo, hi = wilson(h, n)
        pooled.append(dict(group=g, K=k, sensitivity=h/n, ci_lo=lo, ci_hi=hi,
                           ms_per_query=rep_means.mean(), ms_sd=rep_means.std(ddof=1)))
pooled = pd.DataFrame(pooled)

print("\nPOOLED SUMMARY (time error = SD over repeats; sensitivity error = 95% Wilson CI)")
display(pooled.round(4))

real = pooled[pooled.group == "Real genomes (9 rows)"].set_index("K")
print(f"Real genomes: K=15 is {real.ms_per_query[7]/real.ms_per_query[15]:.1f}x faster than K=7, "
      f"sensitivity {real.sensitivity[7]:.3f} -> {real.sensitivity[15]:.3f}")

# ------------------------- FIGURE 1: trade-off -----------------------
fig, ax = plt.subplots(1, 3, figsize=(18, 5))
for g in GROUPS:
    d, c = pooled[pooled.group == g], COLORS[g]
    yerr = [d.sensitivity - d.ci_lo, d.ci_hi - d.sensitivity]
    ax[0].errorbar(d.K, d.ms_per_query, yerr=d.ms_sd, marker="o", capsize=3, color=c, label=g)
    ax[1].errorbar(d.K, d.sensitivity, yerr=yerr, marker="o", capsize=3, color=c, label=g)
    ax[2].errorbar(d.ms_per_query, d.sensitivity, xerr=d.ms_sd, yerr=yerr,
                   marker="o", capsize=3, color=c, label=g)
    for _, r in d.iterrows():
        ax[2].annotate(f"K={int(r.K)}", (r.ms_per_query, r.sensitivity),
                       textcoords="offset points", xytext=(5, 6), fontsize=8, color=c)
for dsn in per.dataset.unique():                     # faint per-dataset lines for context
    d = per[per.dataset == dsn]
    c = COLORS[[g for g, k in GROUPS.items() if k == d.kind.iloc[0]][0]]
    ax[0].plot(d.K, d.ms, color=c, alpha=0.15)
    ax[1].plot(d.K, d.sensitivity, color=c, alpha=0.15)
ax[0].set(xlabel="Word size K", ylabel="Time per query (ms)", title="(a) Execution time vs K")
ax[1].set(xlabel="Word size K", ylabel="Sensitivity", title="(b) Sensitivity vs K")
ax[2].set(xlabel="Time per query (ms)", ylabel="Sensitivity", title="(c) Trade-off: sensitivity vs time")
for a in ax:
    a.grid(alpha=0.3); a.legend(fontsize=8)
plt.tight_layout()
fig.savefig("figure1_tradeoff.png", dpi=150)
if not SHOW_LINE_PLOTS:
    plt.close(fig)                                   # skip the line graphs

# ------------------------- FIGURE 2: dataset heatmaps ----------------
def heat(a, mat, title, cmap, fmt, vmax=None):
    sns.heatmap(mat, ax=a, cmap=cmap, annot=True, fmt=fmt, vmax=vmax,
                annot_kws={"size": 6.5}, linewidths=0.6, linecolor="white",
                cbar_kws={"shrink": 0.65})
    a.set_title(title, fontweight="bold", fontsize=12)
    a.set_xlabel("word size K")
    a.set_ylabel("")
    a.set_yticklabels(LABELS, rotation=0, fontsize=10)
    a.set_xticklabels(K_VALUES, rotation=0)

sens_mat  = per.pivot(index="dataset", columns="K", values="sensitivity").loc[LABELS]
time_mat  = per.pivot(index="dataset", columns="K", values="ms").loc[LABELS]
score_mat = per.pivot(index="dataset", columns="K", values="score").loc[LABELS]

print("\nMEAN TIME PER QUERY (ms) - the numbers plotted in the middle heatmap:")
print(time_mat.round(2).to_string())

fig, ax = plt.subplots(1, 3, figsize=(22, 8))
heat(ax[0], sens_mat,  "Sensitivity",                    "YlGnBu",   ".2f", vmax=1.0)
heat(ax[1], time_mat,  "Mean time per query (ms)",      "rocket_r", ".1f")
heat(ax[2], score_mat, "Mean on-target alignment score", "mako_r",   ".0f")
fig.suptitle("Dataset x word-size heatmaps  [v3: time in ms]", fontsize=15, fontweight="bold")
plt.tight_layout()
fig.savefig("figure2_heatmaps.png", dpi=150)
plt.show()                                           # opens the heatmap window(s)
