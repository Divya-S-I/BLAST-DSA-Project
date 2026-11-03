"""Reproducible, bounded five-step adaptive word-size experiment using local BLAST+.

The test queries are planted, substitution-mutated fragments from a reference.
Ground truth source positions are known. This is a controlled sensitivity test,
not a proof of biological homology or general BLAST sensitivity.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import random
import shutil
import subprocess
import tempfile
import time
from typing import Callable


BLAST_SETTINGS = {
    "task": "blastn", "reward": 2, "penalty": -3,
    "gapopen": 5, "gapextend": 2, "evalue": 1e-3,
    "dust": "no", "soft_masking": "false", "num_threads": 1,
}


def available() -> bool:
    return shutil.which("blastn") is not None and shutil.which("makeblastdb") is not None


def fasta_to_reference(text: str) -> str:
    """Concatenate the sequence lines of a single-FASTA-record input.

    Reject multiple records so known-locus offsets remain meaningful.
    """
    seqs = []
    current = []
    seen_header = False
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if seen_header:
                seqs.append("".join(current))
                current = []
            seen_header = True
        else:
            if not seen_header and not current:
                seen_header = True
            current.append(line.upper())
    seqs.append("".join(current))
    if len(seqs) != 1:
        raise ValueError("Provide one reference FASTA record. Multiple contigs are not supported in this experiment.")
    seq = seqs[0]
    if not seq or any(ch not in "ACGTN" for ch in seq):
        raise ValueError("Reference must contain only A/C/G/T/N bases (one FASTA record).")
    return seq


def hit_probability(k: int, p: float, N: int) -> float:
    """Exact probability of >= 1 length-k all-match run in N iid Bernoulli(p) trials."""
    if k > N:
        return 0.0
    states = [0.0] * k
    states[0] = 1.0
    for _ in range(N):
        nxt = [0.0] * k
        for j, mass in enumerate(states):
            nxt[0] += mass * (1 - p)
            if j + 1 < k:
                nxt[j + 1] += mass * p
        states = nxt
    return max(0.0, min(1.0, 1 - sum(states)))


def choose_k(p_hat: float, length: int, lower: int, upper: int, target: float):
    probabilities = {k: hit_probability(k, p_hat, length) for k in range(lower, upper + 1)}
    eligible = [k for k, s in probabilities.items() if s >= target]
    return (max(eligible) if eligible else lower), bool(eligible), probabilities


def mutate_exact(seq: str, identity: float, rng: random.Random) -> str:
    bases = "ACGT"
    chars = list(seq)
    changes = round(len(seq) * (1.0 - identity))
    for at in rng.sample(range(len(seq)), changes):
        chars[at] = rng.choice([b for b in bases if b != chars[at]])
    return "".join(chars)


def planted_queries(reference: str, query_len: int, count: int, identity: float, seed: int):
    rng = random.Random(seed)
    # Sampling avoids materialising millions of Python integers for whole genomes.
    starts = []
    for _ in range(count):
        for _attempt in range(10000):
            pos = rng.randrange(len(reference) - query_len + 1)
            if "N" not in reference[pos:pos + query_len]:
                starts.append(pos)
                break
        else:
            raise ValueError("Unable to find enough unambiguous A/C/G/T reference regions.")
    return {f"q{i+1}": (mutate_exact(reference[start:start + query_len], identity, rng), start) for i, start in enumerate(starts)}


def write_fasta(path: Path, records: dict[str, str]):
    with path.open("w") as file:
        for name, seq in records.items():
            file.write(f">{name}\n")
            for at in range(0, len(seq), 80):
                file.write(seq[at:at + 80] + "\n")


def run_blast(query_fasta: Path, db_prefix: Path, k: int):
    cmd = ["blastn", "-query", str(query_fasta), "-db", str(db_prefix),
           "-word_size", str(k), "-outfmt", "6 qseqid sstart send pident length bitscore"]
    for key, value in BLAST_SETTINGS.items():
        cmd.extend(["-" + key, str(value)])
    t0 = time.perf_counter()
    result = subprocess.run(cmd, capture_output=True, text=True, check=True, timeout=180)
    return (time.perf_counter() - t0), result.stdout


def parse_hits(output: str):
    for line in output.splitlines():
        if not line.strip():
            continue
        cols = line.split("\t")
        if len(cols) < 6:
            continue
        yield {"qid": cols[0], "start": int(cols[1]), "end": int(cols[2]),
               "pident": float(cols[3]), "length": int(cols[4]), "bitscore": float(cols[5])}


def estimate_identity(pilot_output: str):
    """Length-weighted %identity of highest-bit-score HSP per query."""
    best = {}
    for row in parse_hits(pilot_output):
        q = row["qid"]
        if q not in best or row["bitscore"] > best[q]["bitscore"]:
            best[q] = row
    denom = sum(row["length"] for row in best.values())
    if denom == 0:
        return None, 0
    return sum(row["pident"] * row["length"] for row in best.values()) / denom / 100.0, len(best)


def recovered_query_ids(output: str, queries: dict, min_overlap: int):
    recovered = set()
    for row in parse_hits(output):
        if row["qid"] not in queries:
            continue
        start0 = queries[row["qid"]][1]
        target_l, target_r = start0 + 1, start0 + len(queries[row["qid"]][0])
        hit_l, hit_r = sorted([row["start"], row["end"]])
        overlap = max(0, min(target_r, hit_r) - max(target_l, hit_l) + 1)
        if overlap >= min_overlap:
            recovered.add(row["qid"])
    return recovered


def run_experiment(reference: str, *, query_len=100, count=20, identity=0.90,
                   pilot_n=5, pilot_k=7, lower_k=7, upper_k=15,
                   target=0.90, fixed_ks=(7, 11, 15), seed=42,
                   min_overlap=40, update: Callable[[str], None] | None=None):
    if not available():
        raise RuntimeError("NCBI BLAST+ is not installed. Both 'blastn' and 'makeblastdb' must be on PATH.")
    if not 7 <= lower_k <= upper_k <= 28:
        raise ValueError("Choose k range between 7 and 28 in this experiment.")
    if not 1 <= pilot_n <= count <= 50:
        raise ValueError("Choose 1 <= pilot_n <= query count <= 50.")
    if query_len < max(upper_k, min_overlap) or query_len > 1000:
        raise ValueError("Query length must support the maximum k and minimum overlap, and be <= 1000.")
    if len(reference) < query_len:
        raise ValueError("Reference must be at least as long as the query length.")
    queries = planted_queries(reference, query_len, count, identity, seed)
    rng = random.Random(seed + 1234)
    pilot_ids = rng.sample(list(queries), pilot_n)
    sample = {qid: queries[qid][0] for qid in pilot_ids}
    with tempfile.TemporaryDirectory(prefix="blast_adaptive_") as td:
        wd = Path(td)
        subject_path = wd / "reference.fasta"
        db_prefix = wd / "reference_db"
        all_path = wd / "queries.fasta"
        pilot_path = wd / "pilot.fasta"
        write_fasta(subject_path, {"reference": reference})
        write_fasta(all_path, {qid: val[0] for qid, val in queries.items()})
        write_fasta(pilot_path, sample)
        if update: update("Building a temporary BLAST nucleotide database…")
        build_start = time.perf_counter()
        subprocess.run(["makeblastdb", "-in", str(subject_path), "-dbtype", "nucl", "-out", str(db_prefix)],
                       check=True, capture_output=True, text=True, timeout=180)
        db_build_s = time.perf_counter() - build_start
        if update: update("Step 1/5: Running pilot BLASTN and estimating sequence identity…")
        pilot_s, pilot_out = run_blast(pilot_path, db_prefix, pilot_k)
        p_hat, pilot_detected = estimate_identity(pilot_out)
        # A pre-specified conservative fallback, rather than inventing a pilot estimate.
        used_fallback = p_hat is None
        if used_fallback:
            p_used = 0.70
        else:
            p_used = min(p_hat, 0.999)
        if update: update("Steps 2–3/5: Computing DP probabilities and selecting k…")
        adaptive_k, meets_target, probs = choose_k(p_used, query_len, lower_k, upper_k, target)
        if update: update("Step 4/5: Running final adaptive BLASTN search…")
        adaptive_search_s, adaptive_out = run_blast(all_path, db_prefix, adaptive_k)
        result_rows = [{"Method":"Adaptive k", "k":adaptive_k,
                        "Recovered":len(recovered_query_ids(adaptive_out, queries, min_overlap)),
                        "Queries":len(queries), "Search time (s)":adaptive_search_s,
                        "Total time incl. pilot (s)":pilot_s + adaptive_search_s}]
        if update: update("Step 5/5: Running fixed and random k baselines…")
        for k in sorted(set(fixed_ks)):
            sec, out = run_blast(all_path, db_prefix, k)
            result_rows.append({"Method":f"Fixed k={k}", "k":k,
                                "Recovered":len(recovered_query_ids(out, queries, min_overlap)),
                                "Queries":len(queries), "Search time (s)":sec,
                                "Total time incl. pilot (s)":sec})
        random_k = rng.randint(lower_k, upper_k)
        sec, out = run_blast(all_path, db_prefix, random_k)
        result_rows.append({"Method":"Random k", "k":random_k,
                            "Recovered":len(recovered_query_ids(out, queries, min_overlap)),
                            "Queries":len(queries), "Search time (s)":sec,
                            "Total time incl. pilot (s)":sec})
    for row in result_rows:
        row["Sensitivity"] = row["Recovered"] / row["Queries"]
        row["Time per query (ms)"] = row["Total time incl. pilot (s)"] / row["Queries"] * 1000
        row["Meets measured target"] = row["Sensitivity"] >= target
    return {"rows": result_rows, "p_hat":p_hat, "p_used":p_used,
            "pilot_detected":pilot_detected, "pilot_count":pilot_n,
            "pilot_seconds":pilot_s, "selected_k":adaptive_k,
            "predicted_sensitivity":probs[adaptive_k], "predicted_target_met":meets_target,
            "fallback":used_fallback, "probabilities":probs,
            "db_build_seconds":db_build_s,
            "reference_length":len(reference), "query_count":len(queries),
            "true_identity":identity,"random_k":random_k,
            "min_overlap":min_overlap, "seed":seed,
            "reference_warning":"Synthetic query fragments were generated from the chosen reference and mutated by substitutions only."}
