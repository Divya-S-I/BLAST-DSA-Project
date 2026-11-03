"""Educational ungapped seed-and-extend comparison. Not NCBI BLAST."""
from collections import defaultdict
from dataclasses import dataclass
from time import perf_counter
from typing import List
import random

DNA = set('ACGT')


def parse_fasta(text):
    records=[]; name=None; chunks=[]
    for line in text.splitlines():
        line=line.strip()
        if not line: continue
        if line.startswith('>'):
            if name is not None: records.append((name, ''.join(chunks).upper()))
            name=line[1:].strip() or f'sequence_{len(records)+1}'; chunks=[]
        elif name is None:
            name='sequence_1'; chunks.append(line)
        else: chunks.append(line)
    if name is not None: records.append((name, ''.join(chunks).upper()))
    return [(name,seq) for name,seq in records if seq]


def validate(records):
    for name,seq in records:
        invalid=set(seq)-DNA
        if invalid: raise ValueError(f'{name}: unsupported characters {sorted(invalid)}. Only A/C/G/T accepted. Remove N or ambiguous bases explicitly before analysis.')


def synthetic(kind='perfect', length=180, seed=42):
    rng=random.Random(seed)
    q=''.join(rng.choices('ACGT',k=length))
    if kind=='repetitive': q=('ACGT'*(length//4+1))[:length]
    s=list(q)
    if kind=='mutated':
        for i in range(9,len(s),10):
            s[i]=rng.choice([a for a in 'ACGT' if a!=s[i]])
    if kind=='clustered':
        center=max(0,length//2-5)
        for i in range(center,min(length,center+10)):
            s[i]=rng.choice([a for a in 'ACGT' if a!=s[i]])
    return [('synthetic_query',q)],[('synthetic_subject',''.join(s))]


def lookup_query(query,k):
    out=defaultdict(list)
    for i in range(len(query)-k+1): out[query[i:i+k]].append(i)
    return out


def find_hits(query,subjects,k,max_hits=150000):
    """One shared seed scan for both methods; deterministic scan order, exact hits only."""
    t0=perf_counter(); lookup=lookup_query(query,k)
    hits=[]; total_windows=0
    for sid,(name,sub) in enumerate(subjects):
        for j in range(max(0,len(sub)-k+1)):
            total_windows+=1
            for i in lookup.get(sub[j:j+k],()):
                hits.append((sid,i,j))
                if len(hits)>=max_hits:
                    return hits, {'lookup_seconds':perf_counter()-t0,'windows':total_windows,'truncated':True,'word_keys':len(lookup)}
    return hits, {'lookup_seconds':perf_counter()-t0,'windows':total_windows,'truncated':False,'word_keys':len(lookup)}


@dataclass
class HitResult:
    sid:int
    i:int
    j:int
    diagonal:int
    start_i:int
    end_i:int
    start_j:int
    end_j:int
    score:int
    comparisons:int
    action:str
    reason:str


def extend(query,sub,i,j,k,match,mismatch,xdrop):
    """Independent left/right ungapped X-drop, best score includes seed only once."""
    score=k*match; comparisons=0
    best_left=(0,0); running=0; best=0
    for t in range(1,min(i,j)+1):
        comparisons+=1
        running+=match if query[i-t]==sub[j-t] else mismatch
        if running>best: best=running; best_left=(t,best)
        if best-running>xdrop: break
    left,ls=best_left
    right_bound=min(len(query)-(i+k),len(sub)-(j+k))
    best_right=(0,0); running=0; best=0
    for t in range(1,right_bound+1):
        comparisons+=1
        running+=match if query[i+k+t-1]==sub[j+k+t-1] else mismatch
        if running>best: best=running; best_right=(t,best)
        if best-running>xdrop: break
    right,rs=best_right
    return (i-left, i+k+right, j-left, j+k+right, score+ls+rs, comparisons)


def run_method(query,subjects,hits,k,match=2,mismatch=-1,xdrop=5,optimized=False,examples=80, collect_limit=200):
    """Pruning heuristic: on a diagonal, skip seeds contained in previously accepted best-scoring interval."""
    t0=perf_counter()
    coverage=defaultdict(list)
    count_extended=0; count_skipped=0; comparisons=0
    best=None; trace=[]; results=[]
    for sid,i,j in hits:
        sub=subjects[sid][1]; d=j-i; key=(sid,d)
        covering=None
        if optimized:
            for start,end in coverage[key]:
                if start<=i and i+k<=end:
                    covering=(start,end); break
        if covering is not None:
            count_skipped+=1
            if len(trace)<examples:
                trace.append(HitResult(sid,i,j,d,i,i+k,j,j+k,0,0,'SKIPPED',f'Seed [{i},{i+k}) contained in prior best interval [{covering[0]},{covering[1]}) on diagonal {d}'))
            continue
        ai,bi,aj,bj,score,comp=extend(query,sub,i,j,k,match,mismatch,xdrop)
        comparisons+=comp; count_extended+=1
        hit=HitResult(sid,i,j,d,ai,bi,aj,bj,score,comp,'EXTENDED',f'X-drop extension scored {score}; aligned query [{ai},{bi})')
        if best is None or hit.score>best.score: best=hit
        if len(results)<collect_limit: results.append(hit)
        if len(trace)<examples: trace.append(hit)
        if optimized:
            # merging intervals on same diagonal, sorted not necessary for membership; avoid merging adjacent unrelated regions
            coverage[key].append((ai,bi))
    return {'seconds':perf_counter()-t0,'extended':count_extended,'skipped':count_skipped,'comparisons':comparisons,'best':best,'trace':trace,'results':results,'coverage':coverage}


def compare(query,subjects,k,match=2,mismatch=-1,xdrop=5,max_hits=150000):
    if not (1<=k<=len(query)): raise ValueError('k must be between 1 and the query length')
    if match<=0 or mismatch>=0 or xdrop<0: raise ValueError('Use match > 0, mismatch < 0, X-drop >= 0')
    validate([('query',query)]+subjects)
    hits,scan=find_hits(query,subjects,k,max_hits)
    base=run_method(query,subjects,hits,k,match,mismatch,xdrop,False)
    opt=run_method(query,subjects,hits,k,match,mismatch,xdrop,True)
    return {'hits':hits,'scan':scan,'baseline':base,'optimized':opt}


def recovered(result,subject_idx,start,end):
    """Checks exact known planted interval coverage; requires one retained best result."""
    hit=result['best']
    return bool(hit and hit.sid==subject_idx and hit.start_j<=start and hit.end_j>=end)
