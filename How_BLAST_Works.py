"""Module 1: visual and interactive walkthrough of simplified DNA BLAST."""
from collections import defaultdict
import sys
from pathlib import Path
import pandas as pd
import streamlit as st
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
from theme import init,hero,footer
from redundancy_engine import lookup_query,find_hits,extend
init('How BLAST Works')
hero('01','How BLAST Works','An interactive, step-by-step investigation of the seed-and-extend idea, from DNA strings to high-scoring local alignments.')
st.caption('Modules: How BLAST Works · Redundant Seed Extensions · BLAST vs Smith-Waterman · Adaptive k-mer Selection')
with st.sidebar:
    st.header('🔬 Experiment controls')
    query=st.text_input('DNA query',value='TTACGTACGTAC')
    subject=st.text_input('DNA database sequence',value='GGGTTACGTTCGTACAAAA')
    k=st.slider('Word size k',2,8,3)
    match=st.number_input('Match reward',value=2,min_value=1,max_value=10)
    mismatch=st.number_input('Mismatch penalty',value=-1,min_value=-10,max_value=-1)
    xdrop=st.slider('X-drop',0,12,5)
query=query.upper().replace(' ','');subject=subject.upper().replace(' ','')
if not query or not subject or set(query+subject)-set('ACGT'):
    st.error('Please enter nonempty A/C/G/T DNA strings.');st.stop()
if k>min(len(query),len(subject)):
    st.warning('The chosen k is longer than one of the sequences. Lower k.');st.stop()
steps=['1 · Input & preparation','2 · Generate overlapping k-mers','3 · Build a hash lookup','4 · Scan database & locate seeds','5 · Extend and apply X-drop','6 · Score and compare candidates','7 · Rank the results']
st.subheader('Explore the algorithm')
# Explicit one-step navigation; buttons and direct jumps share session state.
if 'blast_stage' not in st.session_state:
    st.session_state.blast_stage = 0
st.session_state.blast_stage = min(max(int(st.session_state.blast_stage), 0), len(steps)-1)
nav_left,nav_middle,nav_right=st.columns([1,3,1])
with nav_left:
    if st.button('← Previous',disabled=st.session_state.blast_stage==0,use_container_width=True):
        st.session_state.blast_stage-=1
with nav_middle:
    st.markdown(f'<p style="text-align:center;font-weight:700;margin:7px 0">Stage {st.session_state.blast_stage+1} of {len(steps)}: {steps[st.session_state.blast_stage].split(" · ",1)[-1]}</p>',unsafe_allow_html=True)
with nav_right:
    if st.button('Next →',disabled=st.session_state.blast_stage==len(steps)-1,use_container_width=True):
        st.session_state.blast_stage+=1
st.caption('Jump directly to any stage:')
for offset in (0,4):
    part=steps[offset:offset+4]
    cols=st.columns(len(part))
    for z,lab in enumerate(part):
        num=offset+z
        with cols[z]:
            if st.button(f'{num+1}. {lab.split(" · ",1)[-1]}',key=f'stage_jump_{num}',type='primary' if st.session_state.blast_stage==num else 'secondary',use_container_width=True):
                st.session_state.blast_stage=num
idx=st.session_state.blast_stage
lookup=lookup_query(query,k)
hits,scan=find_hits(query,[('subject',subject)],k,max_hits=25000)
hit_rows=[{'Query start':i,'Database start':j,'Word':query[i:i+k],'Diagonal j−i':j-i} for _,i,j in hits]
if idx==0:
    st.markdown('**Goal:** find local similarities between a shorter query and a database. This demo uses exact DNA words and ungapped extensions. Real BLAST has additional modes and statistics.')
    a,b=st.columns(2)
    a.code('QUERY\n'+query);b.code('DATABASE\n'+subject)
    st.info(f'Query length m = {len(query)} · Database length n = {len(subject)} · Each contains overlapping words of length k = {k}.')
elif idx==1:
    st.markdown(f'Slide a window of length **{k}** one position at a time. The query has **{len(query)-k+1}** overlapping words.')
    df=pd.DataFrame([{'Start position':i,'k-mer':query[i:i+k]} for i in range(len(query)-k+1)])
    st.dataframe(df,hide_index=True,use_container_width=True,height=230)
    st.latex(r'\mathrm{Number\ of\ }k\mathrm{\!-mers}=m-k+1')
elif idx==2:
    st.markdown('**Hash map:** key = exact word, value = *all query starting positions*. Duplicates are retained as lists of positions.')
    st.json(dict(lookup))
    st.caption('Typical hash-table lookups are expected O(1) for a fixed-size encoded word. Slicing and hashing a length-k Python string costs O(k).')
elif idx==3:
    st.markdown('Slide the same k-base window across the database; look up each word in the query dictionary. A successful lookup becomes a **seed hit**.')
    database_words=[{'Database start':j,'Word':subject[j:j+k],'Found in lookup?':subject[j:j+k] in lookup,'Matching query starts':', '.join(map(str,lookup.get(subject[j:j+k],[])))} for j in range(len(subject)-k+1)]
    st.dataframe(pd.DataFrame(database_words),use_container_width=True,hide_index=True,height=240)
    st.metric('Matching position pairs (seeds)',len(hits))
    st.caption('A word can occur at several query positions. Every matching pair of coordinates counts as a distinct seed hit.')
elif idx==4:
    if not hits: st.warning('No exact seed hits. Try smaller k or different input.')
    else:
        hit_index=st.slider('Choose a seed',0,min(len(hits)-1,199),0)
        _,i,j=hits[hit_index]
        ai,bi,aj,bj,score,comp=extend(query,subject,i,j,k,int(match),int(mismatch),int(xdrop))
        st.markdown(f'Seed **`{query[i:i+k]}`** begins at query **{i}**, database **{j}**. Extend left and right by comparing aligned characters and keeping the best score in each direction.')
        c1,c2,c3=st.columns(3);c1.metric('Final score',score);c2.metric('Extra comparisons',comp);c3.metric('Diagonal',j-i)
        a=query[ai:bi];b=subject[aj:bj]
        marker=''.join('|' if x==y else '·' for x,y in zip(a,b))
        st.code(f'Query  [{ai}:{bi}]    {a}\n                     {marker}\nSubject[{aj}:{bj}]    {b}',language='text')
        st.latex(r'\mathrm{Stop\ when:}\;S_{best}-S_{current}>X')
        st.caption('This implementation extends both sides separately, retains each direction’s highest-scoring endpoint and combines their scores. It does not support insertions/deletions.')
elif idx==5:
    st.markdown('For each candidate seed, count matches and mismatches along the retained alignment. Larger scores indicate stronger alignments *under the chosen scoring scheme*.')
    rows=[]
    for _,i,j in hits[:1000]:
        ai,bi,aj,bj,score,comp=extend(query,subject,i,j,k,int(match),int(mismatch),int(xdrop))
        rows.append({'Seed':query[i:i+k],'q start':i,'db start':j,'Alignment score':score,'Length':bi-ai,'Comparisons':comp})
    if rows: st.dataframe(pd.DataFrame(rows).sort_values('Alignment score',ascending=False),hide_index=True,use_container_width=True,height=260)
    else: st.warning('No exact seeds available.')
else:
    rows=[]
    for _,i,j in hits[:2000]:
        ai,bi,aj,bj,score,comp=extend(query,subject,i,j,k,int(match),int(mismatch),int(xdrop))
        rows.append({'Query interval':f'[{ai},{bi})','Database interval':f'[{aj},{bj})','Score':score,'Seed':query[i:i+k]})
    if rows:
        st.dataframe(pd.DataFrame(rows).sort_values('Score',ascending=False).drop_duplicates(['Query interval','Database interval']).head(15),hide_index=True,use_container_width=True)
    else: st.warning('No seeds found, so this heuristic has no candidate to report.')
    st.markdown('**Production BLAST adds more:** protein neighbourhood words, two-hit filters in some modes, gapped alignment, bit scores and E-values. This demonstration intentionally focuses on exact-seed, ungapped DNA alignment.')
st.divider()
a,b=st.columns(2)
with a:
    st.markdown('#### Mathematical insight')
    st.latex(r'P(\mathrm{random\ exact\ seed})=4^{-k}')
    st.latex(r'E[h]=\frac{(m-k+1)(n-k+1)}{4^k}')
with b:
    st.markdown('#### What to notice')
    st.markdown('Increase **k** to reduce random seeds, but some true matching regions will no longer contain a long enough exact word. Smaller **X-drop** reduces extension work, but can also stop before later matches.')
footer()
