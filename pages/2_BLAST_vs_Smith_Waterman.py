"""Module 3: analysis of uploaded protein experiment and small live DNA demo."""
import sys
from pathlib import Path
import pandas as pd
import streamlit as st
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from theme import init,hero,footer
from redundancy_engine import compare
init('BLAST vs Smith-Waterman')
hero('03','BLAST vs Smith-Waterman','An evidence-based comparison of heuristic search and exact local alignment: speed, score recovery, and the limits of each method.')
st.markdown('The **protein experiment** below uses the uploaded saved outputs (192 within-gene pairs). The optional **live DNA demo** runs our educational ungapped seed-and-extend algorithm, not NCBI BLASTP.')
@st.cache_data

def load_data():
    b=pd.read_csv(ROOT/'data'/'blastp_optimized_results.csv')
    s=pd.read_csv(ROOT/'data'/'smith_waterman_results.csv')
    keys=['gene','species_a','species_b']
    merged=b.merge(s,on=keys,suffixes=('_blast','_sw'),validate='one_to_one')
    merged['score_ratio']=merged['raw_score_blast']/merged['raw_score_sw']
    merged['score_match']=merged['raw_score_blast'].sub(merged['raw_score_sw']).abs()<1e-8
    merged['speedup']=merged['runtime_seconds_sw']/merged['runtime_seconds_blast']
    return merged
try: df=load_data()
except Exception as e:
    st.error(f'Cannot load the experiment result files: {e}');st.stop()
with st.sidebar:
    st.header('Experiment filters')
    gene=st.selectbox('Protein family',['All genes']+sorted(df.gene.unique().tolist()))
    mode=st.radio('Display',['Saved protein results','Interactive DNA demonstration'])
if mode=='Saved protein results':
    d=df if gene=='All genes' else df[df.gene==gene]
    c1,c2,c3,c4=st.columns(4)
    c1.metric('Protein pairs',f'{len(d):,}')
    c2.metric('BLAST-style total time',f'{d.runtime_seconds_blast.sum():.2f} s')
    c3.metric('Smith-Waterman total time',f'{d.runtime_seconds_sw.sum():.2f} s')
    c4.metric('Total-runtime speedup',f'{d.runtime_seconds_sw.sum()/d.runtime_seconds_blast.sum():.2f}×')
    e1,e2=st.columns(2)
    e1.metric('Mean per-pair score recovery',f'{100*d.score_ratio.mean():.1f}%')
    e2.metric('Identical raw scores',f'{int(d.score_match.sum())} / {len(d)}')
    st.caption('Mean score recovery averages BLAST/SW per comparison; total-runtime speedup uses the ratio of summed run times. Both depend on the supplied implementation and scoring settings.')
    st.subheader('Performance across protein families')
    groups=[]
    for g,x in d.groupby('gene'):
        groups.append({'Gene':g,'Pairs':len(x),'BLAST time (s)':x.runtime_seconds_blast.sum(),'SW time (s)':x.runtime_seconds_sw.sum(),'Speedup':x.runtime_seconds_sw.sum()/x.runtime_seconds_blast.sum(),'Mean recovery %':x.score_ratio.mean()*100,'Exact score ties':int(x.score_match.sum()),'BLAST > SW':int((x.raw_score_blast>x.raw_score_sw).sum())})
    summary=pd.DataFrame(groups)
    st.dataframe(summary.style.format({'BLAST time (s)':'{:.2f}','SW time (s)':'{:.2f}','Speedup':'{:.2f}×','Mean recovery %':'{:.1f}%'}),use_container_width=True,hide_index=True)
    left,right=st.columns(2)
    with left:
        st.markdown('**Total runtime per family (seconds)**')
        st.bar_chart(summary.set_index('Gene')[['BLAST time (s)','SW time (s)']])
    with right:
        st.markdown('**Mean BLAST score as % of Smith-Waterman**')
        st.bar_chart(summary.set_index('Gene')[['Mean recovery %']])
    st.subheader('Inspect a particular sequence pair')
    names=[f'{r.gene} · {r.species_a} vs {r.species_b}' for r in d.itertuples()]
    selected=st.selectbox('Sequence pair',range(len(d)),format_func=lambda i:names[i])
    row=d.iloc[selected]
    c1,c2,c3=st.columns(3)
    c1.metric('BLAST raw score',f'{row.raw_score_blast:.2f}')
    c2.metric('Smith-Waterman raw score',f'{row.raw_score_sw:.2f}')
    c3.metric('Score recovery',f'{100*row.score_ratio:.1f}%')
    st.dataframe(pd.DataFrame([{'Algorithm':'BLAST-style','Runtime (seconds)':row.runtime_seconds_blast,'Alignment length':row.alignment_length_blast,'Percent identity':row.identity_percentage_blast},{'Algorithm':'Smith-Waterman','Runtime (seconds)':row.runtime_seconds_sw,'Alignment length':row.alignment_length_sw,'Percent identity':row.identity_percentage_sw}]),use_container_width=True,hide_index=True)
    st.download_button('Download filtered comparison CSV',d.to_csv(index=False),'blast_vs_smith_waterman.csv','text/csv')
    disagree=d[d.raw_score_blast>d.raw_score_sw]
    if len(disagree):
        st.warning(f'{len(disagree)} pairs in this filter have BLAST-style scores higher than Smith-Waterman. Under the same scoring rules and valid alignment definitions, this should not occur. Review score normalisation, gap penalties, implementation and coordinate conventions before claiming a rigorous accuracy comparison.')
    st.markdown('<div class="note"><strong>Interpretation:</strong> These are saved <em>protein</em> experiment outputs. They compare the recorded implementations, not necessarily equivalent production BLASTP vs an independently validated Smith-Waterman baseline. Raw-score recovery is not the same thing as biological sensitivity or accuracy.</div>',unsafe_allow_html=True)
else:
    st.subheader('Small, live DNA example')
    c1,c2=st.columns(2)
    with c1: q=st.text_input('Query (ACGT only)','ACGTACGT')
    with c2: s=st.text_input('Subject (ACGT only)','TTACGTTCGT')
    k=st.slider('k-mer size',2,8,3)
    match=st.number_input('Match reward',1,10,2)
    mismatch=st.number_input('Mismatch penalty',-10,-1,-1)
    gap=st.number_input('Gap penalty',-10,-1,-2)
    xdrop=st.slider('X-drop threshold',0,12,5)
    def smith_waterman(a,b,match,mismatch,gap):
        H=[[0]*(len(b)+1) for _ in range(len(a)+1)]
        best=(0,0,0)
        for i in range(1,len(a)+1):
            for j in range(1,len(b)+1):
                v=max(0,H[i-1][j-1]+(match if a[i-1]==b[j-1] else mismatch),H[i-1][j]+gap,H[i][j-1]+gap)
                H[i][j]=v
                if v>best[0]:best=(v,i,j)
        return best[0],len(a)*len(b)
    if st.button('Run both algorithms',type='primary'):
        q=q.upper().strip();s=s.upper().strip()
        if len(q)>350 or len(s)>350 or not q or not s or set(q+s)-set('ACGT') or k>min(len(q),len(s)):
            st.error('Use DNA strings with A/C/G/T only, length ≤350 and k no greater than the shorter sequence.')
        else:
            sw_score,cells=smith_waterman(q,s,match,mismatch,gap)
            comparison=compare(q,[('subject',s)],k,match,mismatch,xdrop)
            hit=comparison['baseline']['best']
            a,b=st.columns(2)
            a.metric('Smith-Waterman optimal local score',sw_score)
            b.metric('Seed-and-extend best score',hit.score if hit else 0)
            st.caption(f'Smith-Waterman evaluated {cells:,} DP cells; the BLAST-style baseline found {len(comparison["hits"]):,} seed hits and used {comparison["baseline"]["comparisons"]:,} additional extension comparisons. These counts measure different operations.')
            if hit:
                st.code(f'Query   {q[hit.start_i:hit.end_i]}\nSubject {s[hit.start_j:hit.end_j]}',language='text')
            st.info('Smith-Waterman allows gaps and finds the optimal local score. The seed heuristic may miss alignments, particularly those with no exact k-mer. This is an educational comparison, not a benchmark of production BLAST.')
footer()
