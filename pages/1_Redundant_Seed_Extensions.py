"""Drop into an existing Streamlit project's pages/ directory."""
import sys
from pathlib import Path
import pandas as pd
import streamlit as st

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path: sys.path.insert(0,str(ROOT))
from redundancy_engine import compare,parse_fasta,synthetic

from theme import init, hero, footer
init('Redundant Seed Extensions')
hero('02','Redundant Seed Extensions','Find overlapping seed hits, avoid repeating selected extensions, and measure the trade-off between saved computation and alignment quality.')
st.caption('Educational, ungapped exact-seed BLAST-style algorithm. Compare an unpruned baseline against diagonal + interval suppression. Not a modification of NCBI BLAST.')
with st.expander('How is real BLAST different from our naive baseline?',expanded=False):
    st.markdown("""**Naive baseline:** extend every exact k-mer hit independently, even when multiple hits describe the same region.

**Production BLAST:** already uses efficient hit processing, masking, X-drop and, in suitable modes, multi-hit filtering or other optimisation. It does **not** simply extend every hit without control.

**Our experiment:** reproduces one established idea, interval-based suppression on a shared ungapped diagonal, to measure its effect against a deliberately naive baseline. It is **not** a new improvement over NCBI BLAST and is not guaranteed to preserve every alignment.""")

@st.cache_data(show_spinner=False)
def read_example():
    p=ROOT/'data'
    q=parse_fasta((p/'ecoli_queries.fasta').read_text())
    s=parse_fasta((p/'ecoli_k12_mg1655_reference.fasta').read_text())
    return q,s

with st.sidebar:
    st.header('Input sequences')
    mode=st.radio('Source',['Synthetic (quick demo)','Bundled E. coli FASTA','Upload FASTA files','Paste DNA'],index=0)
    if mode=='Synthetic (quick demo)':
        kind=st.selectbox('Sequence type',['perfect','repetitive','mutated','clustered'])
        synth_len=st.slider('Sequence length',40,1000,160,20)
        queries,subjects=synthetic(kind,synth_len)
    elif mode=='Bundled E. coli FASTA':
        try: queries,subjects=read_example()
        except FileNotFoundError:
            st.error('Bundled FASTA not present. Upload FASTA files instead.');st.stop()
    elif mode=='Upload FASTA files':
        qf=st.file_uploader('Query FASTA',type=['fa','fasta','fna','txt'],key='qf')
        sf=st.file_uploader('Database FASTA (multiple records allowed)',type=['fa','fasta','fna','txt'],key='sf')
        if not qf or not sf: st.info('Upload both files.');st.stop()
        queries=parse_fasta(qf.getvalue().decode('utf-8-sig'))
        subjects=parse_fasta(sf.getvalue().decode('utf-8-sig'))
    else:
        query_text=st.text_area('Query A/C/G/T','ACGTACGTACGT')
        database_text=st.text_area('Database A/C/G/T, or FASTA with multiple records','>subject_1\nTTACGTACGTACGTAA')
        queries=parse_fasta(query_text); subjects=parse_fasta(database_text)
    if not queries or not subjects: st.warning('Provide at least one query and one subject.');st.stop()
    chosen=st.selectbox('Query',list(range(len(queries))),format_func=lambda ix: f'{ix+1}. {queries[ix][0][:35]} ({len(queries[ix][1]):,} bp)')
    query=queries[chosen][1]
    k=st.number_input('k-mer length',min_value=2,max_value=max(2,len(query)),value=min(11,len(query)),step=1)
    xdrop=st.number_input('X-drop',min_value=0,max_value=100,value=5)
    match=st.number_input('Match reward',min_value=1,max_value=20,value=2)
    penalty=st.number_input('Mismatch penalty (negative)',min_value=-20,max_value=-1,value=-1)
    cap=st.number_input('Max collected seed hits',min_value=1000,max_value=200000,value=120000,step=1000)
    st.caption('Very small k on a large genome can produce huge numbers of hits. The cap prevents excessive processing; capped runs are incomplete.')
    go=st.button('▶ Run comparison',type='primary',use_container_width=True)

st.write(f'**Query:** {queries[chosen][0]} · {len(query):,} bp  |  **Database:** {len(subjects)} sequence(s), {sum(len(v) for _,v in subjects):,} bp')
if not go:
    st.info('Choose inputs and click **Run comparison** in the sidebar. The bundled genome is large; try the synthetic example first.')
    st.stop()
try:
    with st.spinner('Scanning k-mers, then running both extension methods…'):
        run=compare(query,subjects,int(k),int(match),int(penalty),int(xdrop),int(cap))
except Exception as exc:
    st.error(str(exc));st.stop()

base=run['baseline']; opt=run['optimized']; scan=run['scan']; hits=run['hits']
if scan['truncated']:
    st.warning('The seed-hit cap was reached: these are **partial search results**, not full-genome results. Increase k or reduce the database size for a complete experiment.')

st.subheader('1. Side-by-side results')
metrics=[('Seed hits discovered',len(hits),len(hits)),('Seeds extended',base['extended'],opt['extended']),('Seeds skipped',base['skipped'],opt['skipped']),('Extension comparisons',base['comparisons'],opt['comparisons']),('Extension time (ms)',round(base['seconds']*1000,3),round(opt['seconds']*1000,3)),('Best alignment score',base['best'].score if base['best'] else None,opt['best'].score if opt['best'] else None)]
st.dataframe(pd.DataFrame(metrics,columns=['Metric','Baseline','Optimised']),hide_index=True,use_container_width=True)
c1,c2,c3=st.columns(3)
saved=base['comparisons']-opt['comparisons']
c1.metric('Character comparisons saved',f'{saved:,}',f'{100*saved/base["comparisons"]:.1f}% fewer' if base['comparisons'] else 'n/a')
c2.metric('Extensions avoided',f'{opt["skipped"]:,}')
c3.metric('Seed lookup + scan time',f'{scan["lookup_seconds"]*1000:.1f} ms')
if (base['best'].score if base['best'] else None)!=(opt['best'].score if opt['best'] else None):
    st.error('⚠ Best alignment score differs. This pruning heuristic lost a better candidate on this input.')
else: st.success('Both methods returned the same highest score for this test. This does NOT prove general correctness.')
st.caption('Extension times exclude the identical shared seed lookup/scan. Total estimated wall times = shared lookup time + each method’s extension time. Small timing differences are noisy; use repeat runs for benchmarks.')

st.subheader('2. What was skipped, and why?')
rows=[]
for h in opt['trace']:
    rows.append({'Subject':subjects[h.sid][0][:26],'Query seed':h.i,'Subject seed':h.j,'Diagonal j−i':h.diagonal,'Action':h.action,'Best query interval':f'[{h.start_i},{h.end_i})','Score':h.score if h.action=='EXTENDED' else '—','Reason':h.reason})
st.dataframe(pd.DataFrame(rows),hide_index=True,use_container_width=True)
st.caption('Preview of the first 80 decisions in database-scan order. Coordinates are zero-based and intervals are end-exclusive. Skipped rows do not get a newly calculated score.')
if opt['best']:
    b=opt['best'];st.markdown(f'**Best optimised alignment:** subject `{subjects[b.sid][0]}` · query `[{b.start_i},{b.end_i})` · subject `[{b.start_j},{b.end_j})` · score **{b.score}**')
    st.code('Query:   '+query[b.start_i:b.end_i][:180]+'\nSubject: '+subjects[b.sid][1][b.start_j:b.end_j][:180],language='text')

st.subheader('3. Interpreting performance')
st.markdown('**Extension comparisons** count the extra character pairs examined after seed discovery. **Seed scan time** is shared by both algorithms. Compare extension time separately, or add shared scan time for an approximate total. Never compare character counts and seed counts as if they had the same units.')
work=pd.DataFrame({'Algorithm':['Baseline','Optimised'],'Extension comparisons':[base['comparisons'],opt['comparisons']], 'Extensions run':[base['extended'],opt['extended']], 'Extension time (ms)':[base['seconds']*1000,opt['seconds']*1000]})
st.dataframe(work,hide_index=True,use_container_width=True)
st.bar_chart(work.set_index('Algorithm')[['Extension comparisons']],color='#087e8b')
st.caption('The bar chart compares only character comparisons, a consistent unit. For tiny inputs, measured milliseconds are noisy. The same best score on one test does not establish overall sensitivity preservation.')

export=pd.DataFrame(metrics,columns=['Metric','Baseline','Optimised']).to_csv(index=False)
st.download_button('⬇ Download comparison CSV',export,'blast_redundancy_comparison.csv',mime='text/csv')

footer()
