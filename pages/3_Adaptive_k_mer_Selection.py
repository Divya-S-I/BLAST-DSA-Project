"""Interactive, reproducible mathematical component of supplied adaptive-k study."""
import sys,random
from pathlib import Path
import pandas as pd
import streamlit as st
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from theme import init,hero,footer
init('Adaptive k-mer Selection')
hero('04','Adaptive k-mer Selection','Choose the largest exact word size that meets a predicted hit-sensitivity target, then study the expected random-hit workload.')
@st.cache_data(show_spinner=False)
def hit_probability(k:int,p:float,N:int)->float:
    if k>N:return 0.0
    f=[0.0]*k; f[0]=1.0
    for _ in range(N):
        new=[0.0]*k
        for j,val in enumerate(f):
            new[0]+=val*(1-p)
            if j+1<k:new[j+1]+=val*p
        f=new
    return max(0.0,min(1.0,1-sum(f)))
def expected_random_hits(k,m,n,both_strands):
    if k>min(m,n):return 0.0
    return (2 if both_strands else 1)*(m-k+1)*(n-k+1)/(4**k)
with st.sidebar:
    st.header('Selection controls')
    p=st.slider('Estimated per-position identity p',0.50,1.00,0.90,0.01)
    N=st.slider('Target region length N',30,500,100,10)
    m=st.number_input('Query length m',min_value=30,max_value=100000,value=100,step=10)
    n=st.number_input('Total database length n',min_value=1000,max_value=100000000,value=4641652,step=1000)
    min_k=st.number_input('Minimum k',min_value=2,max_value=30,value=8)
    max_k=st.number_input('Maximum k',min_value=int(min_k),max_value=40,value=max(28,int(min_k)))
    target=st.slider('Target predicted sensitivity',0.50,0.999,0.90,0.01)
    both_strands=st.checkbox('Approximate both DNA strands',value=True)
ks=list(range(int(min_k),int(max_k)+1))
rows=[{'k':k,'Predicted sensitivity':hit_probability(k,p,N),'Expected random seed hits':expected_random_hits(k,int(m),int(n),both_strands)} for k in ks]
df=pd.DataFrame(rows)
valid=df[df['Predicted sensitivity']>=target]
selected=int(valid['k'].max()) if len(valid) else int(min_k)
selected_row=df.loc[df.k==selected].iloc[0]
c1,c2,c3=st.columns(3)
c1.metric('Selected k',selected)
c2.metric('Predicted sensitivity',f'{selected_row["Predicted sensitivity"]:.1%}')
c3.metric('Expected random position-pair hits',f'{selected_row["Expected random seed hits"]:,.2f}')
if not len(valid):st.warning('No k in the range reaches the chosen target. Fallback uses the minimum k; the target is **not** satisfied.')
st.markdown('### 1. Mathematical selection rule')
st.latex(r'S(k,p,N)=P(\text{at least one run of }k\text{ matches in }N\text{ positions})')
st.latex(r'k^*=\max\{k:S(k,\hat p,N)\geq \tau\}')
st.caption('Dynamic programming tracks the probability of having no k-run so far and the current run length j = 0, …, k−1. A mismatch resets j to zero; a match advances j. Probabilities reaching k leave the no-hit states.')
st.markdown('### 2. Sensitivity versus seed workload')
st.line_chart(df.set_index('k')[['Predicted sensitivity']],color='#087e8b')
st.line_chart(df.set_index('k')[['Expected random seed hits']],color='#d69b38')
st.caption('The random-hit formula assumes independent, uniform DNA; actual genomic repeats and base composition change the count. The sensitivity model assumes independent aligned matches and does not model gaps, BLAST filters or E-values.')
st.dataframe(df.assign(**{'Meets target':df['Predicted sensitivity']>=target}).style.format({'Predicted sensitivity':'{:.2%}','Expected random seed hits':'{:,.3f}'}),hide_index=True,use_container_width=True)
# -------------------- Full experimental strategy -----------------------
st.divider()
st.markdown('### 3. Actual BLASTN experiment: five-step adaptive strategy')
st.write('This is a **real NCBI BLASTN experiment**, separate from the mathematical simulator above. It builds a nucleotide database, creates queries from known locations in your chosen reference, performs a pilot search, chooses k, and measures recovery and execution time against fixed and random k values.')
st.caption('Planted query fragments are mutated by substitutions to a known identity. Measured sensitivity means recovering their known reference location, not proving biological homology.')
from adaptive_blastn_experiment import available as blast_plus_available, fasta_to_reference, run_experiment

if not blast_plus_available():
    st.warning('NCBI BLAST+ is not installed or not on your PATH. The mathematical model above still works. To run this section locally on macOS, install BLAST+ using `brew install blast` (with Homebrew installed), restart Streamlit, and confirm `blastn -version` and `makeblastdb -version` in Terminal.')

with st.form('real_blastn_experiment'):
    ref_choice=st.radio('Reference sequence', ['Bundled E. coli genome', 'Paste your own A/C/G/T sequence', 'Upload a single-record FASTA'], horizontal=True)
    pasted_ref=''
    uploaded_ref=None
    if ref_choice=='Paste your own A/C/G/T sequence':
        pasted_ref=st.text_area('Reference DNA (at least as long as the query)', value='ACGT'*300, height=85)
    elif ref_choice=='Upload a single-record FASTA':
        uploaded_ref=st.file_uploader('Reference FASTA', type=['fa','fasta','fna','txt'], key='adaptive_reference_upload')
    s1,s2,s3=st.columns(3)
    with s1:
        true_identity=st.slider('Planted sequence identity',0.65,1.0,0.90,0.01)
        query_length=st.number_input('Query length (bp)',min_value=50,max_value=500,value=100,step=10)
    with s2:
        num_queries=st.number_input('Queries per method',min_value=5,max_value=40,value=20,step=5)
        pilot_count=st.number_input('Pilot sample queries',min_value=1,max_value=5,value=5)
    with s3:
        real_target=st.slider('Required predicted sensitivity',0.50,0.99,0.90,0.01,key='real_target')
        experiment_seed=st.number_input('Reproducibility seed',min_value=0,max_value=100000,value=42,key='real_seed')
    st.markdown('**Experimental configuration:** pilot k = 7; candidate k = 7–15; fixed baselines k = 7, 11, 15; plus one randomly chosen k. BLASTN uses reward +2, mismatch −3, gap-open 5, gap-extend 2, E-value threshold 0.001, one CPU thread. Both strands are searched.')
    execute=st.form_submit_button('Run actual five-step BLASTN experiment',type='primary',disabled=not blast_plus_available())

if execute:
    try:
        if ref_choice=='Bundled E. coli genome':
            reference=fasta_to_reference((ROOT/'data'/'ecoli_k12_mg1655_reference.fasta').read_text())
        elif ref_choice=='Paste your own A/C/G/T sequence':
            reference=fasta_to_reference(pasted_ref)
        else:
            if uploaded_ref is None:
                raise ValueError('Upload a reference FASTA before starting.')
            if uploaded_ref.size > 12_000_000:
                raise ValueError('Please provide a FASTA smaller than 12 MB for this interactive experiment.')
            reference=fasta_to_reference(uploaded_ref.getvalue().decode('utf-8'))
        if len(reference)<int(query_length):
            raise ValueError('Your reference is shorter than the query length.')
        progress=st.empty()
        def progress_message(message):progress.info(message)
        with st.spinner('Running NCBI BLASTN against the selected reference...'):
            outcome=run_experiment(reference,query_len=int(query_length),count=int(num_queries),
                identity=float(true_identity),pilot_n=int(pilot_count),pilot_k=7,
                lower_k=7,upper_k=15,target=float(real_target),
                fixed_ks=(7,11,15),seed=int(experiment_seed),min_overlap=min(40,int(query_length)),
                update=progress_message)
        progress.empty()
        st.session_state['adaptive_real_result']=outcome
    except Exception as exc:
        st.error(f'Experiment could not complete: {exc}')

if 'adaptive_real_result' in st.session_state:
    outcome=st.session_state['adaptive_real_result']
    st.markdown('#### Results from the actual BLASTN run')
    a,b,c,d=st.columns(4)
    a.metric('Pilot-estimated identity', 'No pilot hits' if outcome['p_hat'] is None else f"{outcome['p_hat']:.1%}")
    b.metric('Selected word size',str(outcome['selected_k']))
    c.metric('Predicted sensitivity',f"{outcome['predicted_sensitivity']:.1%}")
    d.metric('Pilot runtime',f"{outcome['pilot_seconds']*1000:.1f} ms")
    if outcome['fallback']:
        st.warning('The pilot found no alignments. For this test a predefined p=0.70 fallback was used. The estimated identity is unknown, so the sensitivity prediction should not be considered calibrated.')
    if not outcome['predicted_target_met']:
        st.warning('None of the candidate k values met the prediction target. The experiment used the minimum k=7 as a fallback.')
    stats=pd.DataFrame(outcome['rows'])
    st.dataframe(stats[['Method','k','Recovered','Queries','Sensitivity','Time per query (ms)','Meets measured target']].style.format({'Sensitivity':'{:.1%}','Time per query (ms)':'{:.2f}'}), hide_index=True,use_container_width=True)
    t1,t2=st.columns(2)
    with t1:
        st.markdown('**Measured true-location recovery**')
        st.bar_chart(stats.set_index('Method')['Sensitivity'], color='#087e8b')
    with t2:
        st.markdown('**Time per query (including pilot for adaptive k)**')
        st.bar_chart(stats.set_index('Method')['Time per query (ms)'], color='#d69b38')
    pred_df=pd.DataFrame([{'k':kk,'Predicted sensitivity':val} for kk,val in outcome['probabilities'].items()])
    with st.expander('Predicted sensitivity for every candidate k'):
        st.dataframe(pred_df.style.format({'Predicted sensitivity':'{:.2%}'}),hide_index=True,use_container_width=True)
    st.download_button('Download actual BLASTN comparison CSV',stats.to_csv(index=False),file_name='adaptive_k_actual_blastn_results.csv',mime='text/csv')
    st.caption(f"Actual BLASTN measurements from {outcome['query_count']} planted {int(query_length)}-bp queries against {outcome['reference_length']:,} reference bases; random seed {outcome['seed']}. Each method is run once. The fixed/random times exclude database construction; the adaptive time includes pilot + final search. Shared database build: {outcome['db_build_seconds']:.2f}s (excluded from all methods). Pilot estimates may be biased, and the independence-based DP model ignores gaps and complex mutation structure. These results are not a claim that adaptive k always outperforms fixed k.")

footer()
