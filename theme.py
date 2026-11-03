import streamlit as st

CSS = '''
<style>
@import url('https://fonts.googleapis.com/css2?family=DM+Sans:wght@400;500;600;700&family=Space+Grotesk:wght@500;600;700&display=swap');
html,body,[class*="css"],.stApp {font-family:'DM Sans',sans-serif;}
h1,h2,h3 {font-family:'Space Grotesk',sans-serif; letter-spacing:-.035em;}
.stApp {background: #f6f8fc; color:#172539;}
[data-testid="stSidebar"] {background: #edf3f8; border-right: 1px solid #d6e3ed;}
[data-testid="stMetric"] {background:#ffffff; border:1px solid #dce5ee; border-radius:14px; padding:16px; box-shadow: 0 3px 14px rgba(22,46,74,.04);}
.stButton>button[kind="primary"] {background:#176b6a; border-color:#176b6a; color:white; border-radius:10px; font-weight:700;}
.stButton>button {border-radius:10px;}
[data-testid="stExpander"] {background:white; border:1px solid #dce5ee; border-radius:12px;}
.block-container {padding-top:2.1rem; max-width:1240px;}
.hero {background:linear-gradient(115deg,#102c47 0%,#176b6a 76%,#247f84 100%); color:#fff; padding:26px 30px;border-radius:20px;margin:4px 0 24px;}
.hero h1 {color:white; margin:0; font-size:2.15rem;}
.hero p {color:#d4e9e8; margin:9px 0 0; font-size:1rem; max-width:800px;}
.eyebrow {color:#a9e1d7; font-weight:700; text-transform:uppercase;letter-spacing:.18em;font-size:.7rem;}
.note {background:#fff; border-left:4px solid #e0a84d; padding:14px 18px; border-radius:9px; margin:12px 0;}
.tip {background:#e2f3ef; border-left:4px solid #1f8079; padding:14px 18px; border-radius:9px; margin:12px 0;}
.smallprint {color:#586b7e;font-size:.85rem;}
</style>
'''

def init(title):
    st.set_page_config(page_title=f'{title} | BLAST Lab',page_icon='🧬',layout='wide',initial_sidebar_state='expanded')
    st.markdown(CSS,unsafe_allow_html=True)

def hero(number,title,subtitle):
    st.markdown(f'<div class="hero"><div class="eyebrow">BLAST LAB · MODULE {number}</div><h1>{title}</h1><p>{subtitle}</p></div>',unsafe_allow_html=True)

def footer():
    st.divider()
    st.caption('BLAST Lab · Educational algorithms and controlled experiments. BLAST-style prototypes are not the NCBI BLAST executable. Coordinate conventions are documented on each page.')
