import streamlit as st

CSS = r'''
<style>
@import url('https://fonts.googleapis.com/css2?family=DM+Sans:wght@400;500;600;700&family=Space+Grotesk:wght@500;600;700&display=swap');
/* Force a consistent light interface even when browser or OS prefers dark mode. */
:root { color-scheme: light; --blast-ink:#172539; --blast-muted:#506379; --blast-teal:#176b6a; --blast-paper:#f6f8fc; }
html, body, .stApp, [data-testid="stAppViewContainer"], [data-testid="stMain"] {background-color:var(--blast-paper) !important; color:var(--blast-ink) !important; font-family:'DM Sans',sans-serif !important;}
.stApp, .stApp * {font-family:'DM Sans',sans-serif;}
h1,h2,h3,h4,h5,h6 {font-family:'Space Grotesk',sans-serif !important;letter-spacing:-.025em; color:var(--blast-ink) !important;}
[data-testid="stMarkdownContainer"], [data-testid="stMarkdownContainer"] p, [data-testid="stMarkdownContainer"] li, [data-testid="stMarkdownContainer"] strong, [data-testid="stMarkdownContainer"] label, [data-testid="stCaptionContainer"], [data-testid="stWidgetLabel"] {color:var(--blast-ink) !important;}
[data-testid="stCaptionContainer"], .smallprint, [data-testid="stWidgetLabel"] small {color:var(--blast-muted) !important;}
/* Sidebar navigation + control labels */
[data-testid="stSidebar"], [data-testid="stSidebarContent"] {background:#edf3f8 !important;color:var(--blast-ink) !important;border-right:1px solid #d6e3ed;}
[data-testid="stSidebar"] p, [data-testid="stSidebar"] span, [data-testid="stSidebar"] label, [data-testid="stSidebar"] h1, [data-testid="stSidebar"] h2, [data-testid="stSidebar"] h3, [data-testid="stSidebar"] a {color:var(--blast-ink) !important;}
[data-testid="stSidebarNav"] a, [data-testid="stSidebarNav"] a span, [data-testid="stSidebarNav"] a p {color:var(--blast-ink) !important;font-weight:600;}
[data-testid="stSidebarNav"] a[aria-current="page"], [data-testid="stSidebarNav"] a[aria-current="page"] span, [data-testid="stSidebarNav"] a[aria-current="page"] p {color:#0c5553 !important;}
/* All buttons: dark readable text on pale backgrounds; white text on teal primaries. */
.stButton > button, .stDownloadButton > button, [data-testid="stFormSubmitButton"] button, button[kind="secondary"], button[kind="tertiary"] {background:#fff !important;color:#172539 !important;border:1px solid #bdccd9 !important;border-radius:10px !important;font-weight:700 !important;opacity:1 !important;}
.stButton > button *, .stDownloadButton > button *, [data-testid="stFormSubmitButton"] button *, button[kind="secondary"] *, button[kind="tertiary"] * {color:#172539 !important;}
.stButton > button[kind="primary"], [data-testid="stFormSubmitButton"] button[kind="primary"], .stDownloadButton > button[kind="primary"] {background:#176b6a !important;border-color:#176b6a !important;color:#ffffff !important;}
.stButton > button[kind="primary"] *, [data-testid="stFormSubmitButton"] button[kind="primary"] *, .stDownloadButton > button[kind="primary"] * {color:#ffffff !important;}
.stButton > button:hover:not(:disabled), .stDownloadButton > button:hover:not(:disabled) {border-color:#176b6a !important;box-shadow:0 0 0 1px #176b6a !important;}
button:disabled {opacity:.55 !important;cursor:not-allowed !important;}
/* Readable numbers and labels on result cards. */
[data-testid="stMetric"] {background:#ffffff !important;color:#172539 !important;border:1px solid #dce5ee;border-radius:14px;padding:16px;box-shadow:0 3px 14px rgba(22,46,74,.04);}
[data-testid="stMetric"] [data-testid="stMetricLabel"], [data-testid="stMetric"] [data-testid="stMetricLabel"] *, [data-testid="stMetric"] [data-testid="stMetricValue"], [data-testid="stMetric"] [data-testid="stMetricValue"] *, [data-testid="stMetric"] p {color:#172539 !important;opacity:1 !important;}
[data-testid="stMetricDelta"] [data-testid="stMetricDeltaIcon-Up"], [data-testid="stMetricDelta"] [data-testid="stMetricDeltaIcon-Down"] {color:inherit !important;}
/* Inputs and dropdowns are always light, with legible text. */
[data-baseweb="input"], [data-baseweb="textarea"], [data-baseweb="select"] > div, [data-testid="stNumberInput"] input, [data-testid="stTextInput"] input, [data-testid="stTextArea"] textarea {background:#ffffff !important;color:#172539 !important;}
[data-baseweb="input"] input, [data-baseweb="textarea"] textarea, [data-baseweb="select"] span, [data-baseweb="select"] input, [data-baseweb="select"] div {color:#172539 !important;-webkit-text-fill-color:#172539 !important;}
[data-baseweb="popover"] [role="listbox"], [data-baseweb="popover"] [role="option"] {background:#ffffff !important;color:#172539 !important;}
/* Other surfaces and feedback. */
[data-testid="stExpander"], [data-testid="stForm"] {background:#fff !important;color:#172539 !important;border:1px solid #dce5ee;border-radius:12px;}
[data-testid="stExpander"] summary, [data-testid="stExpander"] summary * {color:#172539 !important;}
[data-testid="stAlert"] {color:#172539 !important;}
[data-testid="stAlert"] p, [data-testid="stAlert"] span {color:#172539 !important;}
/* Hero is intentionally dark. Ensure it keeps white text. */
.hero {background:linear-gradient(115deg,#102c47 0%,#176b6a 76%,#247f84 100%);color:#fff !important;padding:26px 30px;border-radius:20px;margin:4px 0 24px;}
.hero h1,.hero h2,.hero p,.hero .eyebrow {color:#fff !important;}
.hero p {color:#e0f1ef !important;margin:9px 0 0;font-size:1rem;max-width:800px;}
.hero .eyebrow {color:#b8eee3 !important;font-weight:700;text-transform:uppercase;letter-spacing:.18em;font-size:.7rem;}
.note {background:#fff;border-left:4px solid #e0a84d;padding:14px 18px;border-radius:9px;margin:12px 0;color:#172539 !important;}
.tip {background:#e2f3ef;border-left:4px solid #1f8079;padding:14px 18px;border-radius:9px;margin:12px 0;color:#172539 !important;}
.block-container {padding-top:2.1rem;max-width:1240px;}
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
