# -*- coding: utf-8 -*-
import os
import sys
import base64
import json
import streamlit as st

# Streamlit Page Config
st.set_page_config(
    page_title="NIK-Predict: NF-κB Inducing Kinase Bioactivity Prediction Platform",
    page_icon="🧪",
    layout="wide",
    initial_sidebar_state="collapsed"
)

# Custom CSS to force full width and remove all Streamlit margins/paddings
st.markdown("""
<style>
#MainMenu {visibility: hidden;}
footer {visibility: hidden;}
header {visibility: hidden;}
.stApp > header {display: none;}
[data-testid="stSidebar"] {display: none;}
.block-container {
    padding: 0rem !important;
    margin: 0rem !important;
    max-width: 100% !important;
}
iframe {
    width: 100% !important;
    border: none !important;
    overflow: auto !important;
}
</style>
""", unsafe_allow_html=True)

# Function to read local file as string/base64
def read_file(relative_path: str) -> str:
    abs_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), relative_path)
    if os.path.exists(abs_path):
        with open(abs_path, "r", encoding="utf-8") as f:
            return f.read()
    return ""

def get_image_base64(relative_path: str) -> str:
    abs_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), relative_path)
    if os.path.exists(abs_path):
        with open(abs_path, "rb") as img_file:
            return base64.b64encode(img_file.read()).decode("utf-8")
    return ""

# Read assets
pup_logo_b64 = get_image_base64(os.path.join("app", "static", "pup_logo.png"))
lab_logo_b64 = get_image_base64(os.path.join("app", "static", "lab_logo_round.jpeg"))
workflow_img_b64 = get_image_base64(os.path.join("app", "static", "nik_workflow.jpeg"))
smiles_drawer_js = read_file(os.path.join("app", "static", "vendor", "smiles-drawer.min.js"))
app_js = read_file(os.path.join("app", "static", "app.js"))
style_css = read_file(os.path.join("app", "static", "style.css"))
index_html = read_file(os.path.join("app", "static", "index.html"))

# Read JSON lookup data dynamically at startup
lookup_file_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "exact_predictions_lookup.json")
exact_predictions_json = "{}"
if os.path.exists(lookup_file_path):
    with open(lookup_file_path, "r", encoding="utf-8") as f:
        lookup_data = json.load(f)
        # Clean any raw newlines inside SVG strings
        for k, v in lookup_data.items():
            if isinstance(v, dict) and "svg" in v and isinstance(v["svg"], str):
                v["svg"] = v["svg"].replace("\r", "").replace("\n", " ")
        exact_predictions_json = json.dumps(lookup_data)

# Construct Single-File Embedded HTML App
# 1. Replace static image paths with Base64 Data URIs
index_html = index_html.replace('/static/pup_logo.png', f'data:image/png;base64,{pup_logo_b64}')
index_html = index_html.replace('/static/lab_logo_round.jpeg', f'data:image/jpeg;base64,{lab_logo_b64}')
index_html = index_html.replace('/static/nik_workflow.jpeg', f'data:image/jpeg;base64,{workflow_img_b64}')

# 2. Embed CSS inline
index_html = index_html.replace('<link rel="stylesheet" href="/static/style.css" />', f'<style>\n{style_css}\n</style>')

# 3. Embed SmilesDrawer JS inline
index_html = index_html.replace('<script src="/static/vendor/smiles-drawer.min.js"></script>', f'<script>\n{smiles_drawer_js}\n</script>')

# 4. Construct JS Engine
embedded_engine_js = f"""
<script src="https://unpkg.com/openchemlib@8.6.0/dist/openchemlib-full.js"></script>
<script>
const EXACT_PREDICTIONS = {exact_predictions_json};
</script>
<script>
{app_js}
</script>
"""

# Replace app.js script tag with full embedded_engine_js
index_html = index_html.replace('<script src="/static/app.js"></script>', embedded_engine_js)

# Render Single-File Embedded HTML App in Streamlit Cloud
st.components.v1.html(index_html, height=2600, scrolling=True)
