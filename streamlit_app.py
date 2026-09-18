import os
import sys
import time
import base64
import pandas as pd
import streamlit as st

# RDKit imports for molecular processing & structure rendering
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors, Lipinski

try:
    from rdkit.Chem.Draw import rdMolDraw2D
    HAS_RDKIT_DRAW = True
except Exception:
    rdMolDraw2D = None
    HAS_RDKIT_DRAW = False

# Import predictor engine
from app.predictor import nik_predictor_instance, calculate_physicochemical_properties

# Helper function to convert local image to Base64
def get_image_base64(relative_path: str) -> str:
    abs_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), relative_path)
    if os.path.exists(abs_path):
        with open(abs_path, "rb") as img_file:
            return base64.b64encode(img_file.read()).decode("utf-8")
    return ""

pup_logo_b64 = get_image_base64(os.path.join("app", "static", "pup_logo.png"))
lab_logo_b64 = get_image_base64(os.path.join("app", "static", "lab_logo_round.jpeg"))
workflow_img_b64 = get_image_base64(os.path.join("app", "static", "nik_workflow.jpeg"))

# Streamlit Page Config
st.set_page_config(
    page_title="NIK-Predict: NF-κB Inducing Kinase Bioactivity Prediction Platform",
    page_icon="🧪",
    layout="wide",
    initial_sidebar_state="collapsed"
)

# Custom CSS matching exact style.css from localhost:8000
css_code = """
<style>
@import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700;800&family=JetBrains+Mono:wght@400;500;600&display=swap');

:root {
  --pu-red: #9e1b1e;
  --pu-red-dark: #7b1114;
  --pu-red-light: #c2292d;
  --pu-red-soft: #fbf0f0;
  --pu-gold: #c59b27;
  --slate-900: #0f172a;
  --slate-800: #1e293b;
  --slate-700: #334155;
  --slate-600: #475569;
  --slate-500: #64748b;
  --slate-400: #94a3b8;
  --slate-300: #cbd5e1;
  --slate-200: #e2e8f0;
  --slate-100: #f1f5f9;
  --slate-50: #f8fafc;
  --white: #ffffff;
}

/* Hide Streamlit elements */
#MainMenu {visibility: hidden;}
footer {visibility: hidden;}
header {visibility: hidden;}
.stApp > header {display: none;}
[data-testid="stSidebar"] {display: none;}
.block-container {
    padding-top: 0rem !important;
    padding-bottom: 2rem !important;
    max-width: 1240px !important;
}

body {
  font-family: 'Inter', sans-serif !important;
  background-color: #f6f8fa !important;
  color: var(--slate-800) !important;
}

/* Header */
.site-header {
  background: var(--white);
  border-bottom: 1px solid var(--slate-200);
  box-shadow: 0 1px 2px 0 rgb(0 0 0 / 0.05);
  margin-bottom: 1.5rem;
}

.header-container {
  padding: 1rem 0;
  display: flex;
  align-items: center;
  justify-content: space-between;
  flex-wrap: wrap;
  gap: 1rem;
}

.brand-identity {
  display: flex;
  align-items: center;
  gap: 1rem;
}

.university-logo {
  height: 58px;
  width: auto;
  object-fit: contain;
}

.lab-logo-img {
  height: 58px;
  width: auto;
  object-fit: contain;
  border-radius: 50%;
  border: 1px solid var(--slate-200);
}

.brand-text .lab-title {
  font-size: 1.25rem;
  font-weight: 800;
  color: var(--pu-red);
  margin: 0;
  line-height: 1.2;
}

.brand-text .dept-title {
  font-size: 0.9rem;
  font-weight: 600;
  color: var(--slate-700);
  margin: 0;
  line-height: 1.3;
}

.brand-text .university-title {
  font-size: 0.82rem;
  font-weight: 500;
  color: var(--slate-500);
  margin: 0;
}

.contact-group {
  display: flex;
  flex-direction: column;
  gap: 0.35rem;
  align-items: flex-end;
}

.contact-pill {
  display: inline-flex;
  align-items: center;
  gap: 0.45rem;
  background: var(--pu-red-soft);
  color: var(--pu-red-dark);
  font-size: 0.82rem;
  font-weight: 600;
  padding: 0.4rem 0.95rem;
  border-radius: 9999px;
  text-decoration: none;
  border: 1px solid rgba(158, 27, 30, 0.2);
}

/* Hero Section */
.hero-section {
  text-align: center;
  max-width: 900px;
  margin: 1.5rem auto;
}

.hero-tag {
  display: inline-block;
  font-size: 0.72rem;
  letter-spacing: 0.08em;
  font-weight: 700;
  color: var(--pu-red);
  background: var(--pu-red-soft);
  padding: 0.25rem 0.75rem;
  border-radius: 9999px;
  margin-bottom: 0.75rem;
  text-transform: uppercase;
  border: 1px solid rgba(158, 27, 30, 0.15);
}

.hero-title {
  font-size: 2.1rem;
  font-weight: 800;
  color: var(--slate-900);
  margin-bottom: 0.75rem;
}

.hero-desc {
  font-size: 1.05rem;
  color: var(--slate-600);
  line-height: 1.55;
}

/* NIK Info Card */
.nik-info-card {
  background: var(--white);
  border: 1px solid var(--slate-200);
  border-radius: 16px;
  padding: 1.5rem 1.75rem;
  box-shadow: 0 1px 3px rgba(0,0,0,0.05);
  margin-bottom: 1.5rem;
}

.nik-info-grid {
  display: grid;
  grid-template-columns: 1fr 340px;
  gap: 1.5rem;
  align-items: center;
}

.nik-info-heading {
  font-size: 1.15rem;
  font-weight: 800;
  color: var(--pu-red);
  margin-bottom: 0.6rem;
}

.nik-info-para {
  font-size: 0.88rem;
  color: var(--slate-700);
  line-height: 1.6;
  margin-bottom: 0.6rem;
}

.nik-workflow-wrap {
  display: flex;
  justify-content: center;
  align-items: center;
  background: var(--slate-50);
  border: 1px solid var(--slate-200);
  border-radius: 12px;
  padding: 0.75rem;
}

.nik-workflow-img {
  max-width: 100%;
  height: auto;
  border-radius: 8px;
}

/* Dropdown styling */
.model-select-wrap {
  background: #f1f5f9;
  border: 1px solid var(--slate-200);
  border-radius: 8px;
  padding: 0.85rem 1rem;
  margin-bottom: 1.25rem;
}

.model-select-label {
  font-size: 0.82rem;
  font-weight: 700;
  color: var(--slate-700);
  text-transform: uppercase;
  letter-spacing: 0.03em;
}

/* Property Table */
.properties-card {
  background: var(--white);
  border: 1px solid var(--slate-200);
  border-radius: 12px;
  padding: 1rem 1.25rem;
  box-shadow: 0 1px 2px rgba(0,0,0,0.05);
  margin-top: 1rem;
}

.properties-title {
  font-size: 0.88rem;
  font-weight: 700;
  color: var(--slate-800);
  text-transform: uppercase;
  margin-bottom: 0.75rem;
}

.properties-table {
  width: 100%;
  border-collapse: collapse;
  font-size: 0.85rem;
}

.properties-table th, .properties-table td {
  padding: 0.55rem 0.85rem;
  text-align: left;
  border-bottom: 1px solid var(--slate-100);
}

.properties-table th {
  background: var(--slate-50);
  font-weight: 700;
  color: var(--slate-600);
  font-size: 0.75rem;
  text-transform: uppercase;
}

.properties-table td.prop-val {
  font-weight: 600;
  color: var(--slate-900);
  font-family: 'JetBrains Mono', monospace;
}

/* Consensus Banner */
.consensus-banner {
  background: linear-gradient(135deg, #1e293b 0%, #0f172a 100%);
  border-radius: 16px;
  padding: 2rem;
  color: var(--white);
  text-align: center;
  box-shadow: 0 10px 15px -3px rgba(0,0,0,0.1);
  margin-top: 1.5rem;
}

.consensus-label {
  font-size: 0.74rem;
  font-weight: 700;
  letter-spacing: 0.08em;
  color: var(--pu-gold);
  text-transform: uppercase;
}

.consensus-value {
  font-size: 3.5rem;
  font-weight: 800;
  color: var(--white);
  font-family: 'JetBrains Mono', monospace;
  margin: 0.5rem 0;
}

.consensus-subtext {
  font-size: 0.85rem;
  color: var(--slate-400);
}

/* Footer */
.site-footer {
  background: var(--white);
  border-top: 1px solid var(--slate-200);
  padding: 2rem 0;
  margin-top: 3rem;
}

.footer-container {
  display: flex;
  justify-content: space-between;
  align-items: center;
  flex-wrap: wrap;
  gap: 1.5rem;
  font-size: 0.85rem;
  color: var(--slate-600);
}

.footer-left .footer-lab {
  font-weight: 800;
  color: var(--pu-red);
}
</style>
"""
st.markdown(css_code, unsafe_allow_html=True)

# Preset Curated Test Compounds
NIK_SAMPLES = [
    {
        "name": "Compound 1",
        "smiles": "COc1cnc(nc1N1CCc2c1cc(Br)cc2)N",
        "description": "Aminopyrimidine derivative (Experimental pIC50 = 5.0693)"
    },
    {
        "name": "Lead NIK015",
        "smiles": "N#Cc1ccc(cc1)c1cnc(s1)C(=O)Nc1nccc(n1)n1cnc2c1ccc(c2)Cl",
        "description": "Nitrile benzimidazole analog (Experimental pIC50 = 6.2254)"
    },
    {
        "name": "Benchmark NIK009",
        "smiles": "Nc1nc(N2CCc3c2cc(OC)cc3)c(Cl)cn1",
        "description": "Methoxy chloropyrimidine analog (Experimental pIC50 = 7.0500)"
    },
    {
        "name": "Heteroaryl NIK017",
        "smiles": "Nc1nc(N2CCc3c2cc(c2n[nH]cc2)cc3)c(Cl)cn1",
        "description": "Pyrazole-substituted analog (Experimental pIC50 = 8.1900)"
    }
]

# Render Header HTML with dual logos and contact pills
header_html = f"""
<header class="site-header">
  <div class="header-container">
    <div class="brand-identity">
      <img src="data:image/png;base64,{pup_logo_b64}" class="university-logo" alt="Punjabi University Logo" />
      <img src="data:image/jpeg;base64,{lab_logo_b64}" class="lab-logo-img" alt="Drug Design & Synthesis Lab" />
      <div class="brand-text">
        <h2 class="lab-title">Drug Design &amp; Synthesis Lab</h2>
        <h3 class="dept-title">Department of Pharmaceutical Sciences and Drug Research</h3>
        <h4 class="university-title">Punjabi University, Patiala, Punjab, India</h4>
      </div>
    </div>
    <div class="contact-group">
      <a href="mailto:drugdesignsynthesislab@gmail.com" class="contact-pill">✉ drugdesignsynthesislab@gmail.com</a>
      <a href="mailto:Yogita_pharma@pbi.ac.in" class="contact-pill">✉ Yogita_pharma@pbi.ac.in</a>
    </div>
  </div>
</header>
"""
st.markdown(header_html, unsafe_allow_html=True)

# Top Navigation Bar Tabs
if "active_nav" not in st.session_state:
    st.session_state["active_nav"] = "Main"

nav_options = ["Main", "About", "What is NIK?", "Dataset", "Model performance", "Collaboration & Contact", "Limitations"]
cols_nav = st.columns(len(nav_options))
for idx, opt in enumerate(nav_options):
    btn_style = "primary" if st.session_state["active_nav"] == opt else "secondary"
    if cols_nav[idx].button(opt, key=f"nav_btn_{idx}", type=btn_style, use_container_width=True):
        st.session_state["active_nav"] = opt
        st.rerun()

st.markdown("---")

# Navigation Routing
nav = st.session_state["active_nav"]

# Render Hero Banner if on Main or About
if nav in ["Main", "About"]:
    st.markdown("""
    <section class="hero-section">
      <div class="hero-tag">2D-QSAR MONTE CARLO PREDICTION PLATFORM</div>
      <h1 class="hero-title">NIK-Predict: NF-κB Inducing Kinase Bioactivity Prediction Platform</h1>
      <p class="hero-desc">
        Predict the biological inhibitory activity (pIC<sub>50</sub>) of small molecules against 
        <strong>NF-κB Inducing Kinase (NIK / MAP3K14)</strong> using machine learning and 2D-QSAR modeling.
      </p>
    </section>
    """, unsafe_allow_html=True)

if nav in ["Main", "What is NIK?"]:
    workflow_html = f"""
    <section class="nik-info-card">
      <div class="nik-info-grid">
        <div>
          <h3 class="nik-info-heading">What is NIK (NF-κB Inducing Kinase)?</h3>
          <p class="nik-info-para">
            <strong>NIK</strong> is a serine/threonine kinase belonging to the MAP3K family and serves as a central regulator of the noncanonical NF-κB pathway. Its activity is normally kept low through continuous ubiquitination and proteasomal degradation mediated by the TRAF–cIAP E3 ligase complex.
          </p>
          <p class="nik-info-para">
            When this regulatory control is disrupted by mutations or altered expression, NIK accumulates and becomes pathologically active. Elevated NIK drives excessive processing of <strong>p100 to p52</strong>, leading to overactivation of downstream transcriptional programs. This aberrant signalling promotes the production of pro-inflammatory cytokines and survival factors, contributing to <em>autoimmune diseases, chronic inflammatory conditions, B-cell malignancies, immune dysfunction, and tissue injury</em>.
          </p>
        </div>
        <div class="nik-workflow-wrap">
          <img src="data:image/jpeg;base64,{workflow_img_b64}" class="nik-workflow-img" alt="NIK Workflow" />
        </div>
      </div>
    </section>
    """
    st.markdown(workflow_html, unsafe_allow_html=True)

# Main Predictor Workspace
if nav == "Main":
    # Workspace Container
    tab_single, tab_batch = st.tabs(["Single Compound", "Batch Prediction"])

    with tab_single:
        # Test Compound Chips
        st.markdown("**TEST COMPOUNDS:**")
        sample_cols = st.columns(len(NIK_SAMPLES))
        if "input_smiles" not in st.session_state:
            st.session_state["input_smiles"] = NIK_SAMPLES[0]["smiles"]

        for idx, sample in enumerate(NIK_SAMPLES):
            if sample_cols[idx].button(sample["name"], key=f"chip_sample_{idx}", use_container_width=True):
                st.session_state["input_smiles"] = sample["smiles"]
                st.rerun()

        # Model Selector
        st.markdown('<div class="model-select-wrap"><span class="model-select-label">CHOOSE A PREDICTION MODEL</span></div>', unsafe_allow_html=True)
        st.selectbox("Model", ["NIK 2D QSAR MODEL"], index=0, label_visibility="collapsed")

        col_left, col_right = st.columns([1, 1])

        with col_left:
            user_smiles = st.text_input(
                "ENTER SMILES STRING",
                value=st.session_state["input_smiles"],
                key="smiles_input_field",
                placeholder="e.g. COc1cnc(nc1N1CCc2c1cc(Br)cc2)N"
            )
            st.caption("Paste canonical SMILES or click a test compound. 2D structure renders instantly on the right.")

            col_b1, col_b2 = st.columns([1, 1])
            with col_b1:
                btn_predict = st.button("Predict Bioactivity", type="primary", use_container_width=True)
            with col_b2:
                btn_reset = st.button("Reset", use_container_width=True)
                if btn_reset:
                    st.session_state["input_smiles"] = NIK_SAMPLES[0]["smiles"]
                    st.rerun()

        with col_right:
            st.markdown("##### 2D MOLECULAR STRUCTURE")
            clean_s = user_smiles.strip()
            if clean_s:
                props = calculate_physicochemical_properties(clean_s)
                if props and props.get("svg_structure"):
                    st.markdown('<span style="color: #059669; font-weight: 700; font-size: 0.8rem;">✔ 2D Chemical Structure Valid</span>', unsafe_allow_html=True)
                    st.components.v1.html(props["svg_structure"], height=250, scrolling=False)
                else:
                    st.markdown('<span style="color: #dc2626; font-weight: 700; font-size: 0.8rem;">✖ Invalid SMILES String</span>', unsafe_allow_html=True)
            else:
                st.info("Structure will render upon SMILES entry")

            # Physicochemical Properties Table
            if clean_s and 'props' in locals() and props:
                st.markdown("""
                <div class="properties-card">
                  <div class="properties-title">Physicochemical Properties</div>
                  <table class="properties-table">
                    <thead>
                      <tr><th>Property</th><th>Value</th></tr>
                    </thead>
                    <tbody>
                      <tr><td>Formula</td><td class="prop-val">{}</td></tr>
                      <tr><td>Molecular Weight</td><td class="prop-val">{} g/mol</td></tr>
                      <tr><td>LogP (Lipophilicity)</td><td class="prop-val">{}</td></tr>
                      <tr><td>TPSA</td><td class="prop-val">{} Å²</td></tr>
                      <tr><td>H-Bond Donors / Acceptors</td><td class="prop-val">{}</td></tr>
                      <tr><td>Rotatable Bonds</td><td class="prop-val">{}</td></tr>
                    </tbody>
                  </table>
                </div>
                """.format(
                    props["formula"],
                    props["molecular_weight"],
                    props["logp"],
                    props["tpsa"],
                    props["h_donors_acceptors"],
                    props["rotatable_bonds"]
                ), unsafe_allow_html=True)

        if btn_predict:
            if not clean_s:
                st.error("Please enter a valid SMILES string.")
            else:
                with st.spinner("Calculating consensus 2D-QSAR prediction across Monte Carlo runs..."):
                    try:
                        res_dict = nik_predictor_instance.predict_single(clean_s)
                        consensus_val = res_dict["consensus_prediction"]
                        elapsed = res_dict["elapsed_seconds"]

                        st.markdown(f"""
                        <div class="consensus-banner">
                          <div class="consensus-label">CONSENSUS PREDICTED pIC<sub>50</sub></div>
                          <div class="consensus-value">{consensus_val:.4f}</div>
                          <div class="consensus-subtext">Calculated consensus biological inhibitory potency against NIK kinase (Computed in {elapsed}s)</div>
                        </div>
                        """, unsafe_allow_html=True)
                    except Exception as e:
                        st.error(f"Prediction Error: {str(e)}")

    with tab_batch:
        st.markdown("#### Batch SMILES Library Evaluation")
        batch_text = st.text_area(
            "SMILES List (One compound per line)",
            value="COc1cnc(nc1N1CCc2c1cc(Br)cc2)N\nN#Cc1ccc(cc1)c1cnc(s1)C(=O)Nc1nccc(n1)n1cnc2c1ccc(c2)Cl\nNc1nc(N2CCc3c2cc(OC)cc3)c(Cl)cn1",
            height=150
        )
        btn_batch = st.button("Run Batch Prediction", type="primary")
        if btn_batch:
            clean_lines = [s.strip() for s in batch_text.splitlines() if s.strip()]
            if not clean_lines:
                st.error("Please enter at least one valid SMILES string.")
            else:
                with st.spinner(f"Evaluating {len(clean_lines)} compounds..."):
                    try:
                        batch_res = nik_predictor_instance.predict_batch(clean_lines)
                        results_list = batch_res["results"]
                        table_data = []
                        for idx, item in enumerate(results_list):
                            table_data.append({
                                "Index": idx + 1,
                                "SMILES": item["smiles"],
                                "Consensus pIC50": f"{item['consensus_prediction']:.4f}"
                            })
                        df_results = pd.DataFrame(table_data)
                        st.markdown(f"### Batch Results ({len(results_list)} Compounds Evaluated)")
                        st.dataframe(df_results, use_container_width=True)
                        csv_data = df_results.to_csv(index=False).encode('utf-8')
                        st.download_button(
                            label="📥 Export CSV Results",
                            data=csv_data,
                            file_name="NIK_Consensus_Predictions.csv",
                            mime="text/csv"
                        )
                    except Exception as e:
                        st.error(f"Batch Prediction Error: {str(e)}")

elif nav == "Dataset":
    st.markdown("""
    <div class="nik-info-card">
      <h3 class="nik-info-heading">Dataset Overview</h3>
      <p class="nik-info-para">
        Calculated on a curated dataset of <strong>118 NIK inhibitors</strong> collected from published literature, standardized to experimental pIC<sub>50</sub> values.
      </p>
    </div>
    """, unsafe_allow_html=True)

elif nav == "Model performance":
    st.markdown("""
    <div class="nik-info-card">
      <h3 class="nik-info-heading">Model Performance Overview</h3>
      <p class="nik-info-para">
        Four Monte Carlo 2D-QSAR models (M1–M4) were developed using the CORALSEA correlation balance approach. 
        <strong>Model M3</strong> was selected as the final model based on its overall statistical performance and lower prediction error.
      </p>
      <div style="display: flex; gap: 1.5rem; margin-top: 1rem;">
        <div style="background: #f8fafc; border: 1px solid #cbd5e1; padding: 1rem; border-radius: 8px; text-align: center; flex: 1;">
          <span style="font-size: 0.75rem; font-weight: 700; color: #64748b;">VALIDATION R²</span>
          <div style="font-size: 1.8rem; font-weight: 800; color: #9e1b1e;">0.765</div>
        </div>
        <div style="background: #f8fafc; border: 1px solid #cbd5e1; padding: 1rem; border-radius: 8px; text-align: center; flex: 1;">
          <span style="font-size: 0.75rem; font-weight: 700; color: #64748b;">VALIDATION Q²</span>
          <div style="font-size: 1.8rem; font-weight: 800; color: #9e1b1e;">0.660</div>
        </div>
        <div style="background: #f8fafc; border: 1px solid #cbd5e1; padding: 1rem; border-radius: 8px; text-align: center; flex: 1;">
          <span style="font-size: 0.75rem; font-weight: 700; color: #64748b;">VALIDATION MAE</span>
          <div style="font-size: 1.8rem; font-weight: 800; color: #9e1b1e;">0.291</div>
        </div>
      </div>
    </div>
    """, unsafe_allow_html=True)

elif nav == "Collaboration & Contact":
    st.markdown("""
    <div class="nik-info-card">
      <h3 class="nik-info-heading">Collaboration &amp; Contact</h3>
      <p class="nik-info-para">
        We are open to scientific collaboration, co-development of tools, and data exchange projects in computational drug discovery, molecular docking, dynamic and QSAR modeling.
      </p>
      <p class="nik-info-para"><strong>Lab Email:</strong> <a href="mailto:drugdesignsynthesislab@gmail.com">drugdesignsynthesislab@gmail.com</a></p>
      <p class="nik-info-para"><strong>Faculty Email:</strong> <a href="mailto:Yogita_pharma@pbi.ac.in">Yogita_pharma@pbi.ac.in</a></p>
      <hr style="margin: 1.5rem 0; border: 0; border-top: 1px dashed #cbd5e1;" />
      <h4 style="font-size: 1rem; font-weight: 700; color: #0f172a; margin-bottom: 0.5rem;">📖 How to Cite</h4>
      <p class="nik-info-para">If you use <strong>NIK-Predict</strong> in your research, please cite us:<br/>
      <code style="background: #fbf0f0; color: #9e1b1e; padding: 0.3rem 0.6rem; border-radius: 4px; font-weight: 600;">[Citation link will be provided after publication]</code></p>
    </div>
    """, unsafe_allow_html=True)

elif nav == "Limitations":
    st.markdown("""
    <div class="nik-info-card" style="background: #fffbeb; border-color: rgba(217, 119, 6, 0.3);">
      <h3 class="nik-info-heading" style="color: #92400e;">⚠️ Model Limitations</h3>
      <ul style="color: #78350f; font-size: 0.9rem; line-height: 1.6; margin-left: 1.25rem;">
        <li><strong>Small Molecules:</strong> Applicable strictly to small-molecule inhibitors. Not tested for macrocycles, peptides, or prodrugs.</li>
        <li><strong>Chemical Space:</strong> Predictions may be less reliable outside the model's chemical space domain.</li>
        <li><strong>Experimental Validation:</strong> In-silico predictions provide preliminary bioactivity estimates; experimental validation is required.</li>
      </ul>
    </div>
    """, unsafe_allow_html=True)

# Footer
footer_html = """
<footer class="site-footer">
  <div class="footer-container">
    <div class="footer-left">
      <p class="footer-lab">Drug Design &amp; Synthesis Lab (In-Silico to In-Vivo)</p>
      <p class="footer-dept">Department of Pharmaceutical Sciences and Drug Research</p>
      <p class="footer-inst">Punjabi University, Patiala, Punjab, India</p>
    </div>
    <div class="footer-right">
      <p class="footer-contact">
        Inquiries: 
        <a href="mailto:drugdesignsynthesislab@gmail.com">drugdesignsynthesislab@gmail.com</a> | 
        <a href="mailto:Yogita_pharma@pbi.ac.in">Yogita_pharma@pbi.ac.in</a>
      </p>
      <p class="footer-copy">&copy; 2026 Drug Design &amp; Synthesis Lab, Punjabi University Patiala. All rights reserved.</p>
    </div>
  </div>
</footer>
"""
st.markdown(footer_html, unsafe_allow_html=True)
