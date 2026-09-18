import os
import sys
import time
import pandas as pd
import streamlit as st
from PIL import Image

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

# Streamlit Page Config
st.set_page_config(
    page_title="NIK-Predict: NF-κB Inducing Kinase Bioactivity Prediction Platform",
    page_icon="🧪",
    layout="wide",
    initial_sidebar_state="expanded"
)

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

# Helper: Render RDKit SVG Drawing in Streamlit
def render_rdkit_svg(smiles_str):
    if not HAS_RDKIT_DRAW or rdMolDraw2D is None:
        return None
    try:
        mol = Chem.MolFromSmiles(smiles_str)
        if mol is None:
            return None
        drawer = rdMolDraw2D.MolDraw2DSVG(400, 260)
        opts = drawer.drawOptions()
        opts.clearBackground = False
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        return drawer.GetDrawingText()
    except Exception:
        return None

# Sidebar Navigation & Lab Branding
st.sidebar.markdown("### 🔬 Drug Design & Synthesis Lab")
st.sidebar.markdown("*In-Silico to In-Vivo*")
st.sidebar.markdown("**Department of Pharmaceutical Sciences and Drug Research**  \n**Punjabi University, Patiala, Punjab, India**")

nav = st.sidebar.radio(
    "Navigation",
    ["Main Predictor", "About", "What is NIK?", "Dataset", "Model Performance", "Collaboration & Contact", "Limitations"]
)

st.sidebar.markdown("---")
st.sidebar.markdown("📧 **Inquiries:**")
st.sidebar.markdown("`drugdesignsynthesislab@gmail.com`")
st.sidebar.markdown("`Yogita_pharma@pbi.ac.in`")

# Header Section with Dual Logos
col_header1, col_header2 = st.columns([1, 4])
with col_header1:
    pup_logo_path = os.path.join("app", "static", "pup_logo.png")
    lab_logo_path = os.path.join("app", "static", "lab_logo_round.jpeg")
    
    if os.path.exists(pup_logo_path):
        st.image(pup_logo_path, width=110)
    if os.path.exists(lab_logo_path):
        st.image(lab_logo_path, width=110)

with col_header2:
    st.markdown("## NIK-Predict: NF-κB Inducing Kinase Bioactivity Prediction Platform")
    st.markdown("##### Drug Design & Synthesis Lab | Punjabi University, Patiala, Punjab, India")
    st.caption("Predict experimental biological inhibitory potency (pIC50) of small molecules targeting NF-κB Inducing Kinase (NIK / MAP3K14)")

st.markdown("---")

# Navigation Routing
if nav == "Main Predictor":
    st.markdown("### 🎯 2D-QSAR Bioactivity Prediction Workspace")
    
    # Model Selector Dropdown
    model_choice = st.selectbox(
        "Choose a prediction model",
        ["NIK 2D QSAR MODEL"],
        index=0
    )
    
    tab_single, tab_batch = st.tabs(["Single Compound Prediction", "Batch Prediction"])
    
    # -------------------------------------------------------------
    # Single Compound Tab
    # -------------------------------------------------------------
    with tab_single:
        st.markdown("#### Test Compounds (Click to load sample):")
        sample_cols = st.columns(len(NIK_SAMPLES))
        
        if "input_smiles" not in st.session_state:
            st.session_state["input_smiles"] = NIK_SAMPLES[0]["smiles"]
            
        for idx, sample in enumerate(NIK_SAMPLES):
            if sample_cols[idx].button(sample["name"], key=f"btn_sample_{idx}"):
                st.session_state["input_smiles"] = sample["smiles"]
                
        user_smiles = st.text_input(
            "Enter candidate SMILES string",
            value=st.session_state["input_smiles"],
            key="smiles_input_field"
        )
        
        col_input_left, col_input_right = st.columns([1, 1])
        
        with col_input_right:
            st.markdown("##### 2D Chemical Structure")
            if user_smiles.strip():
                svg_data = render_rdkit_svg(user_smiles.strip())
                if svg_data:
                    st.components.v1.html(svg_data, height=270, scrolling=False)
                    st.success("2D Chemical Structure Valid")
                else:
                    st.error("Invalid SMILES input string")
            else:
                st.info("Structure will render upon SMILES entry")

        with col_input_left:
            st.markdown("##### Physicochemical Properties")
            if user_smiles.strip():
                props = calculate_physicochemical_properties(user_smiles.strip())
                if props:
                    prop_df = pd.DataFrame([
                        {"Property": "Formula", "Value": props["formula"]},
                        {"Property": "Molecular Weight", "Value": f"{props['molecular_weight']} g/mol"},
                        {"Property": "LogP (Lipophilicity)", "Value": props["logp"]},
                        {"Property": "TPSA", "Value": f"{props['tpsa']} Å²"},
                        {"Property": "H-Bond Donors / Acceptors", "Value": props["h_donors_acceptors"]},
                        {"Property": "Rotatable Bonds", "Value": props["rotatable_bonds"]}
                    ])
                    st.table(prop_df)
                else:
                    st.warning("Could not calculate physicochemical properties for this SMILES.")

        btn_predict = st.button("Predict Bioactivity (pIC50)", type="primary")
        
        if btn_predict:
            if not user_smiles.strip():
                st.error("Please enter a valid SMILES string.")
            else:
                with st.spinner("Calculating consensus 2D-QSAR prediction across Monte Carlo runs..."):
                    try:
                        res_dict = nik_predictor_instance.predict_single(user_smiles.strip())
                        consensus_val = res_dict["consensus_prediction"]
                        elapsed = res_dict["elapsed_seconds"]
                        
                        st.markdown("### Consensus Prediction Result")
                        st.metric(
                            label="CONSENSUS PREDICTED pIC50",
                            value=f"{consensus_val:.4f}",
                            delta=f"Computed in {elapsed}s"
                        )
                        st.caption("Calculated consensus biological inhibitory potency against NIK kinase")
                    except Exception as e:
                        st.error(f"Prediction Error: {str(e)}")

    # -------------------------------------------------------------
    # Batch Prediction Tab
    # -------------------------------------------------------------
    with tab_batch:
        st.markdown("#### Batch SMILES Library Evaluation")
        batch_text = st.text_area(
            "Enter SMILES list (One compound per line)",
            value="COc1cnc(nc1N1CCc2c1cc(Br)cc2)N\nN#Cc1ccc(cc1)c1cnc(s1)C(=O)Nc1nccc(n1)n1cnc2c1ccc(c2)Cl\nNc1nc(N2CCc3c2cc(OC)cc3)c(Cl)cn1",
            height=150
        )
        
        btn_batch = st.button("Run Batch Prediction", type="primary")
        
        if btn_batch:
            clean_lines = [s.strip() for s in batch_text.splitlines() if s.strip()]
            if not clean_lines:
                st.error("Please enter at least one valid SMILES string.")
            else:
                with st.spinner(f"Evaluating {len(clean_lines)} compounds in parallel..."):
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

elif nav == "About":
    st.markdown("### 🏛️ About NIK-Predict")
    st.write(
        "An open-source, high-throughput computational platform for predicting the biological inhibitory activity "
        "(pIC50) of small molecules against NF-κB Inducing Kinase (NIK) using a validated 3-run consensus Monte Carlo 2D-QSAR model."
    )

elif nav == "What is NIK?":
    st.markdown("### 🧬 What is NIK (NF-κB Inducing Kinase)?")
    col1, col2 = st.columns([3, 2])
    with col1:
        st.write(
            "**NIK** is a serine/threonine kinase belonging to the MAP3K family and serves as a central regulator of the noncanonical NF-κB pathway. "
            "Its activity is normally kept low through continuous ubiquitination and proteasomal degradation mediated by the TRAF–cIAP E3 ligase complex."
        )
        st.write(
            "When this regulatory control is disrupted by mutations or altered expression, NIK accumulates and becomes pathologically active. "
            "Elevated NIK drives excessive processing of **p100 to p52**, leading to overactivation of downstream transcriptional programs. "
            "This aberrant signalling promotes the production of pro-inflammatory cytokines and survival factors, contributing to "
            "*autoimmune diseases, chronic inflammatory conditions, B-cell malignancies, immune dysfunction, and tissue injury*."
        )
    with col2:
        workflow_img_path = os.path.join("app", "static", "nik_workflow.jpeg")
        if os.path.exists(workflow_img_path):
            st.image(workflow_img_path, caption="NIK-Predict Workflow Diagram", use_column_width=True)

elif nav == "Dataset":
    st.markdown("### 📊 Dataset Overview")
    st.write(
        "Calculated on a curated dataset of **118 NIK inhibitors** collected from published literature, "
        "standardized to experimental pIC50 values."
    )

elif nav == "Model Performance":
    st.markdown("### 📈 Model Performance Overview")
    st.write(
        "Four Monte Carlo 2D-QSAR models (M1–M4) were developed using the CORALSEA correlation balance approach. "
        "**Model M3** was selected as the final model based on its overall statistical performance and lower prediction error."
    )
    col_m1, col_m2, col_m3 = st.columns(3)
    col_m1.metric("R² Validation", "0.765")
    col_m2.metric("Q² Validation", "0.660")
    col_m3.metric("MAE Validation", "0.291")

elif nav == "Collaboration & Contact":
    st.markdown("### 🤝 Collaboration & Contact")
    st.write(
        "We are open to scientific collaboration, co-development of tools, and data exchange projects in computational drug discovery, "
        "molecular docking, dynamic and QSAR modeling."
    )
    st.markdown("📧 **Laboratory Contact:** `drugdesignsynthesislab@gmail.com`")
    st.markdown("📧 **Faculty Contact:** `Yogita_pharma@pbi.ac.in`")
    st.markdown("---")
    st.markdown("### 📖 How to Cite")
    st.info("If you use NIK-Predict in your research, please cite us: `[Citation link will be provided after publication]`")

elif nav == "Limitations":
    st.markdown("### ⚠️ Model Limitations")
    st.warning("• **Small Molecules:** Applicable strictly to small-molecule inhibitors. Not tested for macrocycles, peptides, or prodrugs.")
    st.warning("• **Chemical Space:** Predictions may be less reliable outside the model's chemical space domain.")
    st.warning("• **Experimental Validation:** In-silico predictions provide preliminary bioactivity estimates; experimental validation is required.")
