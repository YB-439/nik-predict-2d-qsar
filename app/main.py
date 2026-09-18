import os
from typing import List
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles
from fastapi.responses import FileResponse

from app.models import (
    PredictSingleRequest,
    PredictBatchRequest,
    PredictResponse,
    SampleCompound,
    ModelInfoResponse,
)
from app.predictor import nik_predictor_instance, NIK_MODEL_CONFIG

app = FastAPI(
    title="NF-κB Inducing Kinase (NIK) 2D-QSAR Portal",
    description=(
        "Computational regression portal for predicting NIK inhibitory potency (pIC50) "
        "using CORALSEA Monte Carlo 3-run consensus modeling. Developed at Drug Design Synthesis Lab, "
        "Department of Pharmaceutical Sciences and Drug Research, Punjabi University, Patiala."
    ),
    version="1.0.0",
)

# Enable CORS for open-source integration
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Path to static folder
STATIC_DIR = os.path.join(os.path.dirname(__file__), "static")
if os.path.isdir(STATIC_DIR):
    app.mount("/static", StaticFiles(directory=STATIC_DIR), name="static")

NIK_SAMPLES = [
    {
        "id": "NIK001",
        "name": "Compound 1",
        "smiles": "COc1cnc(nc1N1CCc2c1cc(Br)cc2)N",
        "experimental_pIC50": 5.0693,
        "description": "Potent aminopyrimidine NIK inhibitor derivative"
    },
    {
        "id": "NIK015",
        "name": "Lead Compound NIK015",
        "smiles": "N#Cc1ccc(cc1)c1cnc(s1)C(=O)Nc1nccc(n1)n1cnc2c1ccc(c2)Cl",
        "experimental_pIC50": 6.2254,
        "description": "High affinity nitrile-substituted benzimidazole derivative"
    },
    {
        "id": "NIK009",
        "name": "Benchmark Compound NIK009",
        "smiles": "Nc1nc(N2CCc3c2cc(OC)cc3)c(Cl)cn1",
        "experimental_pIC50": 7.0500,
        "description": "Methoxy-substituted chloropyrimidine analog"
    },
    {
        "id": "NIK017",
        "name": "Heteroaryl Analog NIK017",
        "smiles": "Nc1nc(N2CCc3c2cc(c2n[nH]cc2)cc3)c(Cl)cn1",
        "experimental_pIC50": 8.1900,
        "description": "Pyrazole-substituted high-potency analog"
    }
]

@app.get("/")
async def serve_portal():
    """Serves the main web UI."""
    index_file = os.path.join(STATIC_DIR, "index.html")
    if os.path.isfile(index_file):
        return FileResponse(index_file)
    return {"message": "NF-κB Inducing Kinase (NIK) 2D-QSAR Portal API running. UI file not found."}

@app.get("/health")
async def health_check():
    """Health check endpoint."""
    return {
        "status": "ok",
        "target": "NF-κB Inducing Kinase (NIK)",
        "department": "Department of Pharmaceutical Sciences and Drug Research",
        "institution": "Punjabi University, Patiala",
    }

@app.get("/api/info", response_model=ModelInfoResponse)
async def get_model_info():
    """Returns model metadata, author attribution, lab details, and statistical specifications."""
    return ModelInfoResponse(
        model_name="NF-κB Inducing Kinase (NIK) 2D-QSAR CORALSEA Regression Model",
        target="NF-κB Inducing Kinase (NIK / MAP3K14)",
        method="CORALSEA-2024 Monte Carlo Optimization with Balance of Correlation (HFG EC0 + SMILES attributes)",
        department="Department of Pharmaceutical Sciences and Drug Research",
        affiliation="Drug Design Synthesis Lab, Punjabi University, Patiala",
        contact_email="drugdesignsynthesislab@gmail.com, Yogita_pharma@pbi.ac.in",
        citation="Model developed at Drug Design Synthesis Lab, Department of Pharmaceutical Sciences and Drug Research, Punjabi University, Patiala. Faculty Lead: Prof. Yogita Bansal.",
        dataset_summary={
            "total_compounds": 118,
            "experimental_assay": "Luminescent biochemical kinase assay (cell-free IC50 standardized to molar pIC50)",
            "data_splitting": "Las Vegas algorithm (Split-3 selected)",
            "percentage_identity": "≤ 35% overlap across splits ensuring independent validation"
        },
        statistical_validation={
            "selected_model": "Model M3 (Split-3)",
            "r2_validation": 0.765,
            "q2_validation": 0.660,
            "ccc_validation": 0.861,
            "mae_validation": 0.291,
            "rm2_validation": 0.7746,
            "delta_rm2": 0.0849,
            "y_randomization": "10-iteration randomization confirmed model validity (R²_rand << 0.05, no chance correlation)",
            "consensus_scheme": "Exact arithmetic mean across Run-1, Run-2, and Run-3"
        },
        applicability_domain={
            "method": "Defect(SMILES) rare structural attribute penalty",
            "training_mean_defect": 3.99,
            "threshold_defect": 7.98,
            "threshold_rule": "Defect(SMILES) < 7.98 (2 × training mean)"
        },
        limitations=[
            "Applicability Domain: Reliable for molecules within the chemical space of the training set with Defect(SMILES) < 7.98.",
            "Biological Scope: Predicts cell-free biochemical kinase inhibitory potency (pIC50); cellular permeability and pharmacokinetics require experimental confirmation."
        ],
        runs=NIK_MODEL_CONFIG,
        samples=NIK_SAMPLES,
    )

@app.get("/api/samples", response_model=List[SampleCompound])
async def get_sample_compounds():
    """Provides curated sample SMILES compounds from the published dataset for 1-click testing."""
    return NIK_SAMPLES

@app.post("/api/predict", response_model=PredictResponse)
async def predict(request: PredictSingleRequest):
    """
    Main prediction endpoint for a single candidate SMILES string.
    Evaluates across all 3 Monte Carlo runs and returns consensus pIC50, DCW, and Defect(SMILES).
    """
    smiles = request.smiles.strip()
    if not smiles:
        raise HTTPException(status_code=400, detail="SMILES string cannot be empty.")

    try:
        raw_res = nik_predictor_instance.predict_batch([smiles])
        return PredictResponse(
            results=raw_res["results"],
            total_elapsed_seconds=raw_res["total_elapsed_seconds"],
        )
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Prediction error: {str(e)}")

@app.post("/api/predict/batch", response_model=PredictResponse)
async def predict_batch(request: PredictBatchRequest):
    """
    Batch prediction endpoint for multiple SMILES strings.
    Evaluates all SMILES in parallel across Run-1, Run-2, and Run-3.
    """
    if not request.smiles_list:
        raise HTTPException(status_code=400, detail="SMILES list cannot be empty.")

    clean_smiles = [s.strip() for s in request.smiles_list if s.strip()]
    if not clean_smiles:
        raise HTTPException(status_code=400, detail="No non-empty SMILES provided.")

    try:
        raw_res = nik_predictor_instance.predict_batch(clean_smiles)
        return PredictResponse(
            results=raw_res["results"],
            total_elapsed_seconds=raw_res["total_elapsed_seconds"],
        )
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Batch prediction error: {str(e)}")
