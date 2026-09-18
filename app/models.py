from typing import List, Dict, Optional, Any
from pydantic import BaseModel, Field

class PredictSingleRequest(BaseModel):
    smiles: str = Field(..., description="Canonical or valid SMILES string of the candidate molecule", example="COc1cnc(nc1N1CCc2c1cc(Br)cc2)N")

class PredictBatchRequest(BaseModel):
    smiles_list: List[str] = Field(..., description="List of SMILES strings for batch evaluation")

class RunPrediction(BaseModel):
    run: str = Field(..., description="Run identifier (Run-1, Run-2, Run-3)")
    endpoint: float = Field(..., description="Predicted biological inhibitory activity (pIC50)")
    smiles_weight: float = Field(..., description="Descriptor of Correlation Weights (DCW) / SMILES weight")
    defect_smiles: float = Field(..., description="Defect SMILES value (SAk rare attribute penalty)")

class PhysicochemicalProperties(BaseModel):
    formula: str = Field(..., description="Molecular formula")
    molecular_weight: float = Field(..., description="Molecular weight in g/mol")
    logp: float = Field(..., description="Lipophilicity (LogP)")
    tpsa: float = Field(..., description="Topological Polar Surface Area in Å²")
    h_donors_acceptors: str = Field(..., description="H-Bond Donors / Acceptors")
    rotatable_bonds: int = Field(..., description="Number of rotatable bonds")

class SinglePredictionResult(BaseModel):
    smiles: str = Field(..., description="Evaluated SMILES string")
    target: str = Field("NF-κB Inducing Kinase (NIK / MAP3K14)", description="Target kinase")
    consensus_prediction: float = Field(..., description="Consensus average predicted pIC50 across all 3 runs")
    consensus_smiles_weight: float = Field(..., description="Average SMILES weight (DCW) across all 3 runs")
    consensus_defect_smiles: float = Field(..., description="Average Defect SMILES across all 3 runs")
    runs: Dict[str, RunPrediction] = Field(..., description="Breakdown of predictions, weights, and defect per individual run")
    elapsed_seconds: float = Field(..., description="Prediction execution duration in seconds")
    physicochemical_properties: Optional[PhysicochemicalProperties] = None
    svg_structure: Optional[str] = None

class PredictResponse(BaseModel):
    target: str = "NF-κB Inducing Kinase (NIK / MAP3K14)"
    results: List[SinglePredictionResult]
    total_elapsed_seconds: float
    department: str = "Department of Pharmaceutical Sciences and Drug Research"
    lab_affiliation: str = "Drug Design Synthesis Lab, Punjabi University, Patiala"
    contact_email: str = "drugdesignsynthesislab@gmail.com"

class SampleCompound(BaseModel):
    id: str
    name: str
    smiles: str
    experimental_pIC50: Optional[float] = None
    description: str

class ModelInfoResponse(BaseModel):
    model_name: str
    target: str
    method: str
    department: str
    affiliation: str
    contact_email: str
    citation: str
    dataset_summary: Dict[str, Any]
    statistical_validation: Dict[str, Any]
    applicability_domain: Dict[str, Any]
    limitations: List[str]
    runs: Dict[str, Dict[str, Any]]
    samples: List[SampleCompound]
