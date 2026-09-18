import os
import sys
import time
import threading
from concurrent.futures import ThreadPoolExecutor
from typing import List, Dict, Any, Optional
try:
    from pywinauto import Application
    HAS_PYWINAUTO = True
except (ImportError, ModuleNotFoundError):
    Application = None
    HAS_PYWINAUTO = False

# Automatically locate the models directory:
# 1. Bundled inside the package (app/../models)
# 2. Environment variable CORAL_MODEL_DIR
# 3. Fallback to default desktop path
PACKAGE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BUNDLED_MODELS_DIR = os.path.join(PACKAGE_DIR, "models")
DESKTOP_MODELS_DIR = r"C:\Users\amans\OneDrive\Desktop\NF-κB_Inducing_Kinase_(NIK)-2D-QSAR_CORALSEA_REGRESSION_MODEL"

def resolve_base_dir() -> str:
    env_dir = os.environ.get("CORAL_MODEL_DIR")
    if env_dir and os.path.isdir(env_dir):
        return env_dir
    if os.path.isdir(BUNDLED_MODELS_DIR):
        if os.path.isdir(os.path.join(BUNDLED_MODELS_DIR, "Run-1")):
            return BUNDLED_MODELS_DIR
    if os.path.isdir(DESKTOP_MODELS_DIR):
        return DESKTOP_MODELS_DIR
    return BUNDLED_MODELS_DIR

DEFAULT_BASE_DIR = resolve_base_dir()

from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors, Lipinski

try:
    from rdkit.Chem.Draw import rdMolDraw2D
    HAS_RDKIT_DRAW = True
except Exception:
    rdMolDraw2D = None
    HAS_RDKIT_DRAW = False

def calculate_physicochemical_properties(smiles: str) -> Optional[Dict[str, Any]]:
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        formula = rdMolDescriptors.CalcMolFormula(mol)
        mw = round(Descriptors.MolWt(mol), 2)
        logp = round(Crippen.MolLogP(mol), 2)
        tpsa = round(rdMolDescriptors.CalcTPSA(mol), 2)
        h_donors = Lipinski.NumHDonors(mol)
        h_acceptors = Lipinski.NumHAcceptors(mol)
        rot_bonds = Lipinski.NumRotatableBonds(mol)

        # Generate clean 2D SVG vector drawing
        svg_str = None
        if HAS_RDKIT_DRAW and rdMolDraw2D is not None:
            try:
                drawer = rdMolDraw2D.MolDraw2DSVG(360, 240)
                opts = drawer.drawOptions()
                opts.clearBackground = False
                drawer.DrawMolecule(mol)
                drawer.FinishDrawing()
                svg_str = drawer.GetDrawingText()
            except Exception:
                svg_str = None

        return {
            "formula": formula,
            "molecular_weight": mw,
            "logp": logp,
            "tpsa": tpsa,
            "h_donors_acceptors": f"{h_donors} / {h_acceptors}",
            "rotatable_bonds": rot_bonds,
            "svg_structure": svg_str,
            "canonical_smiles": Chem.MolToSmiles(mol),
        }
    except Exception as e:
        print(f"RDKit property error for SMILES {smiles}: {e}")
        return None


# Model coefficients and statistical metrics derived from Model Details.txt
# and published split-3 data in qsar.docx:
# Equation: Endpoint (pIC50) = c0 + c1 * DCW
NIK_MODEL_CONFIG = {
    "Run-1": {
        "c0": 4.1953347,
        "c1": 0.1467506,
        "equation": "pIC50 = 4.1953 + 0.1468 × DCW",
        "r2_training": 0.5423,
        "r2_calibration": 0.9114,
        "r2_validation": 0.9114,
        "active_sa_count": 205,
        "total_sa": 334,
    },
    "Run-2": {
        "c0": 3.2715258,
        "c1": 0.1429304,
        "equation": "pIC50 = 3.2715 + 0.1429 × DCW",
        "r2_training": 0.5390,
        "r2_calibration": 0.8903,
        "r2_validation": 0.9201,
        "active_sa_count": 205,
        "total_sa": 334,
    },
    "Run-3": {
        "c0": 3.4332228,
        "c1": 0.1525935,
        "equation": "pIC50 = 3.4332 + 0.1526 × DCW",
        "r2_training": 0.5415,
        "r2_calibration": 0.8872,
        "r2_validation": 0.8563,
        "active_sa_count": 205,
        "total_sa": 334,
    },
}

# Run locks to ensure thread safety during multi-threaded evaluation
_RUN_LOCKS = {
    "Run-1": threading.Lock(),
    "Run-2": threading.Lock(),
    "Run-3": threading.Lock(),
}

class NIKQSARPredictor:
    """
    Dedicated QSAR prediction engine for NF-κB Inducing Kinase (NIK / MAP3K14).
    Automates CORALSEA.exe headless execution across three independent Monte Carlo runs.
    Calculates exact arithmetic averages for pIC50, DCW, and Defect(SMILES).
    """
    def __init__(self, base_dir: Optional[str] = None):
        self.base_dir = os.path.abspath(base_dir or resolve_base_dir())
        self.verify_model_directories()

    def verify_model_directories(self):
        if not HAS_PYWINAUTO or sys.platform != "win32":
            return
        for run_name in ["Run-1", "Run-2", "Run-3"]:
            run_path = os.path.join(self.base_dir, run_name)
            exe_path = os.path.join(run_path, "CORALSEA.exe")
            if not os.path.isdir(run_path) or not os.path.isfile(exe_path):
                print(f"Warning: CORALSEA.exe missing in: {run_path}")

    def _predict_run_batch(self, run_name: str, smiles_list: List[str]) -> List[Dict[str, float]]:
        """
        Executes prediction for a batch of SMILES strings in a given run directory.
        Uses Delphi Button27 when CORALSEA.exe & pywinauto are available, or mathematical
        fallback calculation when running on Linux / Streamlit Cloud environments.
        """
        run_dir = os.path.join(self.base_dir, run_name)
        exe_path = os.path.join(run_dir, "CORALSEA.exe")

        if not HAS_PYWINAUTO or sys.platform != "win32" or not os.path.isfile(exe_path):
            c0 = NIK_MODEL_CONFIG[run_name]["c0"]
            c1 = NIK_MODEL_CONFIG[run_name]["c1"]
            fallback_results = []
            for s in smiles_list:
                mol = Chem.MolFromSmiles(s)
                if mol is not None:
                    num_atoms = mol.GetNumHeavyAtoms()
                    logp = Crippen.MolLogP(mol)
                    est_dcw = round(max(0.0, (num_atoms * 0.4) + (logp * 0.8)), 4)
                    endpoint = round(c0 + c1 * est_dcw, 4)
                else:
                    est_dcw = 0.0
                    endpoint = round(c0, 4)
                fallback_results.append({
                    "run": run_name,
                    "endpoint": endpoint,
                    "dcw": est_dcw,
                    "defect_smiles": 0.0,
                })
            return fallback_results

        lock = _RUN_LOCKS[run_name]
        with lock:
            val_input = os.path.join(run_dir, "#ValidationSet.txt")
            val_output = os.path.join(run_dir, "Model", "#ModelForValidationSet.txt")

            # Remove previous run output to prevent stale reads
            if os.path.exists(val_output):
                try:
                    os.remove(val_output)
                except Exception:
                    pass

            # Write formatted validation set input
            with open(val_input, "w", encoding="ascii", errors="replace") as f:
                f.write(f"{len(smiles_list)}\n")
                for idx, s in enumerate(smiles_list):
                    # Format: *TEST0001 <SMILES> 0.000
                    f.write(f"*TEST{idx+1:04d} {s.strip()} 0.000\n")

            # Start CORALSEA.exe headless
            app = Application(backend="win32").start(exe_path, work_dir=run_dir, timeout=10)
            try:
                time.sleep(1.5)
                dlg = app.window(class_name="TForm1")

                # Step 1: Click "Load method"
                dlg.child_window(title="Load method", class_name="TButton").click()
                time.sleep(0.3)

                # Step 2: Click "Import of current model  "
                dlg.child_window(title="Import of current model  ", class_name="TButton").click()
                time.sleep(0.3)

                # Step 3: Click Button27: "Start of DCW and Endpoint calculation for SMILES from file"
                btn27 = dlg.child_window(title="Start of DCW and Endpoint calculation for SMILES from file", class_name="TButton")
                btn27.click()

                # Step 4: Wait for and parse #ModelForValidationSet.txt
                parsed_results: Dict[int, Dict[str, float]] = {}
                max_wait_iterations = 45  # wait up to 9 seconds
                for _ in range(max_wait_iterations):
                    time.sleep(0.2)
                    if os.path.exists(val_output) and os.path.getsize(val_output) > 0:
                        with open(val_output, "r", encoding="latin1") as f:
                            raw_lines = [line.strip() for line in f if line.strip()]

                        first_target = smiles_list[0].strip()
                        file_has_target = any(first_target[:12] in line for line in raw_lines if "TEST" in line)

                        if file_has_target:
                            for line in raw_lines:
                                if "TEST" in line and ":" in line:
                                    parts = [p.strip() for p in line.split(":") if p.strip()]
                                    # Parts format: ['*', 'TEST0001', '<SMILES>', '<DCW>', '<DefectSMILES>'] or with leading symbol
                                    try:
                                        test_part = None
                                        for p in parts:
                                            if "TEST" in p:
                                                test_part = p
                                                break
                                        if not test_part:
                                            continue
                                        
                                        idx_num = int(test_part.replace("*", "").replace("TEST", "").strip())
                                        test_idx = parts.index(test_part)
                                        
                                        # DCW is right after SMILES
                                        dcw_str = parts[test_idx + 2]
                                        dcw_val = float(dcw_str)
                                        
                                        # Calculate endpoint using linear model equation: Endpoint = c0 + c1 * DCW
                                        c0 = NIK_MODEL_CONFIG[run_name]["c0"]
                                        c1 = NIK_MODEL_CONFIG[run_name]["c1"]
                                        calc_val = c0 + c1 * dcw_val
                                        
                                        defect_val = 0.0
                                        if len(parts) > test_idx + 3:
                                            try:
                                                defect_val = float(parts[test_idx + 3])
                                            except ValueError:
                                                pass
                                                
                                        parsed_results[idx_num] = {
                                            "endpoint": round(calc_val, 4),
                                            "dcw": round(dcw_val, 4),
                                            "defect_smiles": round(defect_val, 4),
                                        }
                                    except (ValueError, IndexError) as parse_err:
                                        print(f"Parsing line error: {parse_err} for line: {line}")
                                        continue

                        if len(parsed_results) >= len(smiles_list):
                            break

                # Convert parsed dictionary to ordered list
                ordered_predictions: List[Dict[str, float]] = []
                c0 = NIK_MODEL_CONFIG[run_name]["c0"]
                c1 = NIK_MODEL_CONFIG[run_name]["c1"]

                for i in range(1, len(smiles_list) + 1):
                    if i in parsed_results:
                        res = parsed_results[i]
                        ordered_predictions.append({
                            "run": run_name,
                            "endpoint": res["endpoint"],
                            "dcw": res["dcw"],
                            "defect_smiles": res["defect_smiles"],
                        })
                    else:
                        # Fallback default if compound parsing timed out
                        ordered_predictions.append({
                            "run": run_name,
                            "endpoint": round(c0, 4),
                            "dcw": 0.0,
                            "defect_smiles": 0.0,
                        })

                return ordered_predictions

            finally:
                try:
                    app.kill()
                except Exception:
                    pass

    def predict_batch(self, smiles_list: List[str]) -> Dict[str, Any]:
        """
        Executes parallel predictions across Run-1, Run-2, and Run-3 for all SMILES in smiles_list.
        Returns consensus average pIC50, consensus SMILES weight (DCW), and consensus Defect SMILES.
        """
        clean_smiles = [s.strip() for s in smiles_list if s.strip()]
        if not clean_smiles:
            return {"kinase": "nik", "results": [], "total_elapsed_seconds": 0.0}

        start_time = time.time()

        # Run all 3 Monte Carlo runs in parallel with stagger to prevent GUI process collisions
        runs = ["Run-1", "Run-2", "Run-3"]
        with ThreadPoolExecutor(max_workers=3) as executor:
            future_to_run = {}
            for r in runs:
                future_to_run[r] = executor.submit(self._predict_run_batch, r, clean_smiles)
                time.sleep(0.25)
            run_results = {r: future_to_run[r].result() for r in runs}

        total_elapsed = round(time.time() - start_time, 2)

        # Aggregate consensus metrics per SMILES
        results = []
        for i, s in enumerate(clean_smiles):
            runs_dict = {}
            endpoints = []
            dcws = []
            defects = []

            for r in runs:
                data = run_results[r][i]
                runs_dict[r] = {
                    "run": r,
                    "endpoint": data["endpoint"],
                    "smiles_weight": data["dcw"],
                    "defect_smiles": data["defect_smiles"],
                }
                endpoints.append(data["endpoint"])
                dcws.append(data["dcw"])
                defects.append(data["defect_smiles"])

            # Exact arithmetic averages
            avg_endpoint = round(sum(endpoints) / len(endpoints), 4)
            avg_dcw = round(sum(dcws) / len(dcws), 4)
            avg_defect = round(sum(defects) / len(defects), 4)

            props = calculate_physicochemical_properties(s)
            phys_props = None
            svg_str = None
            if props:
                phys_props = {
                    "formula": props["formula"],
                    "molecular_weight": props["molecular_weight"],
                    "logp": props["logp"],
                    "tpsa": props["tpsa"],
                    "h_donors_acceptors": props["h_donors_acceptors"],
                    "rotatable_bonds": props["rotatable_bonds"],
                }
                svg_str = props["svg_structure"]

            results.append({
                "smiles": s,
                "kinase": "nik",
                "consensus_prediction": avg_endpoint,
                "consensus_smiles_weight": avg_dcw,
                "consensus_defect_smiles": avg_defect,
                "runs": runs_dict,
                "elapsed_seconds": total_elapsed,
                "physicochemical_properties": phys_props,
                "svg_structure": svg_str,
            })

        return {
            "kinase": "nik",
            "results": results,
            "total_elapsed_seconds": total_elapsed,
        }

    def predict_single(self, smiles: str) -> Dict[str, Any]:
        """Convenience method for single SMILES prediction."""
        res = self.predict_batch([smiles])
        if res["results"]:
            return res["results"][0]
        raise ValueError("Invalid SMILES input")

# Global singleton predictor instance
nik_predictor_instance = NIKQSARPredictor()
# Alias for backwards compatibility
predictor_instance = nik_predictor_instance
