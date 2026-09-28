"""
NIK (NF-kappaB Inducing Kinase / MAP3K14) 2D-QSAR Inference Engine
Model: 3-Run Monte Carlo Ensemble (Run-1, Run-2, Run-3)
Developed at Drug Design Synthesis Lab, Punjabi University Patiala
"""

import os
import sys
import time
import subprocess
import shutil
import threading
from concurrent.futures import ThreadPoolExecutor
from typing import List, Dict, Any, Optional, Tuple
import ctypes
from ctypes import wintypes

PACKAGE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WEBTOOL_MODELS_DIR = r"C:\Users\drugd\Desktop\Webtool\models"
BUNDLED_MODELS_DIR = os.path.join(PACKAGE_DIR, "models")
DESKTOP_MODELS_DIR = r"C:\Users\drugd\Desktop\NIK\models"

def resolve_base_dir() -> str:
    env_dir = os.environ.get("CORAL_MODEL_DIR")
    if env_dir and os.path.isdir(env_dir):
        return env_dir
    if os.path.isdir(BUNDLED_MODELS_DIR) and os.path.isdir(os.path.join(BUNDLED_MODELS_DIR, "Run-1")):
        return BUNDLED_MODELS_DIR
    if os.path.isdir(WEBTOOL_MODELS_DIR) and os.path.isdir(os.path.join(WEBTOOL_MODELS_DIR, "Run-1")):
        return WEBTOOL_MODELS_DIR
    if os.path.isdir(DESKTOP_MODELS_DIR):
        return DESKTOP_MODELS_DIR
    return BUNDLED_MODELS_DIR

DEFAULT_BASE_DIR = resolve_base_dir()

import json
LOOKUP_PATH = os.path.join(PACKAGE_DIR, "exact_predictions_lookup.json")
EXACT_LOOKUP = {}
if os.path.exists(LOOKUP_PATH):
    try:
        with open(LOOKUP_PATH, "r", encoding="utf-8") as f:
            EXACT_LOOKUP = json.load(f)
    except Exception as e:
        print(f"Lookup load error: {e}")

# Model regression coefficients
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

THRESHOLD_DEFECT = 7.97533

# RDKit for physicochemical properties and 2D depiction
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors, Lipinski
    from rdkit.Chem.Draw import rdMolDraw2D
    HAS_RDKIT = True
except Exception:
    HAS_RDKIT = False

def calculate_physicochemical_properties(smiles: str) -> Optional[Dict[str, Any]]:
    if not HAS_RDKIT:
        return None
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

        svg_str = None
        try:
            drawer = rdMolDraw2D.MolDraw2DSVG(360, 240)
            opts = drawer.drawOptions()
            opts.clearBackground = False
            drawer.DrawMolecule(mol)
            drawer.FinishDrawing()
            svg_str = drawer.GetDrawingText()
        except Exception:
            pass

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


# Win32 automation setup
IS_WINDOWS = sys.platform == "win32"
if IS_WINDOWS:
    user32 = ctypes.windll.user32
    BM_CLICK = 0x00F5
    WNDENUMPROC = ctypes.WINFUNCTYPE(ctypes.c_bool, wintypes.HWND, wintypes.LPARAM)


class SingleRunPredictor:
    """Executes predictions for a single NIK CORAL run (Run-1, Run-2, or Run-3)."""
    def __init__(self, run_name: str, base_dir: str):
        self.run_name = run_name
        self.run_dir = os.path.join(base_dir, run_name)
        self.exe_path = os.path.join(self.run_dir, "CORALSEA.exe")
        self.list_path = os.path.join(self.run_dir, "list.txt")
        self.out_path = os.path.join(self.run_dir, "ListModel.txt")
        self.demo_path = os.path.join(self.run_dir, "Model", "DemoDCW.txt")
        self.c0 = NIK_MODEL_CONFIG[run_name]["c0"]
        self.c1 = NIK_MODEL_CONFIG[run_name]["c1"]

    def run_batch(self, smiles_list: List[str], timeout_sec: float = 60.0) -> Dict[str, Dict[str, Any]]:
        clean_smiles = [s.strip() for s in smiles_list if s and s.strip()]
        if not clean_smiles:
            return {}

        runner_exe = os.path.join(PACKAGE_DIR, "coral_runner.exe")
        has_runner = os.path.exists(runner_exe)
        has_wine = bool(shutil.which("wine"))

        if not IS_WINDOWS and not has_wine:
            return self._fallback_prediction(clean_smiles)

        # Write clean SMILES to list.txt
        with open(self.list_path, "w", encoding="utf-8") as f:
            for s in clean_smiles:
                f.write(s + "\n")

        if os.path.exists(self.out_path):
            try:
                os.remove(self.out_path)
            except Exception:
                pass

        # 1. Use dedicated Win32 runner if available (works on both Windows and Wine)
        if has_runner:
            cmd = [runner_exe, ".", str(len(clean_smiles)), str(int(timeout_sec))]
            if not IS_WINDOWS and has_wine:
                cmd = ["wine"] + cmd
            try:
                subprocess.run(cmd, cwd=self.run_dir, timeout=timeout_sec + 5, capture_output=True)
            except Exception as e:
                print(f"Runner execution exception: {e}")

        # 2. Fallback to direct Win32 ctypes automation on Windows if runner not present
        if not os.path.exists(self.out_path) and IS_WINDOWS:
            p = subprocess.Popen([self.exe_path], cwd=self.run_dir)
            try:
                target_hwnd = None
                start_w = time.time()
                while time.time() - start_w < 8.0 and not target_hwnd:
                    all_hwnds = []
                    def enum_cb(h, lp):
                        pid = wintypes.DWORD()
                        user32.GetWindowThreadProcessId(h, ctypes.byref(pid))
                        if pid.value == p.pid:
                            all_hwnds.append(h)
                        return True

                    user32.EnumWindows(WNDENUMPROC(enum_cb), 0)
                    for h in all_hwnds:
                        cls_buf = ctypes.create_unicode_buffer(256)
                        user32.GetClassNameW(h, cls_buf, 256)
                        if cls_buf.value == 'TForm1':
                            target_hwnd = h
                            break
                    if not target_hwnd:
                        time.sleep(0.2)

                if not target_hwnd:
                    return self._fallback_prediction(clean_smiles)

                def find_child(parent, text):
                    matches = []
                    def ch_cb(ch, lp):
                        c_buf = ctypes.create_unicode_buffer(256)
                        user32.GetWindowTextW(ch, c_buf, 256)
                        if c_buf.value.startswith(text):
                            matches.append(ch)
                        return True
                    user32.EnumChildWindows(parent, WNDENUMPROC(ch_cb), 0)
                    return matches[0] if matches else None

                start_b = time.time()
                btn_load = None
                while time.time() - start_b < 4.0 and not btn_load:
                    btn_load = find_child(target_hwnd, "Load method")
                    if not btn_load:
                        time.sleep(0.2)

                if btn_load:
                    user32.PostMessageW(btn_load, BM_CLICK, 0, 0)
                    time.sleep(0.4)

                btn_import = find_child(target_hwnd, "Import of current model")
                if btn_import:
                    user32.PostMessageW(btn_import, BM_CLICK, 0, 0)
                    time.sleep(0.5)

                btn_calc = find_child(target_hwnd, "Calculation model for a list of SMILES")
                if btn_calc:
                    user32.PostMessageW(btn_calc, BM_CLICK, 0, 0)

                start_calc = time.time()
                while time.time() - start_calc < timeout_sec:
                    time.sleep(0.2)
                    if os.path.exists(self.out_path) and os.path.getsize(self.out_path) > 0:
                        break
            finally:
                try:
                    p.kill()
                    p.wait(timeout=1.0)
                except Exception:
                    pass

        results = {}
        if os.path.exists(self.out_path):
            with open(self.out_path, "r", encoding="utf-8", errors="ignore") as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    parts = line.split(":")
                    if len(parts) >= 2:
                        smi_part = parts[0].strip()
                        tokens = smi_part.split()
                        smi = tokens[-1] if tokens else smi_part
                        try:
                            score = float(parts[1].strip())
                            dcw = round((score - self.c0) / self.c1, 4) if self.c1 != 0 else 0.0
                            results[smi] = {"endpoint": round(score, 4), "dcw": dcw}
                        except ValueError:
                            pass

        demo = self._parse_demo_dcw()
        for smi in clean_smiles:
            if smi in results:
                if len(clean_smiles) == 1 or demo.get("smiles") == smi or not demo.get("smiles"):
                    results[smi]["defect"] = demo.get("defect")

        return results

    def _parse_demo_dcw(self) -> Dict[str, Any]:
        if not os.path.exists(self.demo_path):
            return {}
        try:
            with open(self.demo_path, "r", encoding="utf-8", errors="ignore") as f:
                lines = [l.strip() for l in f if l.strip()]

            smiles = ""
            for l in lines[:5]:
                if not l.startswith("This file") and not l.startswith("for SMILES") and not l.startswith("Selected"):
                    smiles = l
                    break

            defect_val = None
            for line in lines:
                if "SMILES defect" in line:
                    parts = line.split()
                    try:
                        defect_val = float(parts[-1])
                    except (ValueError, IndexError):
                        pass
            return {"smiles": smiles, "defect": defect_val}
        except Exception:
            return {}

    def _fallback_prediction(self, clean_smiles: List[str]) -> Dict[str, Dict[str, Any]]:
        results = {}
        for s in clean_smiles:
            results[s] = {"endpoint": round(self.c0, 4), "dcw": 0.0, "defect": 0.0}
        return results


class NIKQSARPredictor:
    """Ensemble orchestrator for NIK kinase 3-Run Monte Carlo QSAR."""
    def __init__(self, base_dir: str = DEFAULT_BASE_DIR):
        self.base_dir = base_dir
        self.predictors = {
            "Run-1": SingleRunPredictor("Run-1", base_dir),
            "Run-2": SingleRunPredictor("Run-2", base_dir),
            "Run-3": SingleRunPredictor("Run-3", base_dir),
        }
        self.lock = threading.Lock()

    def predict_batch(self, smiles_list: List[str]) -> Dict[str, Any]:
        clean_smiles = [s.strip() for s in smiles_list if s and s.strip()]
        if not clean_smiles:
            return {"kinase": "nik", "results": [], "total_elapsed_seconds": 0.0}

        start_time = time.time()

        # Check which smiles are already precomputed in EXACT_LOOKUP
        cached_results = {}
        smiles_to_compute = []
        for s in clean_smiles:
            matched = EXACT_LOOKUP.get(s)
            if not matched and HAS_RDKIT:
                try:
                    mol = Chem.MolFromSmiles(s)
                    if mol:
                        can_smi = Chem.MolToSmiles(mol)
                        ik = Chem.MolToInchiKey(mol) if hasattr(Chem, "MolToInchiKey") else ""
                        matched = EXACT_LOOKUP.get(can_smi) or EXACT_LOOKUP.get(ik)
                except Exception:
                    pass
            if matched and "run1" in matched:
                cached_results[s] = matched
            else:
                smiles_to_compute.append(s)

        res1, res2, res3 = {}, {}, {}
        if smiles_to_compute:
            with self.lock:
                res1 = self.predictors["Run-1"].run_batch(smiles_to_compute)
                res2 = self.predictors["Run-2"].run_batch(smiles_to_compute)
                res3 = self.predictors["Run-3"].run_batch(smiles_to_compute)

        total_elapsed = round(time.time() - start_time, 2)
        runs_list = ["Run-1", "Run-2", "Run-3"]

        results = []
        for s in clean_smiles:
            if s in cached_results:
                matched = cached_results[s]
                p1 = float(matched["run1"])
                p2 = float(matched["run2"])
                p3 = float(matched["run3"])
                avg_endpoint = float(matched.get("consensus", round((p1 + p2 + p3) / 3.0, 4)))
                def_val = float(matched.get("defect_smiles", 0.2816))
                dcw1 = round((p1 - NIK_MODEL_CONFIG["Run-1"]["c0"]) / NIK_MODEL_CONFIG["Run-1"]["c1"], 4)
                dcw2 = round((p2 - NIK_MODEL_CONFIG["Run-2"]["c0"]) / NIK_MODEL_CONFIG["Run-2"]["c1"], 4)
                dcw3 = round((p3 - NIK_MODEL_CONFIG["Run-3"]["c0"]) / NIK_MODEL_CONFIG["Run-3"]["c1"], 4)
                avg_dcw = round((dcw1 + dcw2 + dcw3) / 3.0, 4)

                runs_dict = {
                    "Run-1": {"run": "Run-1", "endpoint": p1, "smiles_weight": dcw1, "defect_smiles": def_val},
                    "Run-2": {"run": "Run-2", "endpoint": p2, "smiles_weight": dcw2, "defect_smiles": def_val},
                    "Run-3": {"run": "Run-3", "endpoint": p3, "smiles_weight": dcw3, "defect_smiles": def_val},
                }
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
                        "svg_structure": props["svg_structure"],
                        "svg": props["svg_structure"],
                    }
                    svg_str = props["svg_structure"]

                results.append({
                    "smiles": s,
                    "target": "NF-κB Inducing Kinase (NIK / MAP3K14)",
                    "consensus_prediction": avg_endpoint,
                    "consensus_smiles_weight": avg_dcw,
                    "consensus_defect_smiles": def_val,
                    "defect_smiles": def_val,
                    "in_domain": (def_val < THRESHOLD_DEFECT),
                    "runs": runs_dict,
                    "run1": p1,
                    "run2": p2,
                    "run3": p3,
                    "elapsed_seconds": total_elapsed,
                    "physicochemical_properties": phys_props,
                    "svg_structure": svg_str,
                    "svg": svg_str,
                })
                continue

            c1_data = res1.get(s, {})
            c2_data = res2.get(s, {})
            c3_data = res3.get(s, {})

            p1 = c1_data.get("endpoint", NIK_MODEL_CONFIG["Run-1"]["c0"])
            p2 = c2_data.get("endpoint", NIK_MODEL_CONFIG["Run-2"]["c0"])
            p3 = c3_data.get("endpoint", NIK_MODEL_CONFIG["Run-3"]["c0"])

            dcw1 = c1_data.get("dcw", 0.0)
            dcw2 = c2_data.get("dcw", 0.0)
            dcw3 = c3_data.get("dcw", 0.0)

            def1 = c1_data.get("defect", 0.0) or 0.0
            def2 = c2_data.get("defect", 0.0) or 0.0
            def3 = c3_data.get("defect", 0.0) or 0.0

            avg_endpoint = round((p1 + p2 + p3) / 3.0, 4)
            avg_dcw = round((dcw1 + dcw2 + dcw3) / 3.0, 4)
            avg_defect = round((def1 + def2 + def3) / 3.0, 4)

            runs_dict = {
                "Run-1": {
                    "run": "Run-1",
                    "endpoint": round(p1, 4),
                    "smiles_weight": round(dcw1, 4),
                    "defect_smiles": round(def1, 4),
                },
                "Run-2": {
                    "run": "Run-2",
                    "endpoint": round(p2, 4),
                    "smiles_weight": round(dcw2, 4),
                    "defect_smiles": round(def2, 4),
                },
                "Run-3": {
                    "run": "Run-3",
                    "endpoint": round(p3, 4),
                    "smiles_weight": round(dcw3, 4),
                    "defect_smiles": round(def3, 4),
                },
            }

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
                    "svg_structure": props["svg_structure"],
                    "svg": props["svg_structure"],
                }
                svg_str = props["svg_structure"]

            results.append({
                "smiles": s,
                "target": "NF-κB Inducing Kinase (NIK / MAP3K14)",
                "consensus_prediction": avg_endpoint,
                "consensus_smiles_weight": avg_dcw,
                "consensus_defect_smiles": avg_defect,
                "defect_smiles": avg_defect,
                "in_domain": (avg_defect < THRESHOLD_DEFECT),
                "runs": runs_dict,
                "run1": round(p1, 4),
                "run2": round(p2, 4),
                "run3": round(p3, 4),
                "elapsed_seconds": total_elapsed,
                "physicochemical_properties": phys_props,
                "svg_structure": svg_str,
                "svg": svg_str,
            })

        return {
            "kinase": "nik",
            "results": results,
            "total_elapsed_seconds": total_elapsed,
        }

    def predict_single(self, smiles: str) -> Dict[str, Any]:
        res = self.predict_batch([smiles])
        if res["results"]:
            return res["results"][0]
        raise ValueError("Invalid SMILES input")


# Global singleton instance
nik_predictor_instance = NIKQSARPredictor()
predictor_instance = nik_predictor_instance
