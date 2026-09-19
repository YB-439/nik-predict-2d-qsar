# NIK-Predict: NF-κB Inducing Kinase Bioactivity Prediction Platform

[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://share.streamlit.io)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

An open-source, high-throughput **Streamlit** and **FastAPI** platform for predicting the biological inhibitory activity (\(pIC_{50}\)) of molecules against **NF-κB Inducing Kinase (NIK)** using a calibrated 3-run Monte Carlo 2D-QSAR model.

---

## 🏛️ Attribution & Laboratory Details

- **Laboratory**: Drug Design & Synthesis Lab (*In-Silico to In-Vivo*)
- **Department**: Department of Pharmaceutical Sciences and Drug Research
- **Institution**: Punjabi University, Patiala, Punjab, India
- **Contact Emails**: `drugdesignsynthesislab@gmail.com` | `Yogita_pharma@pbi.ac.in`

---

## 🚀 Live Streamlit Cloud Deployment

- **Live Application URL**: [https://nik-ddsl.streamlit.app/](https://nik-ddsl.streamlit.app/)

### Standard Deployment Instructions (`https://share.streamlit.io`):

1. Repository: [`https://github.com/YB-439/nik-predict-2d-qsar`](https://github.com/YB-439/nik-predict-2d-qsar)
2. Main file path: `streamlit_app.py`

---

## 💻 Local Installation & Setup

### Prerequisites
- Python 3.10+ installed on Windows, macOS, or Linux.

### 1. Clone Repository
```bash
git clone https://github.com/YB-439/nik-predict-2d-qsar.git
cd nik-predict-2d-qsar
```

### 2. Install Dependencies
```bash
pip install -r requirements.txt
```

### 3. Launch Streamlit Web App
```bash
streamlit run streamlit_app.py
```
Open your browser to [http://localhost:8501](http://localhost:8501).

---

## 🔬 Model Performance & Specifications

- **Selected Model**: Model M3 (CORALSEA Correlation Balance)
- **Validation Statistics**:
  - $R^2 = 0.765$
  - $Q^2 = 0.660$
  - $	ext{MAE} = 0.291$

---

## 📄 License
This project is open-sourced under the MIT License for academic and drug discovery research.
