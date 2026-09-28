# 🚀 24/7 Cloud Deployment Guide (Option 2: Docker + Wine32)

This guide explains how to deploy the authentic **CORALSEA.exe** Monte Carlo prediction engine on a free 24/7 cloud server (Hugging Face Spaces or Render) using Docker, Wine32, and Xvfb.

---

## 🏗️ Architecture Overview

```
                      +---------------------------------------+
                      |   Client Web Browser (Any Device)     |
                      |   https://nik-ddsl.streamlit.app/     |
                      +---------------------------------------+
                                          |
                                          | (POST /api/predict)
                                          v
                      +---------------------------------------+
                      |  Hugging Face Spaces (Free Docker)    |
                      |  Ubuntu 22.04 + Wine32 + Xvfb        |
                      |                                       |
                      |  FastAPI (Port 7860)                  |
                      |       |                               |
                      |       v                               |
                      |  coral_runner.exe (Win32 Driver)      |
                      |       |                               |
                      |       v (PostMessage BM_CLICK)        |
                      |  CORALSEA.exe (Run-1, Run-2, Run-3)   |
                      |       |                               |
                      |       v (Exact Model Output)          |
                      |  ListModel.txt & DemoDCW.txt          |
                      +---------------------------------------+
```

---

## 📋 Steps to Deploy on Hugging Face Spaces (100% Free, 24/7)

### Step 1: Create a Free Hugging Face Space
1. Log in to [Hugging Face](https://huggingface.co/).
2. Click your profile avatar (top right) → **New Space**.
3. Fill in the details:
   - **Space name**: `nik-qsar-engine`
   - **License**: `mit`
   - **Space SDK**: Select **Docker** → **Blank**
   - **Visibility**: `Public`
4. Click **Create Space**.

---

### Step 2: Push Repository to Hugging Face
In your local project directory (`C:\Users\drugd\Desktop\Webtool\NIK_2D_QSAR\nik_qsar_web`):

```bash
# Add Hugging Face Space as a remote
git remote add space https://huggingface.co/spaces/<YOUR_HF_USERNAME>/nik-qsar-engine

# Push main branch to Hugging Face Space
git push space main
```

*(Alternatively, you can connect your GitHub repository directly to Hugging Face Spaces).*

---

### Step 3: Space Automatically Builds & Launches
Hugging Face will detect `Dockerfile` and execute:
1. Ubuntu 22.04 base installation.
2. 32-bit Wine architecture (`wine32`).
3. Virtual X display server (`Xvfb`).
4. MinGW C compilation of `coral_runner.c`.
5. Pre-configured models (`Run-1`, `Run-2`, `Run-3`) with `CORALSEA.exe`.
6. Startup of FastAPI server on port 7860.

Once the status turns **Running**, your live cloud endpoint will be:
`https://<YOUR_HF_USERNAME>-nik-qsar-engine.hf.space`

---

### Step 4: Connect Streamlit Cloud to the Cloud Engine
1. Go to your Streamlit Cloud dashboard: [https://share.streamlit.io](https://share.streamlit.io).
2. Click the `...` menu next to `nik-ddsl` → **Settings** → **Secrets**.
3. Paste the following configuration:
   ```toml
   CORAL_CLOUD_API = "https://<YOUR_HF_USERNAME>-nik-qsar-engine.hf.space"
   ```
4. Click **Save**.

---

## ✅ Verification
1. Open `https://nik-ddsl.streamlit.app/`.
2. Input any test compound, for example:
   `Nc1nc(N2CCc3c2cc(O)cc3)c(Cl)cn1` (TEST0008)
3. Click **Predict NIK Activity**.
4. The prediction will be computed by the real `CORALSEA.exe` in the cloud container:
   - **Run-1**: `7.1758`
   - **Run-2**: `7.1902`
   - **Run-3**: `7.4105`
   - **Consensus Average**: `7.2588`
   - **Defect SMILES**: `0.2816` (In Domain: Reliable)
