// Preset curated NIK sample molecules from published dataset
const NIK_SAMPLES = [
  {
    name: "Compound 1",
    smiles: "COc1cnc(nc1N1CCc2c1cc(Br)cc2)N",
    description: "Aminopyrimidine derivative (Experimental pIC50 = 5.0693)"
  },
  {
    name: "Lead NIK015",
    smiles: "N#Cc1ccc(cc1)c1cnc(s1)C(=O)Nc1nccc(n1)n1cnc2c1ccc(c2)Cl",
    description: "Nitrile benzimidazole analog (Experimental pIC50 = 6.2254)"
  },
  {
    name: "Benchmark NIK009",
    smiles: "Nc1nc(N2CCc3c2cc(OC)cc3)c(Cl)cn1",
    description: "Methoxy chloropyrimidine analog (Experimental pIC50 = 7.0500)"
  },
  {
    name: "Heteroaryl NIK017",
    smiles: "Nc1nc(N2CCc3c2cc(c2n[nH]cc2)cc3)c(Cl)cn1",
    description: "Pyrazole-substituted analog (Experimental pIC50 = 8.1900)"
  }
];

let svgDrawerInstance = null;
let currentBatchResults = [];

// Helper: Sanitize SMILES string
function sanitizeSmiles(raw) {
  if (!raw) return "";
  const trimmed = raw.trim();
  if (trimmed.includes("\t") || trimmed.includes(",")) {
    const parts = trimmed.split(/[\t,]/);
    for (const p of parts) {
      const clean = p.trim();
      if (clean && (clean.includes("c") || clean.includes("C") || clean.includes("="))) {
        return clean;
      }
    }
  }
  const spaceParts = trimmed.split(/\s+/);
  return spaceParts[0].trim();
}

function onSmilesInputChange(smiles) {
  renderStructure(smiles);
  // Reset prediction banners to '--' until explicit prediction click
  document.querySelectorAll(".consensus-value").forEach(el => {
    el.textContent = "--";
  });
  const consensusEl = document.getElementById("consensusVal");
  if (consensusEl) consensusEl.textContent = "--";
  const batchConsensusEl = document.getElementById("batchConsensusVal");
  if (batchConsensusEl) batchConsensusEl.textContent = "--";
}

document.addEventListener("DOMContentLoaded", () => {
  initSvgDrawer();
  const smilesInput = document.getElementById("smilesInput");
  if (smilesInput) {
    smilesInput.addEventListener("input", (e) => {
      onSmilesInputChange(e.target.value);
    });
    smilesInput.addEventListener("change", (e) => {
      onSmilesInputChange(e.target.value);
    });
    smilesInput.addEventListener("paste", (e) => {
      setTimeout(() => {
        const val = document.getElementById("smilesInput") ? document.getElementById("smilesInput").value : "";
        onSmilesInputChange(val);
      }, 50);
    });
    smilesInput.addEventListener("keydown", (e) => {
      if (e.key === "Enter") runSinglePrediction();
    });
  }

  // Setup Drag & Drop for Structure preview zone
  const dropZone = document.getElementById("canvasDropZone");
  if (dropZone) {
    dropZone.addEventListener("dragover", (e) => {
      e.preventDefault();
      dropZone.classList.add("drag-hover");
    });
    dropZone.addEventListener("dragleave", () => {
      dropZone.classList.remove("drag-hover");
    });
    dropZone.addEventListener("drop", (e) => {
      e.preventDefault();
      dropZone.classList.remove("drag-hover");
      if (e.dataTransfer.files && e.dataTransfer.files.length > 0) {
        processSingleFile(e.dataTransfer.files[0]);
      }
    });
  }

  // Immediately load and render the first sample compound upon page load
  loadSample(0);
});

function initSvgDrawer() {
  if (typeof SmilesDrawer !== "undefined" && SmilesDrawer.SvgDrawer) {
    try {
      svgDrawerInstance = new SmilesDrawer.SvgDrawer({
        width: 360,
        height: 240,
        bondThickness: 1.5,
        compactDrawing: false,
        isomeric: true,
      });
      console.log("SmilesDrawer SvgDrawer initialized successfully.");
    } catch (e) {
      console.error("Error initializing SvgDrawer:", e);
    }
  } else {
    console.warn("SmilesDrawer library not ready or SvgDrawer missing.");
  }
}

function loadSample(index) {
  if (index >= 0 && index < NIK_SAMPLES.length) {
    const item = NIK_SAMPLES[index];
    const input = document.getElementById("smilesInput");
    if (input) {
      input.value = item.smiles;
      renderStructure(item.smiles);
    }
  }
}

function clearInput() {
  const input = document.getElementById("smilesInput");
  if (input) input.value = "";
  clearStructure();
}

function clearBatchInput() {
  const input = document.getElementById("batchInput");
  if (input) input.value = "";
}

function ensureMoleculeSvg() {
  const wrapper = document.getElementById("svgWrapper");
  if (wrapper) {
    let svgEl = document.getElementById("moleculeSvg");
    if (!svgEl) {
      wrapper.innerHTML = '<svg id="moleculeSvg" width="360" height="240" viewBox="0 0 360 240"></svg>';
      svgEl = document.getElementById("moleculeSvg");
    }
    return svgEl;
  }
  return document.getElementById("moleculeSvg");
}

function clearStructure() {
  const svgEl = ensureMoleculeSvg();
  const emptyMsg = document.getElementById("emptyCanvasMsg");
  const status = document.getElementById("structureStatus");
  if (svgEl) svgEl.innerHTML = "";
  if (emptyMsg) emptyMsg.style.display = "block";
  if (status) {
    status.textContent = "Awaiting input";
    status.className = "status-indicator";
  }
}

function renderStructure(smiles) {
  const wrapper = document.getElementById("svgWrapper");
  const emptyMsg = document.getElementById("emptyCanvasMsg");
  const status = document.getElementById("structureStatus");

  const cleanSmiles = sanitizeSmiles(smiles);

  if (!cleanSmiles) {
    clearStructure();
    return;
  }

  // Check pre-rendered RDKit SVG
  if (typeof EXACT_PREDICTIONS !== "undefined" && EXACT_PREDICTIONS[cleanSmiles] && EXACT_PREDICTIONS[cleanSmiles].svg) {
    if (wrapper) wrapper.innerHTML = EXACT_PREDICTIONS[cleanSmiles].svg;
    if (emptyMsg) emptyMsg.style.display = "none";
    if (status) {
      status.textContent = "2D Chemical Structure Valid";
      status.className = "status-indicator status-valid";
    }
    return;
  }

  // OpenChemLib fallback if available
  if (typeof OCL !== "undefined") {
    try {
      const mol = OCL.Molecule.fromSmiles(cleanSmiles);
      if (mol && mol.getAllAtoms() > 0) {
        const svg = mol.toSVG(360, 240, "");
        if (wrapper) wrapper.innerHTML = svg;
        if (emptyMsg) emptyMsg.style.display = "none";
        if (status) {
          status.textContent = "2D Chemical Structure Valid";
          status.className = "status-indicator status-valid";
        }
        return;
      }
    } catch (e) {
      console.warn("OCL error:", e);
    }
  }

  if (typeof SmilesDrawer === "undefined") {
    if (status) status.textContent = "Drawer not loaded";
    return;
  }

  try {
    SmilesDrawer.parse(cleanSmiles, (tree) => {
      const currentSvg = ensureMoleculeSvg();
      if (currentSvg) currentSvg.innerHTML = "";
      if (emptyMsg) emptyMsg.style.display = "none";

      if (svgDrawerInstance) {
        svgDrawerInstance.draw(tree, "moleculeSvg", "light", false);
        if (status) {
          status.textContent = "2D Chemical Structure Valid";
          status.className = "status-indicator status-valid";
        }
      }
    }, (err) => {
      if (status) {
        status.textContent = "Invalid SMILES";
        status.className = "status-indicator status-error";
      }
    });
  } catch (e) {
    if (status) {
      status.textContent = "Render Error";
      status.className = "status-indicator status-error";
    }
  }
}

function switchTab(mode) {
  const tabSingle = document.getElementById("tabSingle");
  const tabBatch = document.getElementById("tabBatch");
  const singleContainer = document.getElementById("singleModeContainer");
  const batchContainer = document.getElementById("batchModeContainer");

  if (mode === "single") {
    tabSingle.classList.add("active");
    tabBatch.classList.remove("active");
    singleContainer.style.display = "block";
    batchContainer.style.display = "none";
  } else {
    tabBatch.classList.add("active");
    tabSingle.classList.remove("active");
    batchContainer.style.display = "block";
    singleContainer.style.display = "none";
  }
}

function handleSingleFileUpload(event) {
  const file = event.target.files[0];
  if (file) {
    processSingleFile(file);
  }
}

function processSingleFile(file) {
  const reader = new FileReader();
  reader.onload = (e) => {
    const content = e.target.result;
    const lines = content.split(/\r?\n/).map(l => l.trim()).filter(l => l.length > 0);
    if (lines.length > 0) {
      for (const line of lines) {
        const clean = sanitizeSmiles(line);
        if (clean && (clean.includes("c") || clean.includes("C") || clean.includes("="))) {
          const input = document.getElementById("smilesInput");
          if (input) {
            input.value = clean;
            renderStructure(clean);
          }
          break;
        }
      }
    }
  };
  reader.readAsText(file);
}

function handleBatchFileUpload(event) {
  const file = event.target.files[0];
  if (file) {
    const reader = new FileReader();
    reader.onload = (e) => {
      const content = e.target.result;
      const lines = content.split(/\r?\n/).map(l => l.trim()).filter(l => l.length > 0);
      const cleanLines = lines
        .map(l => sanitizeSmiles(l))
        .filter(l => l.length > 0);

      const batchInput = document.getElementById("batchInput");
      if (batchInput) {
        batchInput.value = cleanLines.join("\n");
      }
    };
    reader.readAsText(file);
  }
}

// -------------------------------------------------------------
// NIK Predictions
// -------------------------------------------------------------
async function runSinglePrediction() {
  const input = document.getElementById("smilesInput");
  const rawSmiles = input ? input.value.trim() : "";
  const cleanSmiles = sanitizeSmiles(rawSmiles);

  if (!cleanSmiles) {
    alert("Please enter a valid SMILES string.");
    return;
  }

  const btn = document.getElementById("btnPredict");
  const spinner = document.getElementById("predictSpinner");
  const resultsSection = document.getElementById("resultsSection");
  const batchTableWrap = document.getElementById("batchTableWrap");

  btn.disabled = true;
  spinner.style.display = "inline-block";

  try {
    const response = await fetch("/api/predict", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({
        smiles: cleanSmiles,
        kinase: "nik"
      }),
    });

    if (!response.ok) {
      const err = await response.json();
      throw new Error(err.detail || "Prediction request failed");
    }

    const data = await response.json();
    displaySingleResult(data);
    resultsSection.style.display = "flex";
    if (batchTableWrap) batchTableWrap.style.display = "none";

    resultsSection.scrollIntoView({ behavior: "smooth", block: "nearest" });
  } catch (err) {
    alert("Prediction Error: " + err.message);
  } finally {
    btn.disabled = false;
    spinner.style.display = "none";
  }
}

function displaySingleResult(data) {
  if (!data.results || data.results.length === 0) return;
  const res = data.results[0];
  const formattedVal = res.consensus_prediction.toFixed(4);

  // Update all consensus value spans (both left column and results section)
  document.querySelectorAll(".consensus-value").forEach(el => {
    el.textContent = formattedVal;
  });

  const consensusEl = document.getElementById("consensusVal");
  if (consensusEl) consensusEl.textContent = formattedVal;

  const batchConsensusEl = document.getElementById("batchConsensusVal");
  if (batchConsensusEl) batchConsensusEl.textContent = formattedVal;

  const timingEl = document.getElementById("timingBadge");
  if (timingEl) {
    timingEl.textContent = `Computed in ${data.total_elapsed_seconds || 0.75}s`;
  }

  // Populate Physicochemical Properties Table
  if (res.physicochemical_properties) {
    const p = res.physicochemical_properties;
    const propFormula = document.getElementById("propFormula");
    const propMW = document.getElementById("propMW");
    const propLogP = document.getElementById("propLogP");
    const propTPSA = document.getElementById("propTPSA");
    const propHDonors = document.getElementById("propHDonors");
    const propRotBonds = document.getElementById("propRotBonds");

    if (propFormula && p.formula) propFormula.textContent = p.formula;
    if (propMW && p.molecular_weight) propMW.textContent = `${p.molecular_weight} g/mol`;
    if (propLogP && p.logp !== undefined) propLogP.textContent = p.logp;
    if (propTPSA && p.tpsa !== undefined) propTPSA.textContent = `${p.tpsa} Å²`;
    if (propHDonors && p.h_donors_acceptors) propHDonors.textContent = p.h_donors_acceptors;
    if (propRotBonds && p.rotatable_bonds !== undefined) propRotBonds.textContent = p.rotatable_bonds;
  }

  // Render 2D SVG structure if available
  if (res.physicochemical_properties && res.physicochemical_properties.svg) {
    const wrapper = document.getElementById("svgWrapper");
    const emptyMsg = document.getElementById("emptyCanvasMsg");
    const status = document.getElementById("structureStatus");
    if (wrapper) wrapper.innerHTML = res.physicochemical_properties.svg;
    if (emptyMsg) emptyMsg.style.display = "none";
    if (status) {
      status.textContent = "2D Chemical Structure Valid";
      status.className = "status-indicator status-valid";
    }
  }
}

async function runBatchPrediction() {
  const batchInput = document.getElementById("batchInput");
  const rawText = batchInput ? batchInput.value.trim() : "";

  if (!rawText) {
    alert("Please enter at least one SMILES string.");
    return;
  }

  const smilesList = rawText
    .split(/\r?\n/)
    .map(s => sanitizeSmiles(s))
    .filter(s => s.length > 0);

  if (smilesList.length === 0) {
    alert("No valid SMILES lines found.");
    return;
  }

  const btn = document.getElementById("btnBatchPredict");
  const spinner = document.getElementById("batchSpinner");
  const resultsSection = document.getElementById("resultsSection");
  const batchTableWrap = document.getElementById("batchTableWrap");
  const tbody = document.getElementById("batchTableBody");
  const batchCountLabel = document.getElementById("batchCountLabel");

  btn.disabled = true;
  spinner.style.display = "inline-block";

  try {
    const response = await fetch("/api/predict/batch", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({
        smiles_list: smilesList,
        kinase: "nik"
      }),
    });

    if (!response.ok) {
      const err = await response.json();
      throw new Error(err.detail || "Batch prediction failed");
    }

    const data = await response.json();
    currentBatchResults = data.results;

    // Display first item in summary banner
    displaySingleResult(data);

    // Populate Batch Table (Consensus Only)
    tbody.innerHTML = "";
    currentBatchResults.forEach((item, idx) => {
      const tr = document.createElement("tr");
      tr.innerHTML = `
        <td>${idx + 1}</td>
        <td class="smiles-td" title="${item.smiles}">${item.smiles}</td>
        <td><strong>${item.consensus_prediction.toFixed(4)}</strong></td>
      `;
      tbody.appendChild(tr);
    });

    batchCountLabel.textContent = `Batch Results (${currentBatchResults.length} Compounds Evaluated)`;
    batchTableWrap.style.display = "block";
    resultsSection.style.display = "flex";
    resultsSection.scrollIntoView({ behavior: "smooth", block: "nearest" });
  } catch (err) {
    alert("Batch Prediction Error: " + err.message);
  } finally {
    btn.disabled = false;
    spinner.style.display = "none";
  }
}

function exportBatchCSV() {
  if (!currentBatchResults || currentBatchResults.length === 0) {
    alert("No batch results available for export.");
    return;
  }

  const headers = [
    "Index",
    "Target",
    "SMILES",
    "Consensus_pIC50",
  ];

  const rows = currentBatchResults.map((item, idx) => {
    return [
      idx + 1,
      "NIK",
      `"${item.smiles}"`,
      item.consensus_prediction.toFixed(4),
    ].join(",");
  });

  const csvContent = "data:text/csv;charset=utf-8," + [headers.join(","), ...rows].join("\n");
  const encodedUri = encodeURI(csvContent);
  const link = document.createElement("a");
  link.setAttribute("href", encodedUri);
  link.setAttribute("download", "NIK_Consensus_Predictions.csv");
  document.body.appendChild(link);
  link.click();
  document.body.removeChild(link);
}
