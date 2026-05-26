# CropPhenoAI (v1.0.0-Beta)
### An Open-Source Localized System Architecture for High-Throughput Phenomics & Multi-Trait Selection
👉 **[CLICK HERE TO TRY THE LIVE INTERACTIVE APP](https://susmita-stha.shinyapps.io/CropPhenoAI/)** 👈

**Author:** Susmita Shrestha 
**Status:** Operational Architecture Prototype / Open-Source Beta Release

---

## 📋 System Overview
CropPhenoAI is an interactive software infrastructure designed to bridge the gap between computer vision and quantitative plant breeding. Built entirely as a secure, local-first framework using R Shiny and Python via `reticulate`, the platform allows researchers to process specimen imagery alongside complex agronomic field trials on a standard workstation with zero external data-leak risks.

### Core Architecture Modules
1. **Automated Image Phenomics:** Local digital matrix processing utilizing OpenCV contour-tracking loops, supported by an interactive spatial pixel calibration layer.
2. **Research Design Generation:** Automated randomized layout configurations for RCBD and Incomplete Block (Alpha Lattice) designs utilizing `agricolae`.
3. **Spatial Single-Trial Analysis:** Mixed-Model Variance partitioning using Linear Mixed Models (LMM) via Residual Maximum Likelihood (REML) to extract true phenotypic Adjusted Means (BLUEs) across spatial field gradients.
4. **Combined GxE Selection Engine:** Synchronous dual-index screening, calculating the Multi-Trait Genotype-Ideotype Distance Index (MGIDI) and the Stress Tolerance Index (STI) for resilience mapping.

---

## 🛠️ Development Methodology & AI Disclosure
This framework represents a modern, human-in-the-loop co-development engineering workflow:
- **System Architecture & Breeding Logic:** Conceptualized, structurally mapped, and algorithmically designed by the author (defining the targeted crop trait matrix, directional ranking parameters, mixed-model formulas, and index ideotype distances).
- **Code Compilation & Interface Layout:** Accelerated utilizing advanced generative AI modeling pipelines (**Google AI Studio**) to translate biological criteria into operational, fully synthesized reactive R scripts. Because the final infrastructure was compiled via an iterative design and debugging workflow, an exact singular prompt does not exist.

---

## 🚀 Current Development Stage & Validation Roadmap
The current codebase is published as a fully operational **Functional Prototype Framework**. 

- [x] **UI/UX Core Dashboard:** Completed and fully responsive.
- [x] **Statistical Analytics Engines:** Algorithms for REML, ANOVA, STI, and MGIDI are mathematically complete and verified using synthetic data distributions.
- [x] **Local Python Interoperability:** `reticulate` virtual environment generation and OpenCV image matrix parsing logic are completely implemented.
- [ ] **Field & Lab Validation (Current Development Sprint):** The image extraction matrix currently utilizes placeholder biological random generation to maintain reactivity while a manual calibration slider handles pixel scaling. The next milestone involves calibrating the computer vision contours against physical caliper measurements of real wheat spike, maize cob, and rice panicle samples to replace the fallback baselines.

---

## 🔧 Installation & Deployment
To execute this platform locally, initialize your R console and install the following runtime dependencies:

```R
install.packages(c("shiny", "shinydashboard", "DT", "shinycssloaders", 
                   "readxl", "writexl", "tidyverse", "scales", 
                   "agricolae", "lmerTest", "emmeans", "reticulate"))
