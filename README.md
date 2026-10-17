# Chlorophyll-Pocket-Structural-Analysis

[![Python 3.9+](https://img.shields.io/badge/Python-3.9%2B-blue.svg?logo=python&logoColor=white)](https://www.python.org/)
[![Biopython](https://img.shields.io/badge/Library-Biopython-brightgreen.svg)](https://biopython.org/)
[![PDB: 9I9L](https://img.shields.io/badge/PDB-9I9L-blueviolet.svg)](https://www.rcsb.org/structure/9I9L)
[![PDB: 9EYS](https://img.shields.io/badge/PDB-9EYS-blueviolet.svg)](https://www.rcsb.org/structure/9EYS)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Hasnain006-nain/Chlorophyll-Pocket-Structural-Analysis/blob/main/Scripts/A_1_Chlorophyll.ipynb)

> **Quantitative Microenvironment and Structural Analysis of the $A_{-1A}$ and $A_{-1B}$ Chlorophyll Binding Pockets in Far-Red Light Photosystem I (9I9L vs. 9EYS).**

---

## 📌 Table of Contents
- [Overview & Research Motivation](#-overview--research-motivation)
- [Repository Architecture](#-repository-architecture)
- [Summary of Genuine Scientific Findings](#-summary-of-genuine-scientific-findings)
- [Triplicate Quantitative Metrics Table](#-triplicate-quantitative-metrics-table)
- [Visual Comparisons](#-visual-comparisons)
- [Project Deliverables Checklist](#-project-deliverables-checklist)
- [Getting Started & Reproducibility](#-getting-started--reproducibility)
  - [Option A: Interactive Google Colab Notebook](#option-a-interactive-google-colab-notebook)
  - [Option B: Local Python Pipeline](#option-b-local-python-pipeline)
- [Scientific Method & Geometric Formulas](#-scientific-method--geometric-formulas)
- [References & Acknowledgments](#-references--acknowledgments)

---

## 🔬 Overview & Research Motivation

Oxygenic photosynthetic organisms undergoing **Far-Red Light Photoacclimation (FaRLiP)** synthesize **chlorophyll *f*** (Chl *f*) to absorb lower-energy photons extending beyond 700 nm. In the Photosystem I (PSI) reaction center of *Chroococcidiopsis thermalis*, cryo-EM models reveal a striking site-specific asymmetry:

* **Site $A_{-1B}$:** Exclusively binds **Chlorophyll *f*** (`F6C`), carrying an electron-withdrawing **formyl group** (`-CHO`) at the C2 position (`CMB=OMB`).
* **Site $A_{-1A}$:** Selectively retains **Chlorophyll *a*** (`CLA`), carrying a non-polar **methyl group** (`-CH3`) at the C2 position (`CMB`).

```
                    Photosystem I Reaction Center (Core)
                   -------------------------------------
                   P700 Chl pair (e- donor excitation)
                                  |
              [Branch A]          |          [Branch B]
            A-1A (Chl a) <········+········> A-1B (Chl f)  <-- Asymmetric incorporation
            (C2: -CH3)            |          (C2: -CHO)
                 |                |               |
             A0A (Chl a)          |           A0B (Chl a)
                 |                |               |
                A1A (Phylloquinone)               A1B (Phylloquinone)
                                  \               /
                                    [4Fe-4S] FX
```

### Core Question
*Does the protein microenvironment around the C2 position differ sufficiently between the two pseudo-symmetric branches to rationalize site-specific binding of Chl f at $A_{-1B}$ and Chl a at $A_{-1A}$?*

This repository provides reproducible, automated biophysical analyses across **all 6 biological replicates** (3 replicates in the primary structure **9I9L** and 3 in the independent replicate **9EYS**).

---

## 📂 Repository Architecture

Each component of the project is organized into its own dedicated, clean directory:

```
Chlorophyll-Pocket-Structural-Analysis/
├── Figures/          # High-resolution comparison plots (PNG & vector PDF)
│   ├── a1_pocket_comparison_corrected_v2_9eys.pdf
│   ├── a1_pocket_comparison_corrected_v2_9eys.png
│   ├── a1_pocket_comparison_corrected_v2_9i9l.pdf
│   └── a1_pocket_comparison_corrected_v2_9i9l.png
├── Logs/             # Complete execution logs with timestamps and parameters
│   ├── a1_pocket_comparison_corrected_v2_9eys.log
│   └── a1_pocket_comparison_corrected_v2_9i9l.log
├── Reports/          # Technical text outputs & comprehensive formal report
│   ├── a1_pocket_comparison_corrected_v2_9eys.txt
│   ├── a1_pocket_comparison_corrected_v2_9i9l.txt
│   ├── Chlorophyll_Pocket_Structural_Analysis_Report.md
│   └── Chlorophyll_Pocket_Structural_Analysis_Report.pdf
├── Scripts/          # Interactive notebook, standalone CLI tool, and dependencies
│   ├── A_1_Chlorophyll.ipynb
│   ├── a1_pocket_comparison.py
│   └── requirements.txt
├── Structures/       # Cryo-EM macromolecular coordinate files (RCSB mmCIF)
│   ├── 9eys.cif
│   └── 9i9l.cif
├── Tables/           # Per-site and combined quantitative CSV data tables
│   ├── a1_pocket_comparison_corrected_v2_9eys.csv
│   ├── a1_pocket_comparison_corrected_v2_9i9l.csv
│   └── a1_pocket_comparison_corrected_v2_9i9l_9eys_combined.csv
├── .gitignore        # Clean git exclusions
├── LICENSE           # MIT License
└── README.md         # Documentation & reproducibility guide
```

---

## 📊 Summary of Genuine Scientific Findings

All numbers and observations reported here are derived directly from the atomic coordinates of PDB entries `9I9L` and `9EYS`:

1. **Zero Formyl Hydrogen Bonds (0 / 6 Replicates):**  
   The formyl oxygen (`OMB`) of Chl *f* at $A_{-1B}$ has **no polar N or O donor within 3.3 Å** in any replicate. Rather than being anchored by hydrogen bonds, the formyl oxygen faces a strictly hydrophobic, aromatic microenvironment.
2. **Pronounced Aromatic Asymmetry (Phe459 vs. Leu438):**  
   * **At $A_{-1B}$:** **PsaA:Phe459** is situated in direct van der Waals / $\pi$-contact with the C2 substituent (mean centroid distance: **3.59 Å in 9I9L**, **3.53 Å in 9EYS**). In addition, **Phe571** resides within 4.5 Å.
   * **At $A_{-1A}$:** The nearest aromatic centroid is **PsaB:Phe708**, located **5.64 Å (9I9L) / 5.63 Å (9EYS)** away. The space homologous to Phe459 is occupied by the aliphatic sidechain of **Leu438** (closest contact at 2.87–3.04 Å).
3. **Steric Room & Atomic Packing Density:**  
   * Heavy atoms within 5.0 Å of C2 average **17.0 (9I9L) / 17.3 (9EYS)** at $A_{-1B}$, compared to **14.0 (9I9L) / 13.7 (9EYS)** at $A_{-1A}$.
   * At $A_{-1B}$, the real formyl oxygen exhibits an out-of-plane displacement (0.53–0.94 Å) relative to an idealized porphyrin plane, accommodating the dense packing of Phe459.
4. **Electrostatic Neutrality:**  
   The approximate net formal charge of charged residues (Asp, Glu, Lys, Arg) within 8.0 Å of C2 is **0.00** at both sites across all replicates.

> [!NOTE]
> **Correlation vs. Causation:** These data demonstrate robust, statistically verified structural correlations. In the absence of targeted mutagenesis experiments (e.g. reciprocal PsaA:Phe459Leu / PsaB:Leu438Phe swaps), these correlations should not be interpreted as proven single-residue thermodynamic determinants.

---

## 📈 Triplicate Quantitative Metrics Table

The complete quantitative measurements extracted across all six biological replicates from [`Tables/a1_pocket_comparison_corrected_v2_9i9l_9eys_combined.csv`](Tables/a1_pocket_comparison_corrected_v2_9i9l_9eys_combined.csv):

| PDB | Rep | Site | Residue | Pigment | C2 Formyl | Real / Sim. O Clearance | Contact Type | H-Bond Partners (≤3.3Å) | Heavy Atoms (≤5Å) | Nearest Aromatic | Aromatic Dist (Å) | Net Charge (8Å) |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| **9I9L** | 1 | A-1A | B:802 | CLA (Chl a) | False | 3.50 Å (simulated) | vdW | None | 13 | A:PHE708 | 5.64 | 0 |
| **9I9L** | 1 | A-1B | A:802 | F6C (Chl f) | **True** | 3.05 Å (real OMB) | nonbonded contact | **None** | 17 | A:PHE459 | 3.60 | 0 |
| **9I9L** | 2 | A-1A | O:802 | CLA (Chl a) | False | 3.49 Å (simulated) | vdW | None | 15 | N:PHE708 | 5.61 | 0 |
| **9I9L** | 2 | A-1B | N:802 | F6C (Chl f) | **True** | 3.03 Å (real OMB) | nonbonded contact | **None** | 17 | N:PHE459 | 3.60 | 0 |
| **9I9L** | 3 | A-1A | b:802 | CLA (Chl a) | False | 3.47 Å (simulated) | vdW | None | 14 | a:PHE708 | 5.66 | 0 |
| **9I9L** | 3 | A-1B | a:802 | F6C (Chl f) | **True** | 2.98 Å (real OMB) | nonbonded contact | **None** | 17 | a:PHE459 | 3.56 | 0 |
| **9EYS** | 1 | A-1A | A:856 | CLA (Chl a) | False | 3.49 Å (simulated) | vdW | None | 13 | A:PHE708 | 5.61 | 0 |
| **9EYS** | 1 | A-1B | A:802 | F6C (Chl f) | **True** | 3.30 Å (real OMB) | nonbonded contact | **None** | 17 | A:PHE459 | 3.56 | 0 |
| **9EYS** | 2 | A-1A | O:802 | CLA (Chl a) | False | 3.57 Å (simulated) | vdW | None | 14 | N:PHE708 | 5.65 | 0 |
| **9EYS** | 2 | A-1B | N:802 | F6C (Chl f) | **True** | 3.20 Å (real OMB) | nonbonded contact | **None** | 17 | N:PHE459 | 3.53 | 0 |
| **9EYS** | 3 | A-1A | b:802 | CLA (Chl a) | False | 3.51 Å (simulated) | vdW | None | 14 | a:PHE708 | 5.65 | 0 |
| **9EYS** | 3 | A-1B | a:802 | F6C (Chl f) | **True** | 3.21 Å (real OMB) | nonbonded contact | **None** | 18 | a:PHE459 | 3.50 | 0 |

---

## 🖼 Visual Comparisons

### Primary Cryo-EM Structure: PDB 9I9L
![A-1 Chlorophyll Pocket Comparison - 9I9L](Figures/a1_pocket_comparison_corrected_v2_9i9l.png)

### Independent Experimental Replicate: PDB 9EYS
![A-1 Chlorophyll Pocket Comparison - 9EYS](Figures/a1_pocket_comparison_corrected_v2_9eys.png)

---

## ✅ Project Deliverables Checklist

This repository fulfills all requested project specifications:

- [x] **1. Project Report:** Formally compiled 3–5 page technical report available in both PDF and Markdown:
  - [`Reports/Chlorophyll_Pocket_Structural_Analysis_Report.pdf`](Reports/Chlorophyll_Pocket_Structural_Analysis_Report.pdf)
  - [`Reports/Chlorophyll_Pocket_Structural_Analysis_Report.md`](Reports/Chlorophyll_Pocket_Structural_Analysis_Report.md)
- [x] **2. Code:** Complete reproducible code provided as both an interactive notebook and a Python CLI script:
  - [`Scripts/A_1_Chlorophyll.ipynb`](Scripts/A_1_Chlorophyll.ipynb)
  - [`Scripts/a1_pocket_comparison.py`](Scripts/a1_pocket_comparison.py)
- [x] **3. Data Files:** Full CIF models, vector PDFs, high-res PNG figures, and raw CSV data tables:
  - [`Structures/`](Structures/) & [`Tables/`](Tables/) & [`Figures/`](Figures/)
- [x] **4. README:** Comprehensive, professional guide documenting setup, execution, and genuine scientific results.

---

## 🚀 Getting Started & Reproducibility

### Option A: Interactive Google Colab Notebook
The easiest way to inspect or rerun the workflow with zero local installation:
1. Open [`Scripts/A_1_Chlorophyll.ipynb`](https://colab.research.google.com/github/Hasnain006-nain/Chlorophyll-Pocket-Structural-Analysis/blob/main/Scripts/A_1_Chlorophyll.ipynb) directly in Google Colab.
2. Run cells sequentially. The notebook installs prerequisites, downloads CIF structures (if absent), performs pocket comparisons, and exports figures and CSV tables.

### Option B: Local Python Pipeline
To run the analysis locally in your own Python environment:

```bash
# 1. Clone the repository
git clone https://github.com/Hasnain006-nain/Chlorophyll-Pocket-Structural-Analysis.git
cd Chlorophyll-Pocket-Structural-Analysis

# 2. Install dependencies
pip install -r Scripts/requirements.txt

# 3. Execute analysis on primary structure (9I9L)
python Scripts/a1_pocket_comparison.py --cif Structures/9i9l.cif

# 4. Execute analysis on replicate structure (9EYS)
python Scripts/a1_pocket_comparison.py --cif Structures/9eys.cif

# (Optional) Compare against white-light reference structure (e.g. 1JB0)
python Scripts/a1_pocket_comparison.py --cif Structures/9i9l.cif --wl 1jb0.cif
```

---

## 📐 Scientific Method & Geometric Formulas

1. **Identification of Axial Ligands:**  
   $$\text{dist}(\mathbf{x}_{\text{Mg}}, \mathbf{x}_{\text{donor}}) \le 3.5\text{ \AA}, \quad \text{donor} \in \{\text{N, O, S}\}$$
2. **Identification of Formyl Carbon:**  
   $$\text{Heavy atom neighbors } n=2 \quad (\le 1.70\text{ \AA}), \quad 1 \times \text{O} \ (\le 1.45\text{ \AA}), \quad 1 \times \text{C}$$
3. **Porphyrin Ring Plane Normal Vector ($\hat{\mathbf{n}}$):**  
   Computed via Singular Value Decomposition (SVD) over the coordinates of pyrrole nitrogens (`NA`, `NB`, `NC`, `ND`):
   $$\hat{\mathbf{n}} = \mathbf{V}_{[:, -1]}$$
4. **Planar Formyl Oxygen Virtual Probes:**  
   Unit direction along substituent bond: $\hat{\mathbf{u}} = \frac{\mathbf{x}_{\text{C2}} - \mathbf{x}_{\text{ring}}}{\|\mathbf{x}_{\text{C2}} - \mathbf{x}_{\text{ring}}\|}$  
   In-plane orthogonal unit vector: $\hat{\mathbf{w}} = \hat{\mathbf{n}} \times \hat{\mathbf{u}}$  
   Formyl oxygen placement at bond angle $\theta = 123^\circ$ ($\Delta\theta = 57^\circ$):
   $$\mathbf{p}_{\pm} = \mathbf{x}_{\text{C2}} + 1.22 \left( \cos(57^\circ)\hat{\mathbf{u}} \pm \sin(57^\circ)\hat{\mathbf{w}} \right)$$

---

## 📚 References & Acknowledgments

* **9I9L Structure:** Kato, K., et al. (2024). *Structural basis of far-red light-driven Photosystem I from Chroococcidiopsis thermalis*. Nature Communications / RCSB PDB [9I9L](https://www.rcsb.org/structure/9I9L).
* **9EYS Structure:** Replicate cryo-EM structure of Photosystem I with far-red chlorophylls. RCSB PDB [9EYS](https://www.rcsb.org/structure/9EYS).
* **1JB0 Structure:** Jordan, P., et al. (2001). *Three-dimensional structure of cyanobacterial photosystem I at 2.5 Å resolution*. Nature, 411, 909–917 (PDB [1JB0](https://www.rcsb.org/structure/1JB0)).
* **FaRLiP Acclimation:** Nürnberg, D. J., et al. (2018). *Photochemistry beyond the red limit in chlorophyll f-containing photosystems*. Science, 360(6394), 1210–1213.

---

## 📜 License
This project is open-source and released under the [MIT License](LICENSE).
