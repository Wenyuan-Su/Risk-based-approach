# A risk-based approach prioritizes carcinogens from the exposome atlas of indoor dust

<p align="center">
  <img src="https://img.shields.io/badge/Python-Data--driven%20modeling-3776AB?logo=python&logoColor=white" alt="Python modeling">
  <img src="https://img.shields.io/badge/R-SMACH%20analysis-276DC3?logo=r&logoColor=white" alt="R SMACH analysis">
  <img src="https://img.shields.io/badge/HRMS-GC--EI%20%26%20LC--ESI-5A67D8" alt="GC-EI and LC-ESI HRMS">
</p>

---

## 🧭 Overview

This repository contains computational workflows supporting the study:

**A risk-based approach prioritizes carcinogens from the exposome atlas of indoor dust**

The study integrates data-driven carcinogenicity prediction, mass spectrometric analysis, mixture exposure experiments, and retrospective exploration of environmental occurrence. Predicted potential carcinogens guide chemical screening in indoor dust, while representative chemical mixtures are evaluated experimentally to identify contributors to the observed mutagenicity.

The study analyzed **135 indoor dust samples** from seven geographic regions of China and prepared **14 composite extracts** for mixture exposure experiments. The ensemble model predicted **1,984 potential carcinogens**, and the mass spectrometric analysis identified **294 potential carcinogens** in dust.

The computational files are organized into three main modules:

| Module | Description | Language |
|---|---|---|
| 🧠 `1. Data-Driven Modeling/` | PubChem bioassay data collection, endpoint classifiers, ensemble prediction, applicability domain characterization, and structural alert analysis | Python, R |
| 🔬 `2. Mass Spectrometric Analysis/` | HRMS feature processing, EI spectral prediction with NEIMS, and candidate ranking with ChemWalker | Python, R |
| 🧪 `3. Similar Mixture Approach/` | Toxicity-weighted comparison of candidate mixtures with a reference mixture using SMACH distances | R |

Experimental conditions, chemical identities, concentration data, and detailed mixture compositions are described in the article and Supplementary Information.

---

## 📁 Repository Structure

The main analysis files are shown below.

```text
.
├── README.md
├── 1. Data-Driven Modeling/
│   ├── 1. Data Collection/
│   │   ├── PUG-REST.ipynb
│   │   ├── pubchem_assays聚类分析图.R
│   │   ├── pubchem-cid.txt
│   │   ├── bioprofile_long.csv
│   │   └── bioprofile_matrix.csv
│   ├── 2. Model Construction/
│   │   ├── AID_1189.ipynb
│   │   ├── AID_1194.ipynb
│   │   ├── AID_1199.ipynb
│   │   ├── AID_1205.ipynb
│   │   ├── AID_1208.ipynb
│   │   ├── AID_1259407.ipynb
│   │   ├── AID_1259408.ipynb
│   │   ├── AID_1259411.ipynb
│   │   └── Ensemble_model.ipynb
│   ├── 3. Applicability Domain Characterization/
│   │   └── Application_domain.ipynb
│   ├── 4. Structural Alert Analysis/
│   │   ├── 模型可解释性.py
│   │   └── 模型可解释性-高亮原子和化学键.py
│   ├── external_validation_standardized.csv
│   └── model_comparasion.csv
├── 2. Mass Spectrometric Analysis/
│   ├── HRMS_data_processing.R
│   ├── NEIMS/
│   │   ├── README.md
│   │   ├── GC_MS2_predict.ipynb
│   │   └── make_spectra_prediction.py
│   └── ChemWalker/
│       ├── README.md
│       ├── environment.yml
│       └── notebooks/
│           └── run_chemwalker.ipynb
└── 3. Similar Mixture Approach/
    └── SMACH distance calculation.R
```

---

# 🧠 1. Data-Driven Modeling

This module supports carcinogenicity prediction by selecting bioassays associated with known carcinogenicity classifications, modeling their endpoints, and combining their predictions.

### Main workflow

```text
IRIS carcinogenicity classifications
        ↓
Selection of associated PubChem bioassays
        ↓
Molecular fingerprints and bioactivity descriptors
        ↓
Endpoint-specific model development and evaluation
        ↓
Ensemble prediction across eight endpoints
        ↓
Applicability domain and structural alert analysis
        ↓
Potential carcinogens for mass spectrometric screening
```

The manuscript uses three criteria to select PubChem bioassays: at least **10 active probe compounds**, a significant association with carcinogenicity classifications (**Fisher's exact test, P < 0.05**), and a **correct classification rate > 0.65**.

The eight enrolled assays are **AID 1189, 1194, 1199, 1205, 1208, 1259407, 1259408, and 1259411**.

Model development evaluates:

- MACCS keys, extended-connectivity fingerprints, functional-class fingerprints, and Chemical Checker Signaturizer bioactivity descriptors;
- random forest, multilayer perceptron, k-nearest neighbor, support vector machine, and gradient boosting classifiers;
- SMOTE for class imbalance, Boruta for feature selection, and randomized hyperparameter search; and
- five-fold cross-validation and independent external validation.

The ensemble prediction probability is the mean probability across the eight endpoint models. Chemicals with a probability **greater than 0.5** are classified as potential carcinogens. The model is applied to chemicals compiled from **IECSC, TSCA, DSL, and REACH** inventories.

Applicability domain analysis uses structural similarity and local discontinuity to identify unreliable predictions. Structural alert analysis uses Bemis–Murcko scaffold enrichment to identify scaffolds associated with carcinogenicity.

### Main files

| File | Purpose |
|---|---|
| [PUG-REST.ipynb](1.%20Data-Driven%20Modeling/1.%20Data%20Collection/PUG-REST.ipynb) | Retrieve PubChem bioassay outcomes and organize chemical bioactivity profiles |
| `AID_*.ipynb` in `2. Model Construction/` | Develop and evaluate the eight endpoint classifiers |
| [Ensemble_model.ipynb](1.%20Data-Driven%20Modeling/2.%20Model%20Construction/Ensemble_model.ipynb) | Combine endpoint predictions and evaluate ensemble performance |
| [Application_domain.ipynb](1.%20Data-Driven%20Modeling/3.%20Applicability%20Domain%20Characterization/Application_domain.ipynb) | Characterize prediction reliability using structural similarity and local discontinuity |
| [模型可解释性.py](1.%20Data-Driven%20Modeling/4.%20Structural%20Alert%20Analysis/模型可解释性.py) | Analyze carcinogenicity-associated molecular scaffolds |
| [模型可解释性-高亮原子和化学键.py](1.%20Data-Driven%20Modeling/4.%20Structural%20Alert%20Analysis/模型可解释性-高亮原子和化学键.py) | Visualize structural alerts in molecular structures |

### Inputs and outputs

Inputs include carcinogenicity classifications, PubChem bioassay outcomes, molecular structures, endpoint datasets, external validation data, and chemical inventories. Outputs include bioactivity profiles, endpoint and ensemble predictions, performance summaries, applicability domain assessments, and structural alerts.

Detailed procedures are provided in **Methods 4.3** and **Supplementary Notes 1 and 2**.

---

# 🔬 2. Mass Spectrometric Analysis

This module supports the identification of potential carcinogens using complementary **GC-EI HRMS** and **LC-ESI HRMS** analyses. The predicted carcinogens provide the suspect list for screening.

### GC-EI: Compound-to-MS workflow

```text
Potential carcinogen structures
        ↓
NEIMS prediction of EI mass spectra
        +
Prediction of retention indices
        ↓
Predicted spectral library
        ↓
GC-EI feature matching and structural annotation
```

The manuscript combines predicted EI spectra and retention indices in an MSP-format library. GC-EI features are compared with experimental libraries, including **NIST 20** and **MoNA GC**, and with the predicted library. The annotation workflow is evaluated using a **507-compound benchmark dataset**.

### LC-ESI: MS-to-Compound workflow

```text
LC-ESI HRMS features and MS/MS spectra
        ↓
Target, suspect, and nontarget screening
        ↓
Spectral library matching or MetFrag candidate ranking
        ↓
Molecular networking
        ↓
ChemWalker candidate re-ranking
        +
ModiFinder modification-site localization
        ↓
Structural annotation and concentration estimation
```

The manuscript describes feature extraction using **MZmine**, candidate generation with **MetFrag** and **PubChemLite**, and molecular formula generation with **GenForm**. Retained suspect features have an absolute abundance **> 10,000** and a mass error **within 5 ppm**. Molecular networking requires a modified cosine score **> 0.5** and at least **two matched fragment ions**.

ChemWalker propagates structural information from verified seed nodes through connected molecular networks using random walks. ModiFinder supports the investigation of structural analogues by locating modifications from shifted fragment ions. Compounds with reference standards are quantified using calibration curves; other compounds are semi-quantified using structurally similar standards.

### Main files

| File or directory | Purpose |
|---|---|
| [HRMS_data_processing.R](2.%20Mass%20Spectrometric%20Analysis/HRMS_data_processing.R) | HRMS feature processing, suspect screening, formula generation, and structural candidate annotation |
| [NEIMS/GC_MS2_predict.ipynb](2.%20Mass%20Spectrometric%20Analysis/NEIMS/GC_MS2_predict.ipynb) | Run EI spectral prediction for input molecular structures |
| [NEIMS/](2.%20Mass%20Spectrometric%20Analysis/NEIMS/) | NEIMS prediction software and associated files |
| [ChemWalker/notebooks/run_chemwalker.ipynb](2.%20Mass%20Spectrometric%20Analysis/ChemWalker/notebooks/run_chemwalker.ipynb) | Run ChemWalker candidate-ranking analysis |
| [ChemWalker/](2.%20Mass%20Spectrometric%20Analysis/ChemWalker/) | ChemWalker software, notebooks, and environment instructions |

### Inputs and outputs

Inputs include chemical structures, the predicted suspect list, HRMS features, mass spectra, molecular-network data, candidate structures, and verified seed assignments. Outputs include predicted EI spectra, candidate rankings, and structural annotations supporting the dust exposome atlas.

The manuscript also describes retrospective searches of azo-compound MS/MS spectra using **MASST** against public **MassIVE** datasets to explore environmental distribution.

Detailed procedures are provided in **Methods 4.4** and **Supplementary Notes 3–5**.

---

# 🧪 3. Similar Mixture Approach

SMACH assesses whether the chemical composition of a candidate mixture is sufficiently similar to a tested reference mixture to support use of the reference dose–response data.

### Main workflow

```text
Measured concentrations of the selected compounds
        ↓
Concentration fractions for reference and candidate mixtures
        ↓
Toxicity weights derived from toxic equivalency factors
        ↓
Weighted Euclidean distance at the reference BMD
        ↓
Comparison with the study's similarity boundary
```

The study compares the mixtures at the reference **benchmark dose (BMD)**. For candidate mixtures without dose–response data, the reference BMD is used as an approximation. Differences in compound proportions are weighted by their toxic equivalency factors and expressed as a distance.

The similarity boundary is the reference-mixture dose difference associated with an increase in revertant colonies from **two to three times the negative-control level**.

### Main script

[SMACH distance calculation.R](3.%20Similar%20Mixture%20Approach/SMACH%20distance%20calculation.R)

The script calculates SMACH distances using a tab-delimited input file named `MIXN.txt`:

- the `ID` column identifies the `TEF` row, the `NWU` reference row, and candidate samples;
- chemical columns contain the TEFs in the `TEF` row and concentration fractions in the reference and candidate rows; and
- `reference_BMD` specifies the BMD of the tested reference mixture.

The concentration fractions must be prepared from the measured concentrations before using this script, because the script uses the supplied reference and candidate values directly. The BMD and similarity boundary must use the same dose scale.

The output is `4_MIXN_SMACH_distance_results.txt`, containing `ID` and `SMACH_distance` columns. The resulting distances are used for comparison with the study's similarity boundary.

### Reference-mixture experiments

The manuscript reconstructs four class-specific reference mixtures and their combined mixture, **MIX N**, from the measured composition of the **NWU** sample:

| Reference mixture | Number of compounds |
|---|---:|
| Legacy polycyclic aromatic hydrocarbons (PAHs) | 16 |
| Novel polycyclic aromatic compounds (PACs) | 52 |
| Halogenated azo compounds | 14 |
| Nonhalogenated azo analogues | 11 |
| MIX N, combining the four groups | 93 |

Mutagenicity is assessed using a mini-Ames test with **Salmonella typhimurium TA98** and **rat liver S9 metabolic activation**. BMDs are estimated by Bayesian model averaging. The manuscript then uses bioanalytical equivalent concentrations (**BEQbio** and **BEQchem**) to estimate the contributions of reconstructed mixtures to dust mutagenicity.

Detailed procedures and compositions are provided in **Methods 4.6**, **Supplementary Note 7**, and **Supplementary Tables 13 and 14**.

---

# ⚙️ Software Requirements

The repository uses **Python**, **R**, and **Jupyter notebooks**. Key tools and methods named in the manuscript include:

| Analysis | Tools and methods |
|---|---|
| Molecular modeling | Chemical Checker Signaturizer, Boruta, and SMOTE |
| Mass spectrometric processing and annotation | MZmine, NEIMS, MetFrag, GenForm, ChemWalker, and ModiFinder |
| Environmental distribution searches | MASST and GNPS/MassIVE |
| Dose–response analysis | Bayesian BMD model averaging; Prism for dose–response visualization |

Software-specific instructions are provided in the [NEIMS README](2.%20Mass%20Spectrometric%20Analysis/NEIMS/README.md) and [ChemWalker README](2.%20Mass%20Spectrometric%20Analysis/ChemWalker/README.md).

---

# 🚀 Workflow Summary

```text
Data-driven carcinogenicity prediction
        ↓
Potential carcinogens as prior molecular information
        ↓
GC-EI and LC-ESI HRMS screening
        ↓
Dust exposome atlas and chemical prioritization
        ↓
Reference-mixture experiments and dose–response analysis
        ↓
SMACH similarity assessment and mutagenicity attribution
```

Related analyses explore global PAH exposure patterns using SMACH and search public MS/MS datasets using MASST.

---

# 📖 Citation

If you use these workflows, please cite the associated manuscript:

> **Su, W. et al. A risk-based approach prioritizes carcinogens from the exposome atlas of indoor dust.**

### Related methods described in the manuscript

1. **Bertoni, M. et al.** Bioactivity descriptors for uncharacterized chemical compounds. *Nature Communications* **12**, 3932 (2021).
2. **Wei, J. N., Belanger, D., Adams, R. P. & Sculley, D.** Rapid prediction of electron–ionization mass spectrometry using neural networks. *ACS Central Science* **5**, 700–708 (2019).
3. **Ruttkies, C., Schymanski, E. L., Wolf, S., Hollender, J. & Neumann, S.** MetFrag relaunched: incorporating strategies beyond in silico fragmentation. *Journal of Cheminformatics* **8**, 3 (2016).
4. **Marshall, S. et al.** An empirical approach to sufficient similarity: combining exposure data and mixtures toxicology data. *Risk Analysis* **33**, 1582–1595 (2013).
