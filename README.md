# MS2 Scan Analyzer

A Streamlit-based web application for interactive exploration and analysis of MS² (tandem mass spectrometry) data in MGF or mzML formats. Leverage both built-in visualization tools and the GNPS FASTSearch API to compare your spectra against public libraries and download comprehensive results.
---

## Features

- **File Upload & Parsing**  
  - Upload **MGF** or **mzML** files.  
  - Automatic parsing of scan metadata (scan number, spectrum ID, precursor m/z, charge state, SMILES if available).  

- **Scan Metadata Viewer**  
  - Browse all parsed scans in an expandable data table.  
  - View key metadata: Scan Number, Spectrum ID, Precursor m/z, Charge State, SMILES.

- **Interactive Spectrum Visualization**  
  - **Unfiltered Spectrum**: raw m/z vs intensity plot via Plotly.  
  - **Filtered Spectrum**: remove noise peaks (keep top‐6 peaks in ±25 Da window, exclude precursor), with two normalization options:  
    - _Normalized_ (vector‐norm scaling)  
    - _Square‐root normalized_ (enhanced low‐intensity peak visibility).

- **GNPS FASTSearch Integration**  
  - Select any scan for on‐the‐fly FASTSearch queries against GNPS libraries.  
  - Adjustable parameters:  
    - Library selection from common indices (e.g. GNPS, ORNL Bioscales, Massive, MetaRepo, etc.)  
    - Analog search toggle (Yes/No)  
    - Cache usage toggle  
    - Delta mass window (below/above in Da)  
    - Precursor & fragment mass tolerances (Da)  
    - Cosine similarity threshold (0.0–1.0)  

- **GNPS Link Generator**  
  - One‐click generation of a pre‐filled GNPS FASTSearch URL containing your scan’s peak list and parameters.  
  - Facilitates offline or manual FASTSearch exploration.

- **Result Display & Export**  
  - View matched results in an interactive data table: Δ mass, USI, charge, cosine score, matching peaks count, dataset, status.  
  - Download complete FASTSearch results across **all** scans in a single CSV via parallel processing.  

- **Mirror Plot Comparison**  
  - Side‐by‐side “mirror” plots of your uploaded scan vs. any matched USI spectrum.  
  - Compare **Unfiltered** and **Filtered** views to visually assess peak alignment.

---

## Getting Started

### Prerequisites

- Python 3.8+  
- [Streamlit](https://streamlit.io/)  
- Key Python packages:
  ```bash
  pip install streamlit pandas numpy matplotlib plotly requests xlsxwriter spectrum-utils pyteomics
