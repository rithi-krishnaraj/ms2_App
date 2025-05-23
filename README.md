# MS2 Scan Analyzer

A Streamlit-based web application for interactive exploration and analysis of tandem mass spectrometry data in MGF or mzML formats. Use both built-in visualization tools and the GNPS FASTSearch database to compare your spectra against public libraries and download comprehensive results.
---

## Files
- Use ms2_main.py with no multiprocessing for downloadable csv file
- Use ms2_multiprocessing.py for quicker creation of downloadable csv file

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
  - Select any scan for direct FASTSearch queries against GNPS libraries.  
  - Adjustable parameters:  
    - Library selection from common indices ( GNPS, ORNL Bioscales, Massive, MetaRepo, etc.)  
    - Analog search toggle (Yes/No)  
    - Cache usage toggle  
    - Delta mass window (below/above in Da)  
    - Precursor & fragment mass tolerances (Da)  
    - Cosine similarity threshold (0.0–1.0)  

- **GNPS Link Generator**  
  - One‐click generation of a pre‐filled GNPS FASTSearch URL containing your scan’s peak list and parameters.  
  - Facilitates offline or manual FASTSearch exploration.

- **Result Display & Export**  
  - View matched results in an interactive data table: delta mass, USI, charge, cosine score, matching peaks count, dataset, status.  
  - Download complete FASTSearch results across **all** scans in a single CSV via parallel processing (CAUTION: Streamlit does not allow other functionalities to be used while the csv file is being created, and it can take minutes up to hours/days for full creation of the file, depending on if multiprocessing is used or not)  

- **Mirror Plot Comparison**  
  - Side‐by‐side “mirror” plots of the uploaded scan vs. any matching USI spectrum.  
  - Compare **Unfiltered** and **Filtered** views to visually assess peak alignment.

---
## Getting Started

### Requirements

- Python 3.8+  
- [Streamlit](https://streamlit.io/)  
- Key Python packages:
  pip install streamlit pandas numpy matplotlib plotly requests xlsxwriter spectrum-utils pyteomics time json

### Running the App Locally

Follow these steps to download, set up, and launch the MS2 Scan Analyzer on your local computer.

1. Open a terminal and navigate to the folder where you’d like to keep the project, and run the following commands:  
```bash
git clone https://github.com/rithi-krishnaraj/ms2_App.git
cd ms2-scan-analyzer
- For macOS/Linux:
    python3 -m venv .venv
    source .venv/bin/activate
- For Windows:
    python -m venv .venv
    .\.venv\Scripts\Activate.ps1
pip install -r requirements.txt
pip install python-dotenv
streamlit run ms2_main.py (replace ms2_main.py with ms2_multiprocessing for faster csv file creation)

2. Upload an mgf or mzmL file to use the Streamlit App



