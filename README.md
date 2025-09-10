# MS2 Scan Analyzer

This project is a Streamlit application designed for analyzing MS2 scans. It provides a user-friendly interface for uploading mass spectrometry data files, processing the data, and visualizing the results.

## Project Structure

```

ms2_main.py          # Main logic for the Streamlit application
requirements.txt          # List of dependencies
requirements-windows.txt
.github
└── workflows
   └── build_windows.yaml # GitHub Actions workflow for building the app
README.md                 # Project documentation
```

## Installation

To install the required dependencies, run the following command:

```
pip install -r requirements.txt
```

## Running the Application

To run the Streamlit application locally, execute the following command:

```
streamlit run ms2_main.py
```

## Building the Windows Installable Application

To create a Windows installable application, the project includes a GitHub Actions workflow. This workflow is defined in `.github/workflows/build_windows.yml`. It automates the process of setting up the environment, installing dependencies, and building the executable. Download the MS2-Scan-Analyzer.msi file from the latest release, and click MS2-Scan-Analyzer in your downloads to run the app. Bypass warnings and start download by working through steps and choosing a location for the download. Once the installable application is downloaded, click the folder in your selected computer path named MS2-Scan-Analyzer and run the MS2-Scan-Analyzer file to run the Streamlit App. 

## Usage
1. Upload your mass spectrometry data files in either MGF or mzML format.
2. View All Scan Numbers and their MetaData
3. Select a scan number to view the MS2 Spectrum.
4. View Unfiltered and Filtered (Normalized & Square-Root) Spectrums for selected user scan
5. Adjust parameters for GNPS FastSearch as needed
6. Generate GNPS Link to use GNPS FastSearch Dashboard directly with your selected user can
7. View Matching USI Results 
8. Use Metabolomics Resolver Spectrum Viewer to compare selected user scan with a selected matching USI in both unfiltered and filtered mirror plots
8. Generate results and download them in CSV format.


