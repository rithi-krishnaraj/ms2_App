# Streamlit App Builder

This project is a Streamlit application designed for analyzing MS2 scans. It provides a user-friendly interface for uploading mass spectrometry data files, processing the data, and visualizing the results.

## Project Structure

```
streamlit-app-builder
├── src
│   ├── ms2_main.py          # Main logic for the Streamlit application
│   └── assets
│       └── app_icon.ico     # Application icon for Windows
├── requirements.txt          # List of dependencies
├── setup.py                  # Packaging configuration
├── .github
│   └── workflows
│       └── build_windows.yml # GitHub Actions workflow for building the app
└── README.md                 # Project documentation
```

## Installation

To install the required dependencies, run the following command:

```
pip install -r requirements.txt
```

## Running the Application

To run the Streamlit application, execute the following command:

```
streamlit run src/ms2_main.py
```

## Building the Application

To create a Windows installable application, the project includes a GitHub Actions workflow. This workflow is defined in `.github/workflows/build_windows.yml`. It automates the process of setting up the environment, installing dependencies, and building the executable.

## Usage

1. Upload your mass spectrometry data files in either MGF or mzML format.
2. Select the scan number to view the MS2 spectrum.
3. Adjust the parameters for GNPS FASTSearch as needed.
4. Generate results and download them in CSV format.


