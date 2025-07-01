import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 
import spectrum_utils.plot as sup
import spectrum_utils.spectrum as sus
import plotly.tools as tls
import json
import requests
import time
import urllib.parse
from xlsxwriter import Workbook
from multiprocessing import Pool, cpu_count
from pyteomics import mzml
LIBRARIES = ["gnpsdata_index", "ORNL_Bioscales2", "ORNL_Populus_LC_MSMS", "gnpsdata_test_index", "gnpslibrary", "massivedata_index", "massivekb_index", "metabolomicspanrepo_index_latest", "metabolomicspanrepo_index_nightly", "panrepo_2024_11_12"]

def create_mirror_plot(spectrum_top, spectrum_bottom):
    try:
        fig, ax = plt.subplots(figsize=(12, 6))
        sup.mirror(spectrum_top, spectrum_bottom, ax=ax)
        ax.set_title("Top: User Scan, Bottom: Selected USI", fontsize=16)
        st.pyplot(fig)
    except Exception as e:
        st.error(f"Failed to create mirror plot: {e}")

def generate_webAPI(peaks, precursor_mz, charge, library_select, analog_select, cache_select, delta_mass_below, delta_mass_above, pm_tolerance, fragment_tolerance, cosine_threshold):
    base_url = "https://fasst.gnps2.org/search"
    
    try:
        charge = int(charge)
    except (ValueError, TypeError):
        charge = 0

    query_spectrum = {
        "n_peaks": len(peaks),
        "peaks": peaks,
        "precursor_charge": charge,
        "precursor_mz": precursor_mz
    }

    query_spectrum_json = json.dumps(query_spectrum)
    
    payload = {
        "library": library_select,
        "analog": "Yes" if analog_select == "Yes" else "No",
        "cache": "Yes" if cache_select == "Yes" else "No", 
        "lower_delta": delta_mass_below,
        "upper_delta": delta_mass_above,
        "pm_tolerance": pm_tolerance,
        "fragment_tolerance": fragment_tolerance,
        "cosine_threshold": cosine_threshold,
        "query_spectrum": query_spectrum_json,
    }
    
    try:
        response = requests.post(base_url, data=payload)
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as e:
        st.error(f"API request failed: {e}")
        return None

def gnps_input_link(mz, intensity, precursor_mz, charge, library_select, analog_select, cache_select, delta_mass_below, delta_mass_above, pm_tolerance, fragment_tolerance, cosine_threshold):
    base_url = "http://fasst.gnps2.org/fastsearch/?"

    data = {
        "usi1": None,
        "precursor_mz": precursor_mz,
        "charge": charge,
        "library_select": library_select,
        "analog_select": "Yes" if analog_select == "Yes" else "No", 
        "delta_mass_below": delta_mass_below,
        "delta_mass_above": delta_mass_above,
        "pm_tolerance": pm_tolerance,
        "fragment_tolerance": fragment_tolerance,
        "cosine_threshold": cosine_threshold,
        "use_peaks": 1
    }
    
    peak_joining = "%5Cn".join(f"{mz}%5Ct{intensity}" for mz, intensity in zip(mz, intensity))
    
    query_string = "&".join(f"{key}={urllib.parse.quote_plus(str(value))}" for key, value in data.items())
    final_url = f"{base_url}{query_string}#%7B%22peaks%22%3A%20%22{peak_joining}%22%7D"

    return final_url

@st.cache_resource
def peak_filtering(user_scan):
    mz_array = user_scan["m/z data"]
    intensity_array = user_scan["intensity data"]
    filtered_mz = []
    filtered_intensities = []
    pepmass_val = float(user_scan["PEPMASS Number"])

    for i, mz in enumerate(mz_array):
        peak_range = [j for j in range(len(mz_array)) if abs(mz_array[j] - mz) <= 25]
        sorted_range = sorted(peak_range, key=lambda j: intensity_array[j], reverse=True)
        if i in sorted_range[:6]:
            if abs(mz - pepmass_val) > 17:
                filtered_mz.append(mz)
                filtered_intensities.append(intensity_array[i])

    sqrt_data, normalized_data = peak_normalizing(filtered_intensities)
    return filtered_mz, sqrt_data, normalized_data

@st.cache_resource 
def peak_normalizing(filtered_intensities):
    normalized_intensities = np.copy(filtered_intensities) / np.linalg.norm(filtered_intensities)
    sqrt_intensities = np.sqrt(normalized_intensities)
    return sqrt_intensities, normalized_intensities

@st.cache_data
def peak_visual(mzs, intensities, scanNum, pepmass, charge):
    spectrum = sus.MsmsSpectrum(mz=mzs, intensity=intensities, identifier=scanNum, precursor_mz=pepmass, precursor_charge=charge)
    sup.spectrum(spectrum)
    plt.title(f"MS2 Spectrum for Scan {scanNum}")
    plt.xlabel("m/z", fontsize=11)
    plt.ylabel("Intensity", fontsize=11)
    fig = plt.gcf()
    plotly_fig = tls.mpl_to_plotly(fig)
    plotly_fig.update_traces(hoverinfo="x+y")
    st.plotly_chart(plotly_fig)

@st.cache_resource
def read_file(mgf_file):
    scans = [] 
    current_scan = None
    scan_numbers = []

    file = mgf_file.read().decode('utf-8').splitlines()
    for line in file:
        line = line.strip()
        if line == "BEGIN IONS":
            current_scan = {"Scan Number": 0, "Spectrum ID": '', "PEPMASS Number": 0.0, "Charge State": 0, "SMILES ID": '', "peaks": [], "m/z data": [], "intensity data": []}
        elif line == "END IONS":
            if current_scan:
                scans.append(current_scan)
                current_scan = None
        elif current_scan is not None:
            if "=" in line:
                data = line.split('=', 1)
                if len(data) == 2:
                    key, value = data
                    if key == "SCANS":
                        current_scan["Scan Number"] = int(value)
                        scan_numbers.append(int(value))
                    elif key == "SPECTRUMID":
                        current_scan['Spectrum ID'] = str(value)
                    elif key == "PEPMASS":
                        current_scan["PEPMASS Number"] = float(value)
                    elif key == "CHARGE":
                        current_scan["Charge State"] = int(value)
                    elif key == 'SMILES':
                        current_scan["SMILES ID"] = str(value.strip())
                    else:
                        continue
            else:
                try:
                    data2 = line.split()
                    if len(data2) == 2:
                        mz, intensity = data2
                        current_scan["peaks"].append((float(mz), float(intensity)))
                        current_scan["m/z data"].append(float(mz))
                        current_scan["intensity data"].append(float(intensity))
                except ValueError:
                    print(f"Skipping unreadable data in line: '{line}")
                    continue

    return scan_numbers, scans

@st.cache_resource
def read_mzml_file(mzml_file):
    try:
        scans = []
        scan_numbers = []

        with mzml.read(mzml_file) as reader:
            for i, spectrum in enumerate(reader):
                if spectrum.get('ms level') == 1 or spectrum.get('ms level') == 2:
                    scan_id = spectrum.get('id', 'Unknown')
                    scan_number = i + 1 #Default indexing 
                    if 'scan=' in scan_id:
                        scan_part = scan_id.split('scan=')[-1].split()[0]
                        try:
                            scan_number = int(scan_part)
                        except ValueError:
                            pass

                    # Precursor information handling
                    precursor_mz = 0.0
                    charge = 0
                    if 'precursorList' in spectrum:
                        precursor = spectrum['precursorList']['precursor'][0]
                        selected_ion = precursor.get('selectedIonList', {}).get('selectedIon', [{}])[0]
                        precursor_mz = selected_ion.get('selected ion m/z', 0.0)
                        charge = int(selected_ion.get('charge state', 0))

                    # Convert numpy arrays to lists
                    mz_array = spectrum['m/z array'].tolist()
                    intensity_array = spectrum['intensity array'].tolist()

                    scans.append({
                        "Scan Number": scan_number,
                        "Spectrum ID": scan_id,
                        "PEPMASS Number": precursor_mz,
                        "Charge State": charge,
                        "SMILES ID": '',
                        "peaks": list(zip(mz_array, intensity_array)),
                        "m/z data": mz_array,
                        "intensity data": intensity_array,
                    })
                    scan_numbers.append(scan_number)

        return scan_numbers, scans
    except Exception as e:
        st.error(f"Failed to read mzML file: {str(e)}")
        return [], []
    


def process_scan_library(values):
    scan, library = values
    try:
        try:
            charge = int(scan["Charge State"])
        except (ValueError, KeyError):
            charge = 0

        api_response = generate_webAPI(
            peaks=scan["peaks"],
            precursor_mz=scan["PEPMASS Number"],
            charge=charge,
            library_select=library,
            analog_select="No",
            cache_select="Yes",
            delta_mass_below=130,
            delta_mass_above=200,
            pm_tolerance=0.05,
            fragment_tolerance=0.05,
            cosine_threshold=0.7
        )
        if api_response and "results" in api_response:
            for result in api_response["results"]:
                result["Scan Number"] = scan["Scan Number"]
                result["Library"] = library
            return api_response["results"]
    except Exception as e:
        print(f"Error processing scan {scan['Scan Number']} with library {library}: {e}")
    return []

def process_scans_in_parallel(scans):
    all_results = []
    tasks = [(scan, library) for scan in scans for library in LIBRARIES]

    with Pool(processes=1) as pool:
        results = pool.map(process_scan_library, tasks)

    for result in results:
        all_results.extend(result)

    return all_results

if __name__ == "__main__":
    st.title("MS2 Scan Analyzer")

    uploaded_file = st.file_uploader("Choose a file", type=["mgf", "mzML"])
    
    if uploaded_file is not None:
        file_type = uploaded_file.name.split(".")[-1].lower()
        if file_type == "mgf":
            scan_numbers, scan_metadata = read_file(uploaded_file)
        elif file_type == "mzml":
            scan_numbers, scan_metadata = read_mzml_file(uploaded_file)
        else:
            st.error("Unsupported file type. Please upload an MGF or mzML file.")
            scan_numbers, scan_metadata = [], []
        
        # Convert all scan numbers to strings for UI consistency
        scan_numbers = [str(sn) for sn in scan_numbers]
        for scan in scan_metadata:
            scan["Scan Number"] = str(scan["Scan Number"])

        if not scan_metadata:
            st.error("No scans found in the uploaded file.")
        else:
            st.subheader(f"Analyzing {file_type.upper()} File")

            if st.button("Generate GNPS Results for all scans in CSV File"):
                start_time = time.time()
                with st.spinner("Processing scans..."):
                    matching_results = process_scans_in_parallel(scan_metadata)
                end_time = time.time()

                if matching_results:
                    results_df = pd.DataFrame(matching_results)
                    desired_columns = ["Scan Number", "Library", "Delta Mass", "USI", "Charge", "Cosine", "Matching Peaks", "Dataset", "Status"]
                    
                    missing_cols = set(desired_columns) - set(results_df.columns)
                    for col in missing_cols:
                        results_df[col] = np.nan
                    
                    results_df = results_df[desired_columns]
                    csv_data = results_df.to_csv(index=False).encode('utf-8')
                    
                    st.download_button(
                        label="Download GNPS Results in CSV File",
                        data=csv_data,
                        file_name="gnps_results.csv",
                        mime="text/csv"
                    )
                    st.success(f"File created successfully in {end_time-start_time:.2f} seconds!")

            df = pd.DataFrame(scan_metadata, columns=["Scan Number", "Spectrum ID", "PEPMASS Number", "Charge State", "SMILES ID"])
            with st.expander("Show Scan Numbers and Metadata"):
                st.dataframe(df)

            scan_input = st.selectbox("Select Scan Number to view MS2 Spectrum", options=scan_numbers)

            if scan_input:
                user_scan = next((scan for scan in scan_metadata if scan["Scan Number"] == scan_input), None)
                mz_filtered, sqrt_filtered, normal_filtered = peak_filtering(user_scan)
                
                with st.expander("View Spectrum", expanded=False):
                    spectrum_option = st.tabs(["Unfiltered Spectrum", "Filtered Spectrum - Normalized", "Filtered Spectrum - Square Root Normalized"])
                    with spectrum_option[0]:
                        peak_visual(user_scan["m/z data"], user_scan["intensity data"],
                                    str(user_scan["Scan Number"]), user_scan["PEPMASS Number"], user_scan["Charge State"])
                    with spectrum_option[1]:
                        peak_visual(mz_filtered, normal_filtered,
                                    str(user_scan["Scan Number"]), user_scan["PEPMASS Number"], user_scan["Charge State"])
                    with spectrum_option[2]:
                        peak_visual(mz_filtered, sqrt_filtered,
                                    str(user_scan["Scan Number"]), user_scan["PEPMASS Number"], user_scan["Charge State"])
                
                st.subheader(f"GNPS FASTSearch for Scan {user_scan['Scan Number']}")
                precursor_mz = user_scan["PEPMASS Number"]
                charge = user_scan["Charge State"]
                col1, col2, col3 = st.columns(3)
                
                with col1:
                    library_select = st.selectbox("Select Library", LIBRARIES, key="library_select")
                    analog_select = st.selectbox("Analog Search", ["No", "Yes"], key="analog_select")
                with col2:
                    delta_mass_below = st.number_input("Delta Mass Below (Da)", value=130, key="delta_mass_below")
                    delta_mass_above = st.number_input("Delta Mass Above (Da)", value=200, key="delta_mass_above")
                    cache_select = st.selectbox("Use Cache", ["Yes", "No"], key="cache_select")
                with col3:
                    pm_tolerance = st.number_input("PM Tolerance (Da)", min_value=0.0, value=0.05, step=0.05, key="pm_tolerance")
                    fragment_tolerance = st.number_input("Fragment Mass Tolerance (Da)", min_value=0.0, value=0.05, step=0.05, key="fragment_tolerance")
                    cosine_threshold = st.number_input("Cosine Similarity Threshold", min_value=0.0, max_value=1.0, value=0.7, step=0.05, key="cosine_threshold")
                
                peaks = user_scan["peaks"]
                mz = user_scan["m/z data"]
                intensity = user_scan["intensity data"]
                api_response = generate_webAPI(peaks, precursor_mz, charge, library_select, analog_select, cache_select, delta_mass_below, delta_mass_above, pm_tolerance, fragment_tolerance, cosine_threshold)
                

                if st.button("Generate GNPS FASTSearch Link"):
                    gnpsLink = gnps_input_link(mz, intensity, precursor_mz, charge, library_select, analog_select, cache_select, delta_mass_below, delta_mass_above, pm_tolerance, fragment_tolerance, cosine_threshold)
                    st.write("GNPS Input Link Created with Scan Data!")
                    st.markdown(f"[Go to GNPS Link]({gnpsLink})", unsafe_allow_html=True)

                if api_response and "results" in api_response:
                    try: 
                        st.subheader("Matching Results")
                        results_df = pd.DataFrame(api_response["results"])
                        desired_columns = ["Delta Mass", "USI", "Charge", "Cosine", "Matching Peaks", "Dataset", "Status"]
                        results_df = results_df[desired_columns]
                        st.dataframe(results_df)
                    except KeyError as e:
                        st.error(f"No Matching Results Found")
                    try:
                        st.subheader("Metabolomics Resolver Spectrum Viewer")
                        usi_list = results_df["USI"].tolist()
                        selected_usi = st.selectbox("Select a USI for Spectrum Exploration", usi_list, key="usi_selector")
                        
                        if selected_usi:
                            spectrum_tabs = st.tabs(["Unfiltered Spectrum", "Filtered Spectrum"])
                            with spectrum_tabs[0]:
                                try:
                                    spectrum_top = sus.MsmsSpectrum(
                                        mz=user_scan["m/z data"],
                                        intensity=user_scan["intensity data"],
                                        identifier=str(user_scan["Scan Number"]),
                                        precursor_mz=user_scan["PEPMASS Number"],
                                        precursor_charge=user_scan["Charge State"]
                                    )
                                    spectrum_bottom = sus.MsmsSpectrum.from_usi(
                                        selected_usi,
                                        precursor_mz=user_scan["PEPMASS Number"],
                                        precursor_charge=user_scan["Charge State"]
                                    )
                                    
                                    create_mirror_plot(spectrum_top, spectrum_bottom)
                                except Exception as e:
                                    st.error(f"Failed to create unfiltered mirror plot: {e}")

                            with spectrum_tabs[1]:
                                try:
                                    spectrum_top = sus.MsmsSpectrum(
                                        mz=mz_filtered,
                                        intensity=normal_filtered,
                                        identifier=str(user_scan["Scan Number"]),
                                        precursor_mz=user_scan["PEPMASS Number"],
                                        precursor_charge=user_scan["Charge State"]
                                    )

                                    spectrum_bottom = sus.MsmsSpectrum.from_usi(
                                        selected_usi,
                                        precursor_mz = user_scan["PEPMASS Number"],
                                        precursor_charge = user_scan["Charge State"])

                                    usi_scan = {
                                        "m/z data": spectrum_bottom.mz,
                                        "intensity data": spectrum_bottom.intensity,
                                        "Scan Number": spectrum_bottom.identifier,
                                        "PEPMASS Number": spectrum_bottom.precursor_mz,
                                        "Charge State": spectrum_bottom.precursor_charge
                                    }
                                    usi_mz_filtered, _, usi_normal_filtered = peak_filtering(usi_scan)
                                    spectrum_bottom = sus.MsmsSpectrum(
                                        mz=usi_mz_filtered,
                                        intensity=usi_normal_filtered,
                                        identifier=str(usi_scan["Scan Number"]),
                                        precursor_mz=usi_scan["PEPMASS Number"],
                                        precursor_charge=usi_scan["Charge State"]
                                    )
                                    create_mirror_plot(spectrum_top, spectrum_bottom)
                                except Exception as e:
                                    st.error(f"Failed to create filtered mirror plot: {e}")
                    except KeyError as e:
                        st.error(f"Could not retrieve spectrum data: {e}")




