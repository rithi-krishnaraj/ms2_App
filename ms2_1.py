import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 
import spectrum_utils.plot as sup
import spectrum_utils.spectrum as sus
import plotly.tools as tls
import json
import requests
import io
import time
import urllib.parse
from urllib.parse import urlparse, parse_qs, unquote, urlencode
from xlsxwriter import Workbook
from multiprocessing import Pool, cpu_count

LIBRARIES = ["gnpsdata_index", "ORNL_Bioscales2", "ORNL_Populus_LC_MSMS", "gnpsdata_test_index", "gnpslibrary", "massivedata_index", "massivekb_index", "metabolomicspanrepo_index_latest", "metabolomicspanrepo_index_nightly", "panrepo_2024_11_12"]

def create_mirror_plot(spectrum_top, spectrum_bottom):
    try:
        # Create the mirror plot
        fig, ax = plt.subplots(figsize=(12, 6))
        sup.mirror(spectrum_top, spectrum_bottom, ax=ax)
        ax.set_title("Top: User Scan, Bottom: Selected USI", fontsize=16)

        # Display the plot in Streamlit
        st.pyplot(fig)
        

    except Exception as e:
        st.error(f"Failed to create mirror plot: {e}")

def generate_webAPI(peaks, precursor_mz, charge, library_select, analog_select, cache_select, delta_mass_below, delta_mass_above, pm_tolerance, fragment_tolerance, cosine_threshold):
    base_url = "https://fasst.gnps2.org/search"
    
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

def process_scan_library(values):
    scan, library = values
    try:
        api_response = generate_webAPI(
            peaks=scan["peaks"],
            precursor_mz=scan["PEPMASS Number"],
            charge=scan["Charge State"],
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

    with Pool(processes=cpu_count()) as pool:
        results = pool.map(process_scan_library, tasks)

    for result in results:
        all_results.extend(result)

    return all_results

if __name__ == "__main__":
    st.title("MS2 Scan Analyzer")

    mgf_file = st.file_uploader("Choose a file", type="mgf")
    
    if mgf_file is not None:
        scan_numbers, scan_metadata = read_file(mgf_file)
        if not scan_metadata:
            st.error("No scans found in the uploaded MGF file.")
        else:

            if st.button("Generate GNPS Results for all scans in CSV File"):
                # Start the timer
                start_time = time.time()

                # Use st.spinner only in the main thread
                with st.spinner("Processing scans..."):
                    matching_results = process_scans_in_parallel(scan_metadata)

                # Stop the timer
                end_time = time.time()
                elapsed_time = end_time - start_time

                if matching_results:
                    results_df = pd.DataFrame(matching_results)
                    desired_columns = ["Scan Number", "Library", "Delta Mass", "USI", "Charge", "Cosine", "Matching Peaks", "Dataset", "Status"]
                    
                    # Add missing columns (if any) with default NaN values.
                    missing_cols = set(desired_columns) - set(results_df.columns)
                    for col in missing_cols:
                        results_df[col] = np.nan
                    
                    # Reorder the DataFrame columns.
                    results_df = results_df[desired_columns]
                    
                    # Convert the DataFrame to CSV format encoded in UTF-8.
                    csv_data = results_df.to_csv(index=False).encode('utf-8')
                    
                    st.download_button(
                        label="Download GNPS Results in CSV File",
                        data=csv_data,
                        file_name="gnps_results.csv",
                        mime="text/csv"
                    )
                    
                    st.success(f"File created successfully in {elapsed_time:.2f} seconds!")

            df = pd.DataFrame(scan_metadata, columns=["Scan Number", "Spectrum ID", "PEPMASS Number", "Charge State", "SMILES ID"])
            with st.expander("Show Scan Numbers and Metadata"):
                st.dataframe(df)

            scan_input = st.selectbox("Select Scan Number to view MS2 Spectrum", options=scan_numbers)

            if scan_input:
                user_scan = next((scan for scan in scan_metadata if scan["Scan Number"] == scan_input), None)
                mz_filtered, sqrt_filtered, normal_filtered = peak_filtering(user_scan)
                
                with st.expander("View Spectrum", expanded = False):
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
                
                with st.expander(f"GNPS FASTSearch for Scan {user_scan['Scan Number']}"):
                    precursor_mz = user_scan["PEPMASS Number"]
                    charge = user_scan["Charge State"]
                    col1, col2, col3 = st.columns(3)
                    
                    with col1:
                        #peak_select = st.selectbox("Select Peaks", ["Normalized Peaks", "Unfiltered Peaks", "Square Root Peaks"], key="peak_select")
                        library_select = st.selectbox("Select Library", ["gnpsdata_index", "ORNL_Bioscales2", "ORNL_Populus_LC_MSMS", "gnpsdata_test_index", "gnpslibrary", "massivedata_index", "massivekb_index", "metabolomicspanrepo_index_latest", "metabolomicspanrepo_index_nightly", "panrepo_2024_11_12"], key="library_select")
                        analog_select = st.selectbox("Analog Search", ["No", "Yes"], key="analog_select")
                    with col2:
                        delta_mass_below = st.number_input("Delta Mass Below (Da)", value = 130, key="delta_mass_below")
                        delta_mass_above = st.number_input("Delta Mass Above (Da)", value = 200, key="delta_mass_above")
                        cache_select = st.selectbox("Use Cache", ["Yes", "No"], key="cache_select")
                    with col3:
                        pm_tolerance = st.number_input("PM Tolerance (Da)", min_value=0.0, value=0.05, step = 0.05, key="pm_tolerance")
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

                    # Display results if they exist (either from current run or previous run)
                    if api_response and "results" in api_response:
        
                        try: 
                            st.subheader("Matching Results")
                            results_df = pd.DataFrame(api_response["results"])
                            desired_columns = ["Delta Mass", "USI", "Charge", "Cosine", "Matching Peaks", "Dataset", "Status"]
                            results_df = results_df[desired_columns]
                            st.dataframe(results_df)
                            
                            #USI selection for Spectrum Exploration
                            st.subheader("Metabolomics Resolver Spectrum Viewer")
                            usi_list = results_df["USI"].tolist()
                            selected_usi = st.selectbox("Select a USI for Spectrum Exploration", usi_list,  key="usi_selector")
                            
                            if selected_usi:
                                spectrum_tabs = st.tabs(["Unfiltered Spectrum", "Filtered Spectrum"])
                                with spectrum_tabs[0]:
                                    st.write("Mirror Plot: Unfiltered Peaks")
                                    try:
                                        # Convert the Scan Number to a string for the identifier
                                        spectrum_top = sus.MsmsSpectrum(
                                            mz=user_scan["m/z data"],
                                            intensity=user_scan["intensity data"],
                                            identifier=str(user_scan["Scan Number"]),
                                            precursor_mz=user_scan["PEPMASS Number"],
                                            precursor_charge=user_scan["Charge State"]
                                        )

                                        # Fetch the spectrum for the selected USI
                                        spectrum_bottom = sus.MsmsSpectrum.from_usi(
                                            selected_usi,
                                            precursor_mz=user_scan["PEPMASS Number"],
                                            precursor_charge=user_scan["Charge State"]
                                        )
                                        create_mirror_plot(spectrum_top, spectrum_bottom)
                                    except Exception as e:
                                        st.error(f"Failed to create unfiltered mirror plot: {e}")

                                # Filtered Spectrum Tab
                                with spectrum_tabs[1]:
                                    st.write("Mirror Plot: Filtered Peaks")
                                    try:
                                        # Create the filtered top spectrum
                                        spectrum_top = sus.MsmsSpectrum(
                                            mz=mz_filtered,
                                            intensity=normal_filtered,
                                            identifier=str(user_scan["Scan Number"]),
                                            precursor_mz=user_scan["PEPMASS Number"],
                                            precursor_charge=user_scan["Charge State"]
                                        )

                                        # Fetch the spectrum for the selected USI
                                        spectrum_bottom = sus.MsmsSpectrum.from_usi(
                                            selected_usi,
                                            precursor_mz=user_scan["PEPMASS Number"],
                                            precursor_charge=user_scan["Charge State"]
                                        )

                                        # Create a temporary dictionary for the USI spectrum
                                        usi_scan = {
                                            "m/z data": spectrum_bottom.mz,
                                            "intensity data": spectrum_bottom.intensity,
                                            "Scan Number": spectrum_bottom.identifier,
                                            "PEPMASS Number": spectrum_bottom.precursor_mz,
                                            "Charge State": spectrum_bottom.precursor_charge
                                        }

                                        # Filter the USI spectrum using the peak_filtering function
                                        usi_mz_filtered, usi_sqrt_filtered, usi_normal_filtered = peak_filtering(usi_scan)

                                        # Create the filtered bottom spectrum
                                        spectrum_bottom = sus.MsmsSpectrum(
                                            mz=usi_mz_filtered,
                                            intensity=usi_normal_filtered,
                                            identifier=str(usi_scan["Scan Number"]),
                                            precursor_mz=usi_scan["PEPMASS Number"],
                                            precursor_charge=usi_scan["Charge State"]
                                        )

                                        # Create the mirror plot
                                        create_mirror_plot(spectrum_top, spectrum_bottom)
                                    except Exception as e:
                                        st.error(f"Failed to create filtered mirror plot: {e}")
                        except KeyError as e:
                            st.error(f"Could not retrieve spectrum data for the selected USI: {e}")
