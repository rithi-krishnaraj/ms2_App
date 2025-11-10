import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 
import spectrum_utils.plot as sup
import spectrum_utils.spectrum as sus
import plotly.tools as tls
from pyteomics import mgf, mzml
import requests
import json
import time
from tqdm import tqdm

#To read libraries from GNPS
@st.cache_data
def get_databases(host="https://fasst.gnps2.org"): #DONE
    url = f"{host}/libraries"
    try:
        response = requests.get(url, timeout=30)
        response.raise_for_status()
        libraries_data = response.json()
        library_names = [lib['value'] for lib in libraries_data]
        return library_names
    except Exception as e:
        st.error(f"Failed to fetch libraries: {e}")
        return []

#For File MetaData
@st.cache_resource(show_spinner=False)
def read_mgf_file(mgf_file):
    scans = [] 
    current_scan = None
    scan_numbers = []

    file = mgf_file.read().decode('utf-8').splitlines()
    for line in file:
        line = line.strip()
        if line == "BEGIN IONS":
            current_scan = {"Scan Number": 0, "Spectrum ID": '', "Precursor M/Z": 0.0, "Charge State": 0, "SMILES ID": '', "peaks": [], "m/z data": [], "intensity data": []}
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
                        try:
                            current_scan["Scan Number"] = int(value)
                            scan_numbers.append(int(value))
                        except ValueError:
                            st.warning(f"Cannot parse SCANS value: {value}")
                    elif key == "SPECTRUMID":
                        current_scan['Spectrum ID'] = str(value)
                    elif key == "PEPMASS":
                        try:
                            current_scan["Precursor M/Z"] = float(value)
                        except ValueError:
                            st.warning(f"Cannot parse PEPMASS value: {value}")

                    elif key == "CHARGE":
                        try:
                            current_scan["Charge State"] = int(value)
                        except ValueError:
                            st.warning(f"Cannot parse CHARGE value: {value}")
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
                    st.warning(f"Cannot parse peak line: {line}")
                    continue

    return scan_numbers, scans

@st.cache_resource(show_spinner=False)
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

                    precursor_mz = 0.0
                    charge = 0
                    if 'precursorList' in spectrum:
                        precursor = spectrum['precursorList']['precursor'][0]
                        selected_ion = precursor.get('selectedIonList', {}).get('selectedIon', [{}])[0]
                        precursor_mz = selected_ion.get('selected ion m/z', 0.0)
                        try:
                            charge = int(selected_ion.get('charge state', 0))
                        except Exception:
                            charge = 0

                    mz_array = spectrum['m/z array'].tolist()
                    intensity_array = spectrum['intensity array'].tolist()

                    scans.append({
                        "Scan Number": scan_number,
                        "Spectrum ID": scan_id,
                        "Precursor M/Z": precursor_mz,
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

#For Spectrum Visualization 
@st.cache_data(show_spinner=False)
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

@st.cache_resource(show_spinner=False)
def peak_filtering(user_scan):
    mz_array = user_scan["m/z data"]
    intensity_array = user_scan["intensity data"]
    filtered_mz = []
    filtered_intensities = []
    pepmass_val = float(user_scan["Precursor M/Z"])

    for i, mz in enumerate(mz_array):
        peak_range = [j for j in range(len(mz_array)) if abs(mz_array[j] - mz) <= 25]
        sorted_range = sorted(peak_range, key=lambda j: intensity_array[j], reverse=True)
        if i in sorted_range[:6]:
            if abs(mz - pepmass_val) > 17:
                filtered_mz.append(mz)
                filtered_intensities.append(intensity_array[i])

    sqrt_data, normalized_data = peak_normalizing(filtered_intensities)
    return filtered_mz, sqrt_data, normalized_data

@st.cache_resource(show_spinner=False)
def peak_normalizing(filtered_intensities):
    normalized_intensities = np.copy(filtered_intensities) / np.linalg.norm(filtered_intensities)
    sqrt_intensities = np.sqrt(normalized_intensities)
    return sqrt_intensities, normalized_intensities

def create_mirror_plot(spectrum_top, spectrum_bottom): #DONE
    try:
        fig, ax = plt.subplots(figsize=(12, 6))
        sup.mirror(spectrum_top, spectrum_bottom, ax=ax)
        ax.set_title(f"Top: Scan {spectrum_top.identifier}, Bottom: Selected USI", fontsize=16)
        st.pyplot(fig)
    except Exception as e:
        st.error(f"Failed to create mirror plot: {e}")

#For MASST Querying
@st.cache_data(show_spinner=False)
def execute_all_queries_batch(queries):
    output_results_list = []

    # group queries list into batches
    batch_size = 50
    query_batches = [queries[i:i + batch_size] for i in range(0, len(queries), batch_size)]

    for batch in query_batches:
        status_results_list = []
        for query in tqdm(batch):
            try:
                if "mgf" in query["datasource"] or "mzml" in query["datasource"]:
                    status = query_fasst_api_peaks(
                        query["precursor_mz"],
                        query["charge"],
                        query["peaks"],
                        query["database"],
                        analog=query["analog"],
                        precursor_mz_tol=query["precursor_mz_tol"],
                        fragment_mz_tol=query["fragment_mz_tol"],
                        min_cos=query["min_cos"],
                        blocking=False
                    )
                else:
                    status = {"status": "ERROR", "error": "Unsupported datasource"}

                status_results_list.append(status)
            except KeyboardInterrupt:
                raise
            except Exception as e:
                status_results_list.append({"status": "ERROR", "error": str(e)})
                continue

        time.sleep(1)
        for k in range(0, 60):
            pending_count = 0
            for i, status in enumerate(status_results_list):
                if status.get("status") == "DONE" or status.get("status") == "TIMEOUT" or status.get("status") == "ERROR":
                    continue

                try:
                    results = get_results(status, blocking=False)
                    if results == "PENDING":
                        status["status"] = "PENDING"
                        pending_count += 1
                        continue
                    else:
                        status["status"] = "DONE"
                        status["results"] = results.get("results", [])
                except KeyboardInterrupt:
                    raise
                except Exception as e:
                    status["status"] = "ERROR"
                    continue

            if pending_count == 0:
                break

            # waiting
            time.sleep(5)
        
        # Formatting all results
        for status, query in zip(status_results_list, batch):
            if status.get("status") == "DONE":
                results_df = pd.DataFrame(status.get("results", []))
                # Always add Scan Number and Library from the query
                results_df["Scan Number"] = query["query_scan"]
                results_df["Library"] = query["database"]
                output_results_list.append(results_df)

    if len(output_results_list) == 0:
        st.info("No results found")
        return pd.DataFrame()

    output_results_df = pd.concat(output_results_list)

    return output_results_df

@st.cache_resource(show_spinner=False)
def query_fasst_api_peaks(precursor_mz, charge, peaks, database, host="https://api.fasst.gnps2.org", analog=False, precursor_mz_tol=0.05, fragment_mz_tol=0.05, min_cos=0.7, lower_delta=100, upper_delta=100, blocking=True):
    spectrum_query = {
        "peaks": peaks,
        "precursor_mz": precursor_mz,
        "charge": charge
    }

    params = {
        "library": database,
        "query_spectrum": json.dumps(spectrum_query),
        "analog": "Yes" if analog else "No",
        "cache": "No",
        "lower_delta": lower_delta,
        "upper_delta": upper_delta,
        "pm_tolerance": precursor_mz_tol,
        "fragment_tolerance": fragment_mz_tol,
        "cosine_threshold": min_cos
    }

    query_url = f"{host.rstrip('/')}/search"

    try:
        r = requests.post(query_url, json=params, timeout=30)
        try:
            r.raise_for_status()
        except requests.exceptions.RequestException as e:
            st.exception("POST to FASST returned error: %s", e)
            try:
                st.error("Response status: %s", r.status_code)
                st.error("Response text: %s", r.text[:2000])
            except Exception:
                pass
            raise

        task_id = r.json().get("id")

        if not task_id:
            # If the API didn't return an id, return the response for debugging
            st.error("FASST POST response did not include task id. Response: %s", r.text)
            raise Exception("FASST API response missing task id")

        params["task_id"] = task_id
        params["status"] = "PENDING"
        
        if blocking is False:
            return params

        return get_results(params, host=host)

    except KeyError:
        raise

@st.cache_data(show_spinner=False)
def execute_all_queries_sync(queries):
    output_results_list = []

    for query in tqdm(queries):
        try:
            results_dict = query_fasst_api_peaks(
                query["precursor_mz"],
                query["charge"],
                query["peaks"],
                query["database"],
                analog=query["analog"],
                precursor_mz_tol=query["precursor_mz_tol"],
                fragment_mz_tol=query["fragment_mz_tol"],
                min_cos=query["min_cos"],
                lower_delta=query.get("lower_delta", 100),
                upper_delta=query.get("upper_delta", 100),
                blocking=True
            )

            if "status" in results_dict:
                continue

            results_df = pd.DataFrame(results_dict.get("results", []))

            # adding the scan number information
            #results_df["query_scan"] = query["query_scan"]
            results_df["Scan Number"] = query["query_scan"]
            results_df["Library"] = query["database"]
            output_results_list.append(results_df)

        except KeyboardInterrupt:
            raise
        except Exception as e:
            print("Error in Query Execution:", e)
            pass

    if len(output_results_list) == 0:
        print("No results found")
        return pd.DataFrame()

    output_results_df = pd.concat(output_results_list, ignore_index=True)

    return output_results_df

@st.cache_resource
def get_results(query_parameters_dictionary, host="https://api.fasst.gnps2.org", blocking=True):
    task_id = query_parameters_dictionary.get("task_id")
    if not task_id:
        st.error("get_results called without task_id: %s", query_parameters_dictionary)
        raise ValueError("task_id is required to fetch results")
    
    retries_max = 120
    current_retries = 0
    while True:
        try:
            r = requests.get(f"{host.rstrip('/')}/search/result/{task_id}", timeout=30)
            r.raise_for_status()
        except requests.exceptions.RequestException as e:
            st.exception("Error fetching results status from FASST: %s", e)
            if blocking is False:
                return "PENDING"
            time.sleep(1)
            current_retries += 1
            if current_retries >= retries_max:
                raise Exception("Timeout waiting for results from FASST API")
            continue

        # checking if the results are ready
        try:
            j = r.json()
        except Exception as e:
            st.exception("Failed to parse JSON from FASST result: %s", e)
            if blocking is False:
                return "PENDING"
            time.sleep(1)
            current_retries += 1
            if current_retries >= retries_max:
                raise Exception("Timeout waiting for results from FASST API (no JSON)")
            continue

        if "status" in j and j["status"] == "PENDING":
            if blocking is False:
                return "PENDING"
            time.sleep(1)
            current_retries += 1

            if current_retries >= retries_max:
                raise Exception("Timeout waiting for results from FASST API")
            continue

        # results ready
        return j


if __name__ == "__main__":
    st.title("MS2 Scan Analyzer")
    with st.expander("ℹ️ About This App", expanded=False):
        st.markdown("""
        **MS2 Scan Analyzer** is an interactive web application for analyzing mass spectrometry MS2 data files. 
        
        This tool is designed for researchers working with metabolomics and mass spectrometry data, streamlining the process from raw data upload to advanced spectral matching and visualization.
        
        ---          
        The app provides the following features:

        - **📁 File Upload & Metadata Extraction:**  
        Upload `.mgf` or `.mzML` files and automatically extract scan metadata, including precursor m/z, charge state, and spectrum information.

        - **🔍 Scan Selection & Visualization:**  
        Select individual scans to view and compare raw and filtered MS2 spectra using interactive plots.

        - **🔗 GNPS FasstSearch Integration:**  
        Search your scans against GNPS spectral libraries using customizable parameters (library, analog search, mass tolerances, cosine similarity threshold, delta mass ranges).

        - **📊 Batch Search for All Scans:**  
        Run GNPS FasstSearch for all scans in your file, with real-time progress bar and estimated time remaining.

        - **🔬 Matching USI Results & Metabolomics Resolver:**  
        View matching USIs for selected scans, and visualize mirror plots comparing your scan to library spectra.

        - **📄 Downloadable Results:**  
        Download search results for individual scans or all scans as CSV files.

        - **💾 Session Persistence:**  
        Results are cached and reused across reruns for faster interaction.

        """)

    uploaded_file = st.file_uploader("Choose a file", type=["mgf", "mzML"])
    try:
        LIBRARIES = get_databases()
    except Exception:
        LIBRARIES = ["NORMAN", "ORNL_Bioscales2", "ORNL_Populus_LC_MSMS", "gnpsdata_index", "gnpsdata_test_index", "gnpsibrary", "massivedata_index", "massivekb_index", "metabolomicspanrepo_index_latest", "metabolomicspanrepo_index_nightly", "panrepo_2024_11_12", "panrepo_2025_07_09", "panrepo_2025_08_06", "ptfi2_index"]

    if uploaded_file is not None:
        file_type = uploaded_file.name.split(".")[-1].lower()
        file_name = uploaded_file.name.split(".")[0]
        
        with st.spinner("Reading File and Extracting Scan Metadata..."):
            if file_type == "mgf":
                scan_numbers, scan_metadata = read_mgf_file(uploaded_file)
            elif file_type == "mzml":
                scan_numbers, scan_metadata = read_mzml_file(uploaded_file)
            else:
                st.error("Unsupported file type. Please upload an MGF or mzML file.")
                scan_numbers, scan_metadata = [], []

        scan_numbers = [str(sn) for sn in scan_numbers]
        for scan in scan_metadata:
            scan["Scan Number"] = str(scan["Scan Number"])
        
        if not scan_metadata:
            st.error("No scans found in the uploaded file.")
        else:
            def reset_scan():
                st.session_state.clicked1 = False
                st.session_state.clicked2 = False
                st.session_state.clicked3 = False
                st.session_state.scan_csv_data = None
                st.session_state.all_csv_data = None

            scan_input = st.selectbox("Select Scan Number to Run Analysis", options=[""] + ["ALL SCANS"] + scan_numbers, on_change=reset_scan)

            #File Information
            st.header(f"{file_type.upper()} File Information ")

            with st.spinner("Loading Scan Metadata..."):
                df = pd.DataFrame(scan_metadata, columns=["Scan Number", "Spectrum ID", "Precursor M/Z", "Charge State", "SMILES ID"])

            #Scan Number Selection
            with st.expander("All Scan Metadata", expanded=False):
                st.write("Total Scans Found: ", len(scan_metadata))
                table = st.dataframe(df, hide_index=True)

             #Spectrum Visualization
            if (scan_input != "ALL SCANS") and (scan_input != ""):
                if scan_input:
                    user_scan = next((scan for scan in scan_metadata if scan["Scan Number"] == scan_input), None)
                    mz_filtered, sqrt_filtered, normal_filtered = peak_filtering(user_scan)
                    with st.expander(f"View Spectrums for Scan {user_scan['Scan Number']}", expanded=False):
                        show_plots = st.checkbox("Load Spectrum Plots")
                        spectrum_option = st.tabs(["Unfiltered Spectrum", "Filtered Spectrum - Normalized", "Filtered Spectrum - Square Root Normalized"])
                        if show_plots:    
                            with spectrum_option[0]:
                                peak_visual(user_scan["m/z data"], user_scan["intensity data"], str(user_scan["Scan Number"]), user_scan["Precursor M/Z"], user_scan["Charge State"])
                            with spectrum_option[1]:
                                peak_visual(mz_filtered, normal_filtered, str(user_scan["Scan Number"]), user_scan["Precursor M/Z"], user_scan["Charge State"])
                            with spectrum_option[2]:
                                peak_visual(mz_filtered, sqrt_filtered, str(user_scan["Scan Number"]), user_scan["Precursor M/Z"], user_scan["Charge State"])

            #GNPS FasstSearch Integration
            if scan_input:
                if scan_input != "ALL SCANS":
                    st.header(f"GNPS FasstSearch for Scan {scan_input}")
                else:
                    st.header("GNPS FasstSearch for All Scans")

                def reset_clicks():
                    st.session_state.clicked1 = False
                    st.session_state.clicked2 = False
                    st.session_state.clicked3 = False
                    st.session_state.scan_csv_data = None
                    st.session_state.all_csv_data = None    
                
                #st.info("After adjusting parameters, click the button below to run GNPS FASSTSearch. Depending on the number of scans in the uploaded file, this process may take several minutes.")

                library_select = st.selectbox("Select Library", LIBRARIES, key="library_select", index=3, on_change=reset_clicks)
                col1, col2 = st.columns(2)
                with col1:
                    analog_select = st.selectbox("Analog Search", ["No", "Yes"], key="analog_select", on_change=reset_clicks)
                    delta_mass_below = st.number_input("Delta Mass Below (Da)", value=130, min_value=0, max_value=300, key="delta_mass_below", on_change=reset_clicks)
                    delta_mass_above = st.number_input("Delta Mass Above (Da)", value=200, min_value=0, max_value=300, key="delta_mass_above", on_change=reset_clicks)
                    if delta_mass_below > delta_mass_above:
                        st.error("Delta Mass Below cannot be greater than Delta Mass Above.")
                with col2:
                    pm_tolerance = st.number_input("PM Tolerance (Da)", min_value=0.0, value=0.05, step=0.05, max_value=0.4, key="pm_tolerance", on_change=reset_clicks)
                    fragment_tolerance = st.number_input("Fragment Mass Tolerance (Da)", min_value=0.0, value=0.05, step=0.05, max_value=0.4, key="fragment_tolerance", on_change=reset_clicks)
                    cosine_threshold = st.number_input("Cosine Similarity Threshold", min_value=0.0, max_value=1.0, value=0.7, step=0.05, key="cosine_threshold", on_change=reset_clicks)

                analog_select = "Yes" if analog_select == "Yes" else "No"
                
                def click_button1():
                    st.session_state.clicked1 = True
                    st.session_state.results_df = None  # Reset so we know to regenerate

                if 'clicked1' not in st.session_state:
                    st.session_state.clicked1 = False 

                if scan_input == "ALL SCANS":           
                    st.button(f"View GNPS FASST Search Results for All Scans", on_click=click_button1)
                else: 
                    st.button(f"View GNPS FASST Search Results for Scan {scan_input}", on_click=click_button1)
                
                if st.session_state.clicked1:
                    try:
                        st.subheader("Matching USI Results")
                        # Check if results_df is already in session_state
                        if "results_df" in st.session_state and st.session_state.results_df is not None:
                            results_df = st.session_state.results_df
                        else:
                            if scan_input != "ALL SCANS":
                                with st.spinner("Running GNPS FASSTSearch...."):
                                    user_scan = next((scan for scan in scan_metadata if scan["Scan Number"] == scan_input), None)
                                    if user_scan is None:
                                        st.error("Selected scan not found.")
                                    else:
                                        precursor_mz = user_scan["Precursor M/Z"]
                                        charge = user_scan["Charge State"]
                                        peaks = user_scan["peaks"]
                                        api_response = query_fasst_api_peaks(precursor_mz, charge, peaks, library_select, analog=(analog_select == "Yes"), precursor_mz_tol=pm_tolerance, fragment_mz_tol=fragment_tolerance, min_cos=cosine_threshold, lower_delta=delta_mass_below, upper_delta=delta_mass_above, blocking=True)
                                        results_df = pd.DataFrame(api_response.get("results", []), columns=["Delta Mass", "USI", "Charge", "Cosine", "Matching Peaks", "Dataset", "Status"])
                            else: 
                                progress_text = "Processing all scans. Please wait..."
                                my_bar = st.progress(0, text=progress_text)
                                eta = st.empty()
                                scan_status = st.empty()
                                start_time = time.time()

                                # Build queries for all scans
                                queries = []
                                for scan in scan_metadata:
                                    queries.append({
                                        "precursor_mz": scan["Precursor M/Z"],
                                        "charge": scan["Charge State"],
                                        "peaks": scan["peaks"],
                                        "database": library_select,
                                        "analog": (analog_select == "Yes"),
                                        "precursor_mz_tol": pm_tolerance,
                                        "fragment_mz_tol": fragment_tolerance,
                                        "min_cos": cosine_threshold,
                                        "lower_delta": delta_mass_below,
                                        "upper_delta": delta_mass_above,
                                        "query_scan": scan["Scan Number"],
                                        "datasource": file_type
                                    })

                                total_scans = len(queries)
                                batch_size = 50
                                results_list = []
                                scans_processed = 0

                                if total_scans > batch_size: 
                                    for batch_idx, batch_start in enumerate(range(0, total_scans, batch_size)):
                                        batch_end = min(batch_start + batch_size, total_scans)
                                        batch = queries[batch_start:batch_end]
                                        #scan_status.text(f"Processing scans {batch_start + 1} to {batch_end} of {total_scans}")
                                        batch_start_time = time.time()
                                        batch_results = execute_all_queries_batch(batch)
                                        batch_duration = time.time() - batch_start_time

                                        scans_processed += (batch_end - batch_start)
                                        elapsed_time = time.time() - start_time
                                        avg_time_per_scan = elapsed_time / scans_processed if scans_processed else 0
                                        scans_left = total_scans - scans_processed
                                        eta_seconds = int(avg_time_per_scan * scans_left)
                                        eta.text(f"Estimated time remaining: {eta_seconds // 60} min {eta_seconds % 60} sec")

                                        my_bar.progress(scans_processed / total_scans, text=f"Processed {scans_processed} of {total_scans} scans")
                                        if not batch_results.empty:
                                            results_list.append(batch_results)
                                else: 
                                    for idx, query in enumerate(queries):
                                        single_result = execute_all_queries_sync([query])
                                        scans_processed += 1
                                        elapsed_time = time.time() - start_time
                                        avg_time_per_scan = elapsed_time / scans_processed if scans_processed else 0
                                        scans_left = total_scans - scans_processed
                                        eta_seconds = int(avg_time_per_scan * scans_left)
                                        eta.text(f"Estimated time remaining: {eta_seconds // 60} min {eta_seconds % 60} sec")
                                        my_bar.progress(scans_processed / total_scans, text=f"Processed {scans_processed} of {total_scans} scans")
                                        if not single_result.empty:
                                            results_list.append(single_result)

                                my_bar.empty()
                                eta.empty()
                                scan_status.empty()

                                if results_list is not None and len(results_list) > 0:
                                    results_df = pd.concat(results_list, ignore_index=True)
                                    desired_columns = [
                                        "Scan Number", "Library", "Delta Mass", "USI", "Charge", "Cosine", "Matching Peaks", "Dataset", "Status"
                                    ]
                                    missing_cols = set(desired_columns) - set(results_df.columns)
                                    for col in missing_cols:
                                        results_df[col] = np.nan
                                    results_df = results_df[desired_columns]
                                else:
                                    st.error("No matching results found for any scan in the selected library.")

                        st.session_state.results_df = results_df  # Save to session_state

                        # Use results_df from session_state for display and downloads
                        results_df = st.session_state.get("results_df", None)

                        if results_df is not None and not results_df.empty:
                            with st.expander(f"View Matching USI Results for {scan_input}", expanded=False):
                                st.write("Total Matching USIs Found: ", len(results_df))
                                results = st.dataframe(results_df, hide_index=True)
                            
                            st.subheader("Metabolomics Resolver Spectrum Viewer")
                            if scan_input != "ALL SCANS":
                                selected_usi = st.selectbox("Select Matching USI", options=[""] + list(results_df["USI"].unique()))
                            else:
                                col1, col2 = st.columns(2)
                                with col1:
                                    st.selectbox("Select a Scan Number", options=[""] + list(results_df["Scan Number"].unique()), key="select_scan_number")
                                with col2:
                                    selected_scan_number = st.session_state.get("select_scan_number")
                                    matching_usis = results_df[results_df["Scan Number"] == selected_scan_number]["USI"].unique().tolist()
                                    st.selectbox("Select Matching USI", options=[""] + matching_usis, key="selected_usi")

                            if scan_input == "ALL SCANS":
                                #scan_input = st.session_state.get("scan_input")
                                user_scan = next((scan for scan in scan_metadata if scan["Scan Number"] == selected_scan_number), None)
                                if user_scan is not None:
                                    mz_filtered, _, normal_filtered = peak_filtering(user_scan)
                                selected_usi = st.session_state.get("selected_usi")
                                
                            with st.spinner("Loading Metabolomics Resolver Spectrums..."): 
                                if selected_usi:
                                    spectrum_tabs = st.tabs(["Unfiltered Spectrum", "Filtered Spectrum"])
                                    with spectrum_tabs[0]:
                                        try:
                                            spectrum_top = sus.MsmsSpectrum(
                                                mz=user_scan["m/z data"],
                                                intensity=user_scan["intensity data"],
                                                identifier=str(user_scan["Scan Number"]),
                                                precursor_mz=user_scan["Precursor M/Z"],
                                                precursor_charge=user_scan["Charge State"]
                                            )
                                            spectrum_bottom = sus.MsmsSpectrum.from_usi(
                                                selected_usi,
                                                precursor_mz=user_scan["Precursor M/Z"],
                                                precursor_charge=user_scan["Charge State"]
                                            )
                                            
                                            create_mirror_plot(spectrum_top, spectrum_bottom)
                                        except Exception as e:
                                            st.error(f"Currently unable to create mirror plot, please try again later.")

                                    with spectrum_tabs[1]:
                                        try:
                                            spectrum_top = sus.MsmsSpectrum(
                                                mz=mz_filtered,
                                                intensity=normal_filtered,
                                                identifier=str(user_scan["Scan Number"]),
                                                precursor_mz=user_scan["Precursor M/Z"],
                                                precursor_charge=user_scan["Charge State"]
                                            )

                                            spectrum_bottom = sus.MsmsSpectrum.from_usi(
                                                selected_usi,
                                                precursor_mz=user_scan["Precursor M/Z"],
                                                precursor_charge=user_scan["Charge State"]
                                            )

                                            usi_scan = {
                                                "m/z data": spectrum_bottom.mz,
                                                "intensity data": spectrum_bottom.intensity,
                                                "Scan Number": spectrum_bottom.identifier,
                                                "Precursor M/Z": spectrum_bottom.precursor_mz,
                                                "Charge State": spectrum_bottom.precursor_charge
                                            }
                                            usi_mz_filtered, _, usi_normal_filtered = peak_filtering(usi_scan)
                                            spectrum_bottom = sus.MsmsSpectrum(
                                                mz=usi_mz_filtered,
                                                intensity=usi_normal_filtered,
                                                identifier=str(usi_scan["Scan Number"]),
                                                precursor_mz=usi_scan["Precursor M/Z"],
                                                precursor_charge=usi_scan["Charge State"]
                                            )
                                            create_mirror_plot(spectrum_top, spectrum_bottom)
                                        except Exception as e:
                                            st.error(f"Currently unable to create mirror plot, please try again later.")
                                    # else:
                                    #     st.info("Select a Scan Number and USI to view the Metabolomics Resolver Spectrum Viewer.")
                    
                            #Downloadable CSV Files
                            st.subheader("Download CSV Files")
                            st.write("Click the button below to download GNPS FASSTSearch Results with the above parameters for the selected scan(s) in a CSV file.")

                            def click_button2():
                                st.session_state.clicked2 = True
                                st.session_state.scan_csv_data = None  # Reset so we know to regenerate

                            def click_button3():
                                st.session_state.clicked3 = True
                                st.session_state.all_csv_data = None  # Reset so we know to regenerate

                            if 'clicked2' not in st.session_state:
                                st.session_state.clicked2 = False
                            if 'clicked3' not in st.session_state:
                                st.session_state.clicked3 = False

                            if scan_input != "ALL SCANS":
                                if results_df is not None and not results_df.empty: 
                                    selected_scan = next((scan for scan in scan_metadata if scan["Scan Number"] == scan_input), None)
                                    if selected_scan is not None:
                                        results_df["Library"] = library_select
                                        results_df["Scan Number"] = selected_scan["Scan Number"]
                                        desired_columns = ["Scan Number", "Library", "Delta Mass", "USI", "Charge", "Cosine", "Matching Peaks", "Dataset", "Status"]
                                        missing_cols = set(desired_columns) - set(results_df.columns)
                                        for col in missing_cols:
                                            results_df[col] = np.nan
                                        results_df = results_df[desired_columns]
                                        csv_data = results_df.to_csv(index=False).encode('utf-8')
                                        st.session_state.scan_csv_data = csv_data
                                        st.session_state.scan_csv_filename = f"{uploaded_file.name.split('.')[0]}_matchingUSIs_scan{scan_input}.csv"
                                    
                                    if st.session_state.get("scan_csv_data"):
                                        st.download_button(
                                            label=f"Download GNPS Results for Scan {scan_input} in CSV File",
                                            data=st.session_state.scan_csv_data,
                                            file_name=st.session_state.scan_csv_filename,
                                            mime="text/csv"
                                        )
                                else: 
                                    st.error("No results found, cannot generate CSV file.")
                            else:
                                if results_df is not None:
                                    csv_data = results_df.to_csv(index=False).encode('utf-8')
                                    st.session_state.all_csv_data = csv_data
                                    st.session_state.all_csv_filename = f"{uploaded_file.name.split('.')[0]}_matchingUSIs_allScans.csv"

                                    st.download_button(
                                        label="Download GNPS Results for All Scans in CSV File",
                                        data=st.session_state.all_csv_data,
                                        file_name=st.session_state.all_csv_filename,
                                        mime="text/csv"
                                    )
                                else: 
                                    st.error("No results found, cannot generate CSV file.")
                        else:
                            st.error("No matching results found. Please adjust your parameters and try again.")

                    except Exception as e:
                        st.error("No matching results found. Please adjust your parameters and try again.")

                