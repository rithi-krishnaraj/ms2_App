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
        ax.set_title("Top: User Scan, Bottom: Selected USI", fontsize=16)
        st.pyplot(fig)
    except Exception as e:
        st.error(f"Failed to create mirror plot: {e}")

#For MASST Querying
@st.cache_resource(show_spinner=False)
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
        for status in status_results_list:
            if status.get("status") == "DONE":
                results_df = pd.DataFrame(status.get("results", []))

                # adding the scan number information
                if "query_scan" in status:
                    results_df["query_scan"] = status["query_scan"]
                if "usi" in status:
                    results_df["query_usi"] = status["usi"]

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

            scan_input = st.selectbox("Select Scan Number to Run Analysis", options=[""] + scan_numbers, on_change=reset_scan)

            #File Information
            st.header(f"{file_type.upper()} File Information ")

            with st.spinner("Loading Scan Metadata..."):
                df = pd.DataFrame(scan_metadata, columns=["Scan Number", "Spectrum ID", "Precursor M/Z", "Charge State", "SMILES ID"])

            #Scan Number Selection
            with st.expander("All Scan Metadata", expanded=False):
                st.write("Total Scans Found: ", len(scan_metadata))
                table = st.dataframe(df, hide_index=True)

             #Spectrum Visualization
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

            #GNPS FaastSearch Integration
            if scan_input:
                st.header(f"GNPS FastSearch for Scan {scan_input}")

                def reset_clicks():
                    st.session_state.clicked1 = False
                    st.session_state.clicked2 = False
                    st.session_state.clicked3 = False
                    st.session_state.scan_csv_data = None
                    st.session_state.all_csv_data = None    

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

                if 'clicked1' not in st.session_state:
                    st.session_state.clicked1 = False            

                st.button(f"View GNPS FAST Search Results for Scan {scan_input}", on_click=click_button1)

                if st.session_state.clicked1:
                    try:
                        st.subheader("Matching Results")
                        with st.spinner("Running GNPS FASTSearch... This may take a few minutes."):
                            user_scan = next((scan for scan in scan_metadata if scan["Scan Number"] == scan_input), None)
                            if user_scan is None:
                                st.error("Selected scan not found.")
                            else:
                                precursor_mz = user_scan["Precursor M/Z"]
                                charge = user_scan["Charge State"]
                                peaks = user_scan["peaks"]
                                api_response = query_fasst_api_peaks(precursor_mz, charge, peaks, library_select, analog=(analog_select == "Yes"), precursor_mz_tol=pm_tolerance, fragment_mz_tol=fragment_tolerance, min_cos=cosine_threshold, lower_delta=delta_mass_below, upper_delta=delta_mass_above, blocking=True)
                                results_df = pd.DataFrame(api_response.get("results", []), columns=["Delta Mass", "USI", "Charge", "Cosine", "Matching Peaks", "Dataset", "Status"])
                        with st.expander(f"View Matching USI Results for Scan {scan_input}", expanded=False):
                            st.write("Total Matching USIs Found: ", len(results_df))
                            selected_usi = st.selectbox("Select Matching USI for Metabolomics Resolver Spectrum Viewer", options=[""] + list(results_df["USI"].unique()))
                            results = st.dataframe(results_df, hide_index=True)

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
                                        st.error(f"Mirror plot creation failed, selected USI may be from the uploaded file. - Try again with a different USI.")

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
                                        st.error(f"Mirror plot creation failed. Try again with a different USI.")
                    except KeyError as e:
                        st.error("No matching results found. Please adjust your parameters and try again.")
                #Downloadable CSV Files
                st.subheader("Download CSV Files")
                st.write("Click the buttons below to run complete GNPS FASTSearch against all libraries with the above parameters for the selected scan or all scans. Download results in a CSV file. Note: Depending on the number of scans in the uploaded file, this process may take several minutes to hours to complete.")
                col1, col2 = st.columns(2)

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

                with col1:
                    st.button(f"Generate Results for Scan {scan_input}", on_click=click_button2)
                    # Only generate if button clicked and not already generated
                    if st.session_state.clicked2 and st.session_state.get("scan_csv_data") is None:
                        start_time = time.time()
                        progress_text = "Processing selected scan against all libraries. Please wait..."
                        my_bar = st.progress(0, text=progress_text)
                        time_left_placeholder = st.empty()

                        selected_scan = next((scan for scan in scan_metadata if scan["Scan Number"] == scan_input), None)
                        all_results = []
                        total_libs = len(LIBRARIES)

                        def update_progress(completed, total):
                            percent = int(completed / total * 100)
                            elapsed = time.time() - start_time
                            if completed > 0:
                                est_total = elapsed / completed * total
                                est_left = est_total - elapsed
                                time_left_placeholder.info(f"Estimated time left: {est_left:.1f} seconds")
                            my_bar.progress(percent, text=progress_text)

                        if selected_scan is not None:
                            if file_type == "mgf":
                                custom = "mgf"
                            else:
                                custom = "mzml"

                            queries = []
                            for library in LIBRARIES:
                                query = { "datasource": custom, "precursor_mz": selected_scan["Precursor M/Z"], "charge": selected_scan["Charge State"], "peaks": selected_scan["peaks"], "database": library, "analog": (analog_select == "Yes"), "precursor_mz_tol": pm_tolerance, "fragment_mz_tol": fragment_tolerance, "min_cos": cosine_threshold, "lower_delta": delta_mass_below, "upper_delta": delta_mass_above, "query_scan": selected_scan["Scan Number"] }
                                queries.append(query)

                            results_list = []
                            for idx, query in enumerate(queries):
                                try:
                                    result = query_fasst_api_peaks(query["precursor_mz"], query["charge"], query["peaks"], query["database"], analog=query["analog"], precursor_mz_tol=query["precursor_mz_tol"], fragment_mz_tol=query["fragment_mz_tol"], min_cos=query["min_cos"], lower_delta=query["lower_delta"], upper_delta=query["upper_delta"], blocking=True)
                                    df = pd.DataFrame(result.get("results", []))
                                    if not df.empty:
                                        df["Library"] = query["database"]
                                        df["Scan Number"] = query["query_scan"]
                                        results_list.append(df)
                                except Exception as e:
                                    continue
                                update_progress(idx + 1, total_libs)

                            my_bar.empty()
                            time_left_placeholder.empty()

                            if results_list:
                                results_df = pd.concat(results_list, ignore_index=True)
                                desired_columns = ["Scan Number", "Library", "Delta Mass", "USI", "Charge", "Cosine", "Matching Peaks", "Dataset", "Status"]
                                missing_cols = set(desired_columns) - set(results_df.columns)
                                for col in missing_cols:
                                    results_df[col] = np.nan
                                results_df = results_df[desired_columns]
                                csv_data = results_df.to_csv(index=False).encode('utf-8')
                                st.session_state.scan_csv_data = csv_data
                                st.session_state.scan_csv_filename = f"{uploaded_file.name.split('.')[0]}_matchingUSIs_scan{scan_input}.csv"
                                st.success(f"File created successfully in {time.time()-start_time:.2f} seconds!")
                            else:
                                st.session_state.scan_csv_data = None
                                st.error("No matching results found for the selected scan in any library with the current parameters.")
                        else:
                            st.session_state.scan_csv_data = None
                            st.error("Selected scan not found.")

                    # Show download button if CSV is available
                    if st.session_state.get("scan_csv_data"):
                        st.download_button(
                            label="Download GNPS Results in CSV File",
                            data=st.session_state.scan_csv_data,
                            file_name=st.session_state.scan_csv_filename,
                            mime="text/csv"
                        )

                with col2:
                    st.button("Generate Results for All Scans", on_click=click_button3)
                    if st.session_state.clicked3 and st.session_state.get("all_csv_data") is None:
                        start_time = time.time()
                        progress_text = "Processing all scans. Please wait..."
                        my_bar = st.progress(0, text=progress_text)
                        time_left_placeholder = st.empty()

                        datasource = "mgf" if file_type == "mgf" else "mzml"
                        queries = []
                        for scan in scan_metadata:
                            for library in LIBRARIES:
                                query = {
                                    "datasource": datasource,
                                    "precursor_mz": scan["Precursor M/Z"],
                                    "charge": scan["Charge State"],
                                    "peaks": scan["peaks"],
                                    "database": library,
                                    "analog": (analog_select == "Yes"),
                                    "precursor_mz_tol": pm_tolerance,
                                    "fragment_mz_tol": fragment_tolerance,
                                    "min_cos": cosine_threshold,
                                    "lower_delta": delta_mass_below,
                                    "upper_delta": delta_mass_above,
                                    "query_scan": scan["Scan Number"]
                                }
                                queries.append(query)

                        total_tasks = len(queries)
                        batch_size = 50
                        results_list = []
                        for batch_start in range(0, total_tasks, batch_size):
                            batch_queries = queries[batch_start:batch_start + batch_size]
                            batch_results_df = execute_all_queries_batch(batch_queries)
                            if not batch_results_df.empty:
                                results_list.append(batch_results_df)
                            completed = min(batch_start + batch_size, total_tasks)
                            percent = int(completed / total_tasks * 100)
                            elapsed = time.time() - start_time
                            if completed > 0:
                                est_total = elapsed / completed * total_tasks
                                est_left = est_total - elapsed
                                time_left_placeholder.info(f"Estimated time left: {est_left:.1f} seconds")
                            my_bar.progress(percent, text=progress_text)

                        my_bar.empty()
                        time_left_placeholder.empty()

                        if results_list:
                            results_df = pd.concat(results_list, ignore_index=True)
                            if "query_scan" in results_df.columns:
                                results_df.rename(columns={"query_scan": "Scan Number"}, inplace=True)
                            if "query_usi" in results_df.columns:
                                results_df.rename(columns={"query_usi": "USI"}, inplace=True)
                            desired_columns = ["Scan Number", "Library", "Delta Mass", "USI", "Charge", "Cosine", "Matching Peaks", "Dataset", "Status"]
                            missing_cols = set(desired_columns) - set(results_df.columns)
                            for col in missing_cols:
                                results_df[col] = np.nan
                            results_df = results_df[desired_columns]
                            csv_data = results_df.to_csv(index=False).encode('utf-8')
                            st.session_state.all_csv_data = csv_data
                            st.session_state.all_csv_filename = f"{uploaded_file.name.split('.')[0]}_matchingUSIs_allScans.csv"
                            st.success(f"File created successfully in {time.time()-start_time:.2f} seconds!")
                        else:
                            st.session_state.all_csv_data = None
                            st.error("No matching results found for any scan in any library.")

                    if st.session_state.get("all_csv_data"):
                        st.download_button(
                            label="Download GNPS Results in CSV File",
                            data=st.session_state.all_csv_data,
                            file_name=st.session_state.all_csv_filename,
                            mime="text/csv"
                        )
