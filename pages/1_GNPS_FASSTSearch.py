import streamlit as st
import src.spectrum_visuals as sv
import pandas as pd
import numpy as np
import src.gnps_fasstsearch as gfs

try:
    LIBRARIES = gfs.get_databases()
except Exception:
    LIBRARIES = ["NORMAN", "ORNL_Bioscales2", "ORNL_Populus_LC_MSMS", "gnpsdata_index", "gnpsdata_test_index", "gnpsibrary", "massivedata_index", "massivekb_index", "metabolomicspanrepo_index_latest", "metabolomicspanrepo_index_nightly", "panrepo_2024_11_12", "panrepo_2025_07_09", "panrepo_2025_08_06", "ptfi2_index"]

st.session_state["LIBRARIES"] = LIBRARIES
LIBRARIES = st.session_state["LIBRARIES"]


if "file" in st.session_state and st.session_state["file"] is not None:
    uploaded_file = st.session_state["file"]
    file_type = uploaded_file.name.split(".")[-1].lower()
    file_name = uploaded_file.name.split(".")[0]

    st.header(f"{file_type.upper()} File Information ")
    
    with st.spinner("Reading File and Extracting Scan Metadata..."):
        if file_type == "mgf":
            scan_numbers, scan_metadata = sv.read_mgf_file(uploaded_file)
        elif file_type == "mzml":
            scan_numbers, scan_metadata = sv.read_mzml_file(uploaded_file)
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

        with st.spinner("Loading Scan Metadata..."):
            df = pd.DataFrame(scan_metadata, columns=["Scan Number", "Spectrum ID", "Precursor M/Z", "Charge State", "SMILES ID"])

        #Scan Number Selection
        with st.expander("All Scan Metadata", expanded=False):
            st.write("Total Scans Found: ", len(scan_metadata))
            table = st.dataframe(df, hide_index=True)

        scan_input = st.selectbox("Select Scan Number to Run Analysis", options=[""] + ["ALL SCANS"] + scan_numbers, on_change=reset_scan)
        st.session_state["scan input"] = scan_input
        scan_input = st.session_state["scan input"]

        #Spectrum Visualization
        if (scan_input != "ALL SCANS") and (scan_input != ""):
            if scan_input:
                user_scan = next((scan for scan in scan_metadata if scan["Scan Number"] == scan_input), None)
                mz_filtered, sqrt_filtered, normal_filtered = sv.peak_filtering(user_scan)
                with st.expander(f"View Spectrums for Scan {user_scan['Scan Number']}", expanded=False):
                    show_plots = st.checkbox("Load Spectrum Plots")
                    spectrum_option = st.tabs(["Unfiltered Spectrum", "Filtered Spectrum - Normalized", "Filtered Spectrum - Square Root Normalized"])
                    if show_plots:    
                        with spectrum_option[0]:
                            sv.peak_visual(user_scan["m/z data"], user_scan["intensity data"], str(user_scan["Scan Number"]), user_scan["Precursor M/Z"], user_scan["Charge State"])
                        with spectrum_option[1]:
                            sv.peak_visual(mz_filtered, normal_filtered, str(user_scan["Scan Number"]), user_scan["Precursor M/Z"], user_scan["Charge State"])
                        with spectrum_option[2]:
                            sv.peak_visual(mz_filtered, sqrt_filtered, str(user_scan["Scan Number"]), user_scan["Precursor M/Z"], user_scan["Charge State"])

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
                                    api_response = gfs.query_fasst_api_peaks(precursor_mz, charge, peaks, library_select, analog=(analog_select == "Yes"), precursor_mz_tol=pm_tolerance, fragment_mz_tol=fragment_tolerance, min_cos=cosine_threshold, lower_delta=delta_mass_below, upper_delta=delta_mass_above, blocking=True)
                                    results_df = pd.DataFrame(api_response.get("results", []), columns=["Delta Mass", "USI", "Charge", "Cosine", "Matching Peaks", "Dataset", "Status"])
                        else: 
                            progress_text = "Processing all scans. Please wait..."
                            my_bar = st.progress(0, text=progress_text)
                            eta = st.empty()
                            scan_status = st.empty()
                            start_time = gfs.time.time()

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
                                    batch_start_time = gfs.time.time()
                                    batch_results = gfs.execute_all_queries_batch(batch)
                                    batch_duration = gfs.time.time() - batch_start_time

                                    scans_processed += (batch_end - batch_start)
                                    elapsed_time = gfs.time.time() - start_time
                                    avg_time_per_scan = elapsed_time / scans_processed if scans_processed else 0
                                    scans_left = total_scans - scans_processed
                                    eta_seconds = int(avg_time_per_scan * scans_left)
                                    eta.text(f"Estimated time remaining: {eta_seconds // 60} min {eta_seconds % 60} sec")

                                    my_bar.progress(scans_processed / total_scans, text=f"Processed {scans_processed} of {total_scans} scans")
                                    if not batch_results.empty:
                                        results_list.append(batch_results)
                            else: 
                                for idx, query in enumerate(queries):
                                    single_result = gfs.execute_all_queries_sync([query])
                                    scans_processed += 1
                                    elapsed_time = gfs.time.time() - start_time
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
                                    results_df[col] = gfs.np.nan
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
                                mz_filtered, _, normal_filtered = sv.peak_filtering(user_scan)
                            selected_usi = st.session_state.get("selected_usi")
                            
                        with st.spinner("Loading Metabolomics Resolver Spectrums..."): 
                            if selected_usi:
                                spectrum_tabs = st.tabs(["Unfiltered Spectrum", "Filtered Spectrum"])
                                with spectrum_tabs[0]:
                                    try:
                                        spectrum_top = sv.sus.MsmsSpectrum(
                                            mz=user_scan["m/z data"],
                                            intensity=user_scan["intensity data"],
                                            identifier=str(user_scan["Scan Number"]),
                                            precursor_mz=user_scan["Precursor M/Z"],
                                            precursor_charge=user_scan["Charge State"]
                                        )
                                        spectrum_bottom = sv.sus.MsmsSpectrum.from_usi(
                                            selected_usi,
                                            precursor_mz=user_scan["Precursor M/Z"],
                                            precursor_charge=user_scan["Charge State"]
                                        )
                                        
                                        gfs.create_mirror_plot(spectrum_top, spectrum_bottom)
                                    except Exception as e:
                                        st.error(f"Currently unable to create mirror plot, please try again later.")

                                with spectrum_tabs[1]:
                                    try:
                                        spectrum_top = sv.sus.MsmsSpectrum(
                                            mz=mz_filtered,
                                            intensity=normal_filtered,
                                            identifier=str(user_scan["Scan Number"]),
                                            precursor_mz=user_scan["Precursor M/Z"],
                                            precursor_charge=user_scan["Charge State"]
                                        )

                                        spectrum_bottom = sv.sus.MsmsSpectrum.from_usi(
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
                                        usi_mz_filtered, _, usi_normal_filtered = sv.peak_filtering(usi_scan)
                                        spectrum_bottom = sv.sus.MsmsSpectrum(
                                            mz=usi_mz_filtered,
                                            intensity=usi_normal_filtered,
                                            identifier=str(usi_scan["Scan Number"]),
                                            precursor_mz=usi_scan["Precursor M/Z"],
                                            precursor_charge=usi_scan["Charge State"]
                                        )
                                        gfs.create_mirror_plot(spectrum_top, spectrum_bottom)
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
else:
    st.warning("Please upload a file in the 'File Upload' page to use GNPS FASSTSearch.")