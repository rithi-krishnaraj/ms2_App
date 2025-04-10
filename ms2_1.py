#output file that organizes matching results for each scan
import streamlit as st  # allows for building web application
import matplotlib.pyplot as plt  # allows for plotting
# from matplotlib.ticker import ScalarFormatter
import pandas as pd  # allows for data handling
import numpy as np  # allows for numerical operations
import spectrum_utils.plot as sup
import spectrum_utils.spectrum as sus
import plotly.tools as tls
import requests
import json
import io
import urllib.parse
from urllib.parse import urlparse, parse_qs, unquote, urlencode
from xlsxwriter import Workbook
LIBRARIES = ["gnpsdata_index", "ORNL_Bioscales2", "ORNL_Populus_LC_MSMS", "gnpsdata_test_index", "gnpslibrary", "massivedata_index", "massivekb_index", "metabolomicspanrepo_index_latest", "metabolomicspanrepo_index_nightly", "panrepo_2024_11_12"]

@st.cache_data
def fetch_spectrum_data(usi):
    url = f"https://metabolomics-usi.gnps2.org/json/?usi1="

    usi = usi.replace(":", "%3A")
    url += usi

    try:
        response = requests.get(url)
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as e:
        st.error(f"Failed to fetch spectrum data: {e}")
        return None

# @st.cache_data
# def create_mirror_plot(scan_mz, scan_intensity, scan_precursor_mz, scan_charge, scan_id,
#                        usi_mz, usi_intensity, usi_precursor_mz, usi_charge, usi_id):
#     # Create spectrum objects
#     scan_spectrum = sus.MsmsSpectrum(
#         identifier=f"Scan {scan_id}",
#         precursor_mz=scan_precursor_mz,
#         precursor_charge=scan_charge,
#         mz=scan_mz,
#         intensity=scan_intensity
#     )
#     usi_spectrum = sus.MsmsSpectrum(
#         identifier=usi_id,
#         precursor_mz=usi_precursor_mz,
#         precursor_charge=usi_charge,
#         mz=usi_mz,
#         intensity = usi_intensity
#     )

#     # Create figure and axis
#     fig, ax = plt.subplots(figsize=(10, 6))

#     # Generate the mirror plot using spectrum_utils.plot.mirror()
#     sup.mirror(spec_top=scan_spectrum, spec_bottom=usi_spectrum, ax=ax)

#     ax.xaxis.set_major_formatter(ScalarFormatter())
#     ax.yaxis.set_major_formatter(ScalarFormatter())
#     # Convert to Plotly for Streamlit
#     try:
#         plotly_fig = tls.mpl_to_plotly(fig)
#         plotly_fig.update_traces(hoverinfo="x+y")
#     except Exception as e:
#         st.error(f"Failed to create mirror plot: {e}")
#         plotly.fig = None

#     plt.close(fig)

#     return plotly_fig
    
@st.cache_data
def create_link(peaks, precursor_mz, charge, library_select): #creates json decoded
    base_url = "https://fasst.gnps2.org/search"
    query_spectrum = {
        "n_peaks":len(peaks),
        "peaks":peaks,
        "precursor_charge":charge,
        "precursor_mz":precursor_mz
    }

    query_spectrum_json = json.dumps(query_spectrum)
    
    payload = {
        "usi": None, 
        "library": library_select,
        "query_spectrum": query_spectrum_json,
    }
    
    # Manually construct the query string.
    query_string = "&".join(f"{key}={urllib.parse.quote_plus(str(value), safe='[]{}:,')}" for key, value in payload.items())
    query_string = query_string.replace("+", "")
    full_url = f"{base_url}?{query_string}"
    
    return full_url

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

@st.cache_data
def generate_webAPI(peaks, precursor_mz, charge, library_select, analog_select, cache_select, delta_mass_below, delta_mass_above, pm_tolerance, fragment_tolerance, cosine_threshold):

    base_url = "https://fasst.gnps2.org/search"
    
    # Build the query spectrum in the expected JSON format.
    query_spectrum = {
        "n_peaks": len(peaks),
        "peaks": peaks,
        "precursor_charge": charge,
        "precursor_mz": precursor_mz
    }

    query_spectrum_json = json.dumps(query_spectrum)
    
    # Build the payload following the CLI code logic.
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
        return response.json() # Return the API response as a dictionary
    except requests.exceptions.RequestException as e:
        st.error(f"API request failed: {e}")
        return None

@st.cache_data
def peak_filtering(user_scan):
    mz_array = user_scan["m/z data"]
    intensity_array = user_scan["intensity data"]
    filtered_mz = []
    filtered_intensities = []
    pepmass_val = float(user_scan["PEPMASS Number"])

    # Basic peak filtering.
    for i, mz in enumerate(mz_array):
        peak_range = [j for j in range(len(mz_array)) if abs(mz_array[j] - mz) <= 25]
        sorted_range = sorted(peak_range, key=lambda j: intensity_array[j], reverse=True)
        if i in sorted_range[:6]:
            if abs(mz - pepmass_val) > 17:
                filtered_mz.append(mz)
                filtered_intensities.append(intensity_array[i])

    sqrt_data, normalized_data = peak_normalizing(filtered_intensities)
    return filtered_mz, sqrt_data, normalized_data

@st.cache_data  
def peak_normalizing(filtered_intensities):
    normalized_intensities = np.copy(filtered_intensities) / np.linalg.norm(filtered_intensities)
    sqrt_intensities = np.sqrt(normalized_intensities)
    return sqrt_intensities, normalized_intensities

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

@st.cache_data
def read_mgf_file(mgf_file):
    scans = [] 
    current_scan = None
    scan_numbers = []
    all_results = []

    file = mgf_file.read().decode('utf-8').splitlines()
    for line in file:
        line = line.strip()
        if line == "BEGIN IONS":
            current_scan = {"Scan Number": 0, "Spectrum ID": '', "PEPMASS Number": 0.0, "Charge State": 0, "SMILES ID": '', "peaks": [], "m/z data": [], "intensity data": []}
        elif line == "END IONS":
            if current_scan:
                scans.append(current_scan)

                for library in LIBRARIES:
                    api_response = generate_webAPI(
                        peaks = current_scan["peaks"],
                        precursor_mz=current_scan["PEPMASS Number"],
                        charge=current_scan["Charge State"],
                        library_select=library,
                        analog_select="No",  # Default value
                        cache_select="Yes",  # Default value
                        delta_mass_below=130,
                        delta_mass_above=200,
                        pm_tolerance=0.05,
                        fragment_tolerance=0.05,
                        cosine_threshold=0.7
                    )
                    if api_response and "results" in api_response:
                        for result in api_response["results"]:
                            # Add scan number and library to each result
                            result["Scan Number"] = current_scan["Scan Number"]
                            result["Library"] = library
                            all_results.append(result)
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


    
    return scans, scan_numbers, all_results

if __name__ == "__main__":
    st.title("MS2 Scan Analyzer")  # App title

    mgf_file = st.file_uploader("Choose a file", type="mgf")
    
    if mgf_file is not None:
        scans, scan_nums, matching_results = read_mgf_file(mgf_file)
        if not scans:
            st.error("No scans found in the uploaded MGF file.")
        else:

        #Download button for all matching libraries
            if matching_results:
                results_df = pd.DataFrame(matching_results)
                desired_columns = ["Scan Number", "Library", "Delta Mass", "USI", "Charge", "Cosine", "Matching Peaks", "Dataset", "Status"]
                results_df = results_df[desired_columns]

                # Save the results to an Excel file
                output = io.BytesIO()
                with pd.ExcelWriter(output, engine="xlsxwriter") as writer:
                    results_df.to_excel(writer, index=False, sheet_name="GNPS Results")

                    # Group by scan and library for better organization
                    grouped = results_df.groupby("Scan Number")
                    for scan_number, group in grouped:
                        sheet_name = f"Scan{scan_number}"
                        group.to_excel(writer, index=False, sheet_name=sheet_name)

                #Excel file as a downloadable file
                st.download_button(
                    label="Download GNPS Results",
                    data=output.getvalue(),
                    file_name="gnps_results.xlsx",
                    mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
                )    

            # Create a DataFrame showing scan metadata.
            df = pd.DataFrame(scans, columns=["Scan Number", "Spectrum ID", "PEPMASS Number", "Charge State", "SMILES ID"])
            
            # Display the DataFrame in an expander.
            with st.expander("Show Scan Numbers and Metadata"):
                st.dataframe(df)

            # Dropdown menu for selecting a scan number.
            scan_input = st.selectbox("Select Scan Number to view MS2 Spectrum", options=scan_nums)
                        
            if scan_input:
                user_scan = next((scan for scan in scans if scan["Scan Number"] == scan_input), None)
                if user_scan is None:
                    st.error("Selected scan not found.")
                else:
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
                
                with st.expander(f"GNPS FASTSearch for Scan {user_scan['Scan Number']}", expanded = False):
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
                    
                    if st.button("Generate GNPS FASTSearch API"):
                        # Choose peaks based on user selection.
                        # if peak_select == "Normalized Peaks":
                        #     peaks = [(mz_filtered[i], normal_filtered[i]) for i in range(len(mz_filtered))]
                        # elif peak_select == "Square Root Peaks":
                        #     peaks = [(mz_filtered[i], sqrt_filtered[i]) for i in range(len(mz_filtered))]
                        # else:
                        # if len(user_scan["peaks"]) > 100:
                        #     peaks = list(format_peaks(user_scan["peaks"]))
                        # else:
                        #     peaks = user_scan["peaks"]

                        peaks = user_scan["peaks"]
                        mz = user_scan["m/z data"]
                        intensity = user_scan["intensity data"]
                        # Call the FASTSearch API using the incorporated code.
                        api_response = generate_webAPI(peaks, precursor_mz, charge, library_select, analog_select, cache_select, delta_mass_below, delta_mass_above, pm_tolerance, fragment_tolerance, cosine_threshold)
                        #Takes user to GNPS FASTSEARCH
                        gnpsLink = gnps_input_link(mz, intensity, precursor_mz, charge, library_select, analog_select, cache_select, delta_mass_below, delta_mass_above, pm_tolerance, fragment_tolerance, cosine_threshold)
                        
        
                        
                        # Display results if they exist (either from current run or previous run)
                        if api_response and "results" in api_response:
                            st.write("GNPS Input Link Created with Scan Data!")
                            st.markdown(f"[Go to GNPS Link]({gnpsLink})", unsafe_allow_html=True)
                            try: 
                                st.subheader("Matching Results")
                                results_df = pd.DataFrame(api_response["results"])
                                desired_columns = ["Delta Mass", "USI", "Charge", "Cosine", "Matching Peaks", "Dataset", "Status"]
                                results_df = results_df[desired_columns]
                                st.dataframe(results_df)
                                
                                # #USI selection for Spectrum Exploration
                                # st.subheader("Metabolomics Resolver Spectrum Viewer")
                                # usi_list = results_df["USI"].tolist()
                                # selected_usi = st.selectbox("Select a USI for Spectrum Exploration", usi_list,  key="usi_selector")
                                
                                # if selected_usi:
                                #     spectrum_data = fetch_spectrum_data(selected_usi)


                                #     if spectrum_data and "peaks" in spectrum_data:
                                #         spectrum_tabs = st.tabs(["Unfiltered Spectrum", "Filtered Spectrum"])
                                #         usi_mz = []
                                #         usi_intensity = []

                                #         for peak in spectrum_data["peaks"]:
                                #             usi_mz.append(peak[0])
                                #             usi_intensity.append(peak[1])

                                #         # usi_mz_filtered, usi_normal_filtered, usi_sqrt_filtered = peak_filtering()
                                #         with spectrum_tabs[0]:
                                #             plotly_fig = create_mirror_plot(
                                #             scan_mz = user_scan["m/z data"],
                                #             scan_intensity = user_scan["intensity data"],
                                #             scan_charge = user_scan["Charge State"],
                                #             scan_precursor_mz = user_scan["PEPMASS Number"],
                                #             scan_id = user_scan["Scan Number"],
                                #             usi_mz = usi_mz,
                                #             usi_intensity = usi_intensity,
                                #             usi_charge = spectrum_data["precursor_charge"],
                                #             usi_precursor_mz = spectrum_data["precursor_mz"],
                                #             usi_id = str(selected_usi)
                                #             )
                                #             st.plotly_chart(plotly_fig)
                                #         # with spectrum_tabs[1]:
                                #         #     scan_mz = mz_filtered
                                #         #     scan_intensity = normal_filtered
                                #         #     scan_charge = user_scan["Charge State"] 
                                #         #     scan_pepmass = user_scan["PEPMASS Number"]
                                #         #     scan_id = user_scan["Scan Number"]
                                #         #     usi_mz = usi_mz_filtered
                                #         #     usi_intensity = usi_normal_filtered
                                #         #     usi_charge = spectrum_data['precursor_charge']"]
                                #         #     usi_pepmass = spectrum_data["precursor_mz"]
                                #         #     usi_id = selected_usi
                                #         #     create_mirror_plot(scan_mz, scan_intensity, scan_charge, scan_pepmass, usi_mz, usi_intensity, usi_charge, usi_pepmass)
                                #     else:
                                #         st.error("Failed to fetch spectrum data")
                            except KeyError as e:
                                st.error(f"Could not retrieve spectrum data for the selected USI: {e}")     
                        else:
                            st.error("NO RESULT")
                                                       