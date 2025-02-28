import streamlit as st  # allows for building web application
import matplotlib.pyplot as plt  # allows for plotting
import pandas as pd  # allows for data handling
import numpy as np  # allows for numerical operations
import spectrum_utils.plot as sup
import spectrum_utils.spectrum as sus
import plotly.tools as tls
import requests
import json
import urllib.parse
from urllib.parse import urlparse, parse_qs, unquote, urlencode

# def format_peaks(peaks):
#     #converts user_scan["peaks"] into a string where each peak is formatted as "m/z<TAB>intensity"
#     #only the top `max_peaks` (sorted by intensity) are included
#     sorted_peaks = sorted(peaks, key=lambda x: x[1], reverse=True)[:100]
#     return ((mz, intensity) for mz, intensity in sorted_peaks)

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

@st.cache_data
def build_url(peaks, precursor_mz, charge, library_select, analog_select, cache_select, delta_mass_below, delta_mass_above, pm_tolerance, fragment_tolerance, cosine_threshold):
    base_url = "https://fasst.gnps2.org/search"
    
    # Build the query spectrum in the expected JSON format.
    query_spectrum = {
        "n_peaks":len(peaks),
        "peaks":peaks,
        "precursor_charge":charge,
        "precursor_mz":precursor_mz
    }

    query_spectrum_json = json.dumps(query_spectrum)
    
    # Build the payload following the CLI code logic.
    data = {
        "usi": None,
        "library": "gnpslibrary",
        "analog": "Yes" if analog_select == "Yes" else "No",
        "cache": "Yes" if cache_select == "Yes" else "No", 
        "lower_delta": delta_mass_below,
        "upper_delta": delta_mass_above,
        "pm_tolerance": pm_tolerance,
        "fragment_tolerance": fragment_tolerance,
        "cosine_threshold": cosine_threshold,
        "query_spectrum": query_spectrum_json,
    }

    return {'url': base_url, 'data': data}

@st.cache_data
def print_POST(req):
    print('{}\n{}\r\n{}\r\n\r\n{}'.format(
        '-----------START-----------',
        req.method + ' ' + req.url,
        '\r\n'.join('{}: {}'.format(k, v) for k, v in req.headers.items()),
        req.body,
    ))

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

    # query_string = urllib.parse.urlencode(payload)
    # full_url = f"{base_url}?{query_string}"
    
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
    scans = []  # list to store parsed scans
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
    return scans, scan_numbers

if __name__ == "__main__":
    st.title("MS2 Scan")  # App title

    mgf_file = st.file_uploader("Choose a file", type="mgf")
    
    if mgf_file is not None:
        scans, scan_nums = read_mgf_file(mgf_file)
        if not scans:
            st.error("No scans found in the uploaded MGF file.")
        else:
            # Create a DataFrame showing scan metadata.
            df = pd.DataFrame(scans, columns=["Scan Number", "Spectrum ID", "PEPMASS Number", "Charge State", "SMILES ID"])
            
            # Dropdown menu for selecting a scan number.
            scan_input = st.selectbox("Select Scan Number to view MS2 Spectrum", options=scan_nums)
            
            # Display the DataFrame in an expander.
            with st.expander("Show Scan Numbers and Metadata"):
                st.dataframe(df)
            
            if scan_input:
                user_scan = next((scan for scan in scans if scan["Scan Number"] == scan_input), None)
                if user_scan is None:
                    st.error("Selected scan not found.")
                else:
                    mz_filtered, sqrt_filtered, normal_filtered = peak_filtering(user_scan)
                    spectrum_option = st.selectbox("Choose Spectrum Type", ["Unfiltered Spectrum", "Filtered Spectrum - Normalized", "Filtered Spectrum - Square Root Normalized"])
                    
                    if st.button("View Spectrum"):
                        if spectrum_option == "Unfiltered Spectrum":
                            peak_visual(user_scan["m/z data"], user_scan["intensity data"],
                                        str(user_scan["Scan Number"]), user_scan["PEPMASS Number"], user_scan["Charge State"])
                        elif spectrum_option == "Filtered Spectrum - Normalized":
                            peak_visual(mz_filtered, normal_filtered,
                                        str(user_scan["Scan Number"]), user_scan["PEPMASS Number"], user_scan["Charge State"])
                        else:
                            peak_visual(mz_filtered, sqrt_filtered,
                                        str(user_scan["Scan Number"]), user_scan["PEPMASS Number"], user_scan["Charge State"])
                
                with st.expander(f"Generate GNPS FASTSearch for Scan {user_scan['Scan Number']}"):
                    precursor_mz = user_scan["PEPMASS Number"]
                    charge = user_scan["Charge State"]
                    col1, col2, col3 = st.columns(3)
                    
                    with col1:
                        #peak_select = st.selectbox("Select Peaks", ["Normalized Peaks", "Unfiltered Peaks", "Square Root Peaks"], key="peak_select")
                        library_select = st.selectbox("Select Library", ["gnpsdata_index", "ORNL_Bioscales2", "ORNL_Populus_LC_MSMS", "gnpsdata_test_index", "gnpslibrary", "massivedata_index", "massivekb_index", "metabolomicspanrepo_index_latest", "metabolomicspanrepo_index_nightly", "panrepo_2024_11_12"], key="library_select")
                        analog_select = st.selectbox("Analog Search", ["No", "Yes"], key="analog_select")
                    with col2:
                        delta_mass_below = st.number_input("Delta Mass Below (Da)", min_value=0.0, key="delta_mass_below")
                        delta_mass_above = st.number_input("Delta Mass Above (Da)", min_value=0.0, key="delta_mass_above")
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
                        
                        # Call the FASTSearch API using the incorporated code.
                        api_response = generate_webAPI(peaks, precursor_mz, charge, library_select, analog_select, cache_select, delta_mass_below, delta_mass_above, pm_tolerance, fragment_tolerance, cosine_threshold)
                        url_response = build_url(peaks, precursor_mz, charge, library_select, analog_select, cache_select, delta_mass_below, delta_mass_above, pm_tolerance, fragment_tolerance, cosine_threshold)

                        # st.write(url_response)
                        url = create_link(peaks, precursor_mz, charge, library_select)

                        if api_response:
                            if "results" in api_response:
                                try: 
                                    st.subheader("Matching Results")
                                    results_df = pd.DataFrame(api_response["results"])
                                    desired_columns = ["Delta Mass", "USI", "Charge", "Cosine", "Matching Peaks", "Dataset", "Status"]
                                    results_df = results_df[desired_columns]   
                                    st.dataframe(results_df)
                                except KeyError:
                                    print("Some expected columns are missing. Available columns:", results_df.columns)
                                # st.write(api_response["results"])
                                output = []
                                
                                for match in api_response["results"]:
                                    try:
                                        response = requests.post(url_response['url'], data=url_response['data'])
                                        json_response = response.json()
                                        if response.status_code == 200:
                                            result = {
                                                'Delta Mass': match['Delta Mass'],
                                                'GNPSLibraryAccession': match.get('GNPSLibraryAccession', 'Not Available'),
                                                'USI': match['USI'],
                                                'Charge': match['Charge'],
                                                'Cosine': match['Cosine'],
                                                'Matching Peaks': match['Matching Peaks'],
                                                'Unit Delta Mass': match['Unit Delta Mass'],
                                                'Dataset': match['Dataset'],
                                                'Status': match['Status'],
                                                #'Adduct': match['Adduct'],
                                                #'CompoundName': match['CompoundName'],
                                                'Query Filename': match['Query Filename'],
                                                'Query Scan': match['Query Scan'],
                                                'Index UnitPM': match['Index UnitPM'],
                                                'Index IdxInUnitPM': match['Index IdxInUnitPM'],
                                                'Filtered Input Spectrum Path': match['Filtered Input Spectrum Path'],
                                            }
                                            output.append(result)
                                        else:
                                            continue
                                    except Exception as e:
                                        print(e)

                                json_data = json.dumps({'results': output})

                                if len(user_scan["peaks"]) > 100:
                                    st.write("URL is too long to access, view text url and response-request data here instead.")
                                    st.write(json_data)
                                    st.markdown(url)
                                else:
                                    st.markdown(url)

                                # if len(user_scan["peaks"]) > 100:
                                #     st.write("Number of peaks is too high, using 100 most intense peaks")
                                
                                # st.markdown(url)
                            else:
                                st.write("NO RESULTS")
                                
