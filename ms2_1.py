import streamlit as st  # allows for building web application
import matplotlib.pyplot as plt  # allows for plotting
import pandas as pd  # allows for data handling
import numpy as np  # allows for numerical operations
import spectrum_utils.plot as sup
import spectrum_utils.spectrum as sus
import plotly.tools as tls
import requests
import json
from urllib.parse import urlparse, parse_qs, unquote

def format_peaks(peaks, max_peaks=100):
    #converts user_scan["peaks"] into a string where each peak is formatted as "m/z<TAB>intensity"
    #only the top `max_peaks` (sorted by intensity) are included
    sorted_peaks = sorted(peaks, key=lambda x: x[1], reverse=True)[:max_peaks]
    return "\n".join(f"{mz}\t{intensity}" for mz, intensity in sorted_peaks)

def generate_webAPI(peaks, precursor_mz, charge, library_select, analog_select, delta_mass_below, delta_mass_above, pm_tolerance, fragment_tolerance, cosine_threshold):

    base_url = "https://fasst.gnps2.org/fastsearch/#"
    
    # gets top 100 peaks from user scan sorted by intensity
    peaks_str = format_peaks(peaks, max_peaks=100)
    
    #query parameters 
    params = {
        "usi1": "None",
        "precursor_mz": precursor_mz,
        "charge": charge,
        "library_select": library_select,
        "analog_select": analog_select,
        "delta_mass_below": delta_mass_below,
        "delta_mass_above": delta_mass_above,
        "pm_tolerance": pm_tolerance,
        "fragment_tolerance": fragment_tolerance,
        "cosine_threshold": cosine_threshold,
        "use_peaks": 1,
        "peaks": peaks_str
    }
    
    try:
        response = requests.get(base_url, params=params, timeout=30)
        response.raise_for_status()
        return response.url
    except requests.exceptions.RequestException as e:
        st.error(f"API request failed: {e}")
        return None

def decode_gnps_response(url):
    parsed = urlparse(url)
    
    # Get query parameters; parse_qs returns lists so we flatten them.
    query_params = {k: v[0] for k, v in parse_qs(parsed.query).items()}
    
    # Get the fragment, which should be a percent-encoded JSON string.
    fragment_data = {}
    if parsed.fragment:
        decoded_fragment = unquote(parsed.fragment)
        try:
            fragment_data = json.loads(decoded_fragment)
        except json.JSONDecodeError as e:
            print("Error decoding JSON fragment:", e)
    
    # Optionally parse the peaks string into a list of (m/z, intensity) pairs.
    peaks_list = []
    if "peaks" in fragment_data:
        peaks_str = fragment_data["peaks"]
        # Split into lines and then each line into components (assuming tab-delimited)
        for line in peaks_str.strip().split("\n"):
            if line:
                parts = line.split("\t")
                if len(parts) >= 2:
                    try:
                        mz = float(parts[0])
                        intensity = float(parts[1])
                        peaks_list.append((mz, intensity))
                    except ValueError:
                        # If conversion fails, skip or handle as needed.
                        continue
        # Store the parsed peaks list into our dictionary.
        fragment_data["peaks_list"] = peaks_list

    # Combine both dictionaries (query params and fragment data)
    decoded_data = {**query_params, **fragment_data}
    return decoded_data

@st.cache_data
def peak_filtering(user_scan): #DONE
    mz_array = user_scan["m/z data"]
    intensity_array = user_scan["intensity data"]
    filtered_mz = []
    filtered_intensities = []
    pepmass_val = float(user_scan["PEPMASS Number"])

    #basic peak filtering
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
def peak_normalizing(filtered_intensities): #DONE

    normalized_intensities = np.copy(filtered_intensities) / np.linalg.norm(filtered_intensities)
    sqrt_intensities = np.sqrt(normalized_intensities)

    return sqrt_intensities, normalized_intensities

def peak_visual(mzs, intensities, scanNum, pepmass, charge): #DONE
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
def read_mgf_file(mgf_file): #DONE
    scans = [] #list to store parsed scans
    current_scan = None
    scan_numbers = []

    file = mgf_file.read().decode('utf-8').splitlines()
    for line in file: #for each line in the file
        line = line.strip()
        if line == "BEGIN IONS": #beginning of scan 
            current_scan = {"Scan Number": 0, "Spectrum ID": '', "PEPMASS Number": 0.0, "Charge State": 0, "SMILES ID": '', "peaks": [], "m/z data": [], "intensity data": []} #initializes new scan with keys
        elif line == "END IONS": #end of scan
            if current_scan:
                scans.append(current_scan) #adding current scan to total scans
                current_scan = None #Reseting for next scan in mgf file
        elif current_scan is not None: #if scan has begun
            if "=" in line:
                data = line.split('=', 1) #limits line split to first '='
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
            else: #must be peak data
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

    mgf_file = st.file_uploader("Choose a file", type="mgf")  # allows user to upload file
    
    if mgf_file is not None:
        scans, scan_nums = read_mgf_file(mgf_file)
        if not scans:
            st.error("No scans found in the uploaded MGF file.")
        else:
            # Create a DataFrame showing scan metadata
            df = pd.DataFrame(scans, columns=["Scan Number", "Spectrum ID", "PEPMASS Number", "Charge State", "SMILES ID"])
            
            # Dropdown menu for selecting a scan number
            scan_input = st.selectbox("Select Scan Number to view MS2 Spectrum", options=scan_nums)
            
            # Display the DataFrame in an expander
            with st.expander("Show Scan Numbers and Metadata"):
                st.dataframe(df)
            
            if scan_input:
                # Find the scan corresponding to the selected scan number
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
                        peak_select = st.selectbox("Select Peaks", ["Normalized Peaks", "Unfiltered Peaks", "Square Root Peaks"],key="peak_select")
                        library_select = st.selectbox("Select Library", ["gnpsdata_index", "other_library"], key="library_select")
                        analog_select = st.selectbox("Analog Search", ["No", "Yes"], key="analog_select")
                    with col2:
                        delta_mass_below = st.number_input("Delta Mass Below (Da)", min_value=0.0, key="delta_mass_below")
                        delta_mass_above = st.number_input("Delta Mass Above (Da)", min_value=0.0, key="delta_mass_above")
                    with col3:
                        pm_tolerance = st.number_input("PM Tolerance (Da)", min_value=0.0, value=0.05, key="pm_tolerance")
                        fragment_tolerance = st.number_input("Fragment Mass Tolerance (Da)", min_value=0.0, value=0.05, key="fragment_tolerance")
                        cosine_threshold = st.number_input("Cosine Similarity Threshold", min_value=0.0, max_value=1.0, value=0.7, key="cosine_threshold")
                    
                    if st.button("Generate GNPS FASTSearch API"):
                        # Choose peaks based on user selection
                        if peak_select == "Normalized Peaks":
                            peaks = [(mz_filtered[i], normal_filtered[i]) for i in range(len(mz_filtered))]
                        elif peak_select == "Square Root Peaks":
                            peaks = [(mz_filtered[i], sqrt_filtered[i]) for i in range(len(mz_filtered))]
                        else:
                            peaks = user_scan["peaks"]
                        
                        # Call the FASTSearch API
                        api_response = generate_webAPI(peaks, precursor_mz, charge, library_select, analog_select, delta_mass_below, delta_mass_above, pm_tolerance, fragment_tolerance, cosine_threshold)
                        st.write("Request URL:", api_response)
                        decoded_results = decode_gnps_response(api_response)
                        st.write("Decoded GNPS FASTSearch Results:")
                        st.write(decoded_results)