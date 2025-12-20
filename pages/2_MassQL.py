import streamlit as st
import pandas as pd
from massql import msql_parser
from massql import msql_engine
import os
from pathlib import Path
import requests
import io

# Page configuration
st.set_page_config(
    page_title="MassQL - Mass Spec Data Search",
    page_icon="🔬",
    layout="wide"
)

st.title("🔬 MassQL - Mass Spectrometry Data Search")

# Sidebar for file selection
st.sidebar.header("Data Source Configuration")
st.sidebar.markdown("""(Minimize Browser Window to View System Dialogs)""")
# GitHub ZIP URL input
zip_url = st.sidebar.text_input(
    "Sample ZIP URL",
    value="https://github.com/rithi-krishnaraj/ms2_App/blob/main/mgf_mzmlSample.zip?raw=true",
)

# Prepare the ZIP bytes once user clicks
zip_bytes = None
if st.sidebar.button("Fetch Sample ZIP"):
    try:
        with st.spinner("Downloading ZIP from GitHub..."):
            resp = requests.get(zip_url, timeout=60)
            resp.raise_for_status()
            zip_bytes = resp.content
        st.sidebar.success("ZIP Found. You can now download it.")
    except Exception as e:
        st.sidebar.error(f"Failed to fetch ZIP: {e}")

# Download button uses bytes or a BytesIO buffer
if zip_bytes:
    st.sidebar.download_button(
        label="📥 Download Sample Data Files",
        data=io.BytesIO(zip_bytes),  # or pass zip_bytes directly
        file_name="sample_data.zip",
        mime="application/zip"
    )

# File uploader for multiple files
uploaded_files = st.sidebar.file_uploader(
    "Upload mass spec files",
    type=['mzml', 'mzxml', 'mgf', 'msp'],
    accept_multiple_files=True,
    help="Select one or more mass spectrometry data files"
)

# Process uploaded files
if uploaded_files:
    st.sidebar.success(f"📁 {len(uploaded_files)} file(s) uploaded")
    
    # Show uploaded file names
    with st.sidebar.expander("View uploaded files"):
        for file in uploaded_files:
            st.write(f"✓ {file.name}")
    
    # Store uploaded files for processing
    selected_files = uploaded_files
else:
    selected_files = []

# Main content area
st.header("MassQL Query")
# Information section
with st.expander("ℹ️ About MassQL"):
    st.markdown("""
    **MassQL** is a query language for mass spectrometry data that allows you to:
    
    - Search for specific ions and fragments
    - Filter spectra by retention time, m/z values, and intensity
    - Perform complex queries across multiple data files
    - Extract and analyze specific spectral features
    
    **Supported File Formats:**
    - mzML
    - mzXML
    - MGF
    - MSP
    
    **Learn More:** [MassQL Documentation](https://mwang87.github.io/MassQueryLanguage_Documentation/)
    """)

# Query input
query = st.text_area(
    "Enter your MassQL query",
    value="QUERY scaninfo(MS2DATA)",
    height=150,
    help="Enter a MassQL query to search your mass spec data"
)

# Example queries
with st.expander("📖 Example Queries"):
    st.code("""
# Find MS2 spectra with a specific product ion
QUERY scaninfo(MS2DATA) WHERE MS2PROD=226.1:TOLERANCEPPM=10

# Find precursor ions in a mass range
QUERY scaninfo(MS2DATA) WHERE MS2PREC=100-500

# Find spectra with retention time filter
QUERY scaninfo(MS2DATA) WHERE RT=1-5 AND MS2PROD=100:TOLERANCEPPM=20

# Complex query with multiple conditions
QUERY scaninfo(MS2DATA) WHERE MS2PROD=226.1:TOLERANCEPPM=10 AND MS2PREC=300-400
    """, language="sql")

# Search button
col1, col2, col3 = st.columns([1, 1, 2])

with col1:
    search_button = st.button("🔍 Run Query", type="primary", use_container_width=True)

with col2:
    clear_button = st.button("🗑️ Clear Results", use_container_width=True)

# Results section
if search_button and selected_files and query:
    st.header("Search Results")
    
    with st.spinner("Searching mass spec data..."):
        try:
            results_container = st.container()
            
            for uploaded_file in selected_files:
                with results_container:
                    st.subheader(f"📄 {uploaded_file.name}")
                    
                    try:
                        # Save uploaded file to temporary location
                        temp_path = Path(f"temp_{uploaded_file.name}")
                        temp_path.write_bytes(uploaded_file.read())
                        
                        # Execute query directly without parsing separately
                        results = msql_engine.process_query(
                            query,  # Pass query string directly
                            str(temp_path)
                        )
                        
                        # Check if results exist (results is already a DataFrame)
                        if results is not None and not results.empty:
                            st.success(f"Found {len(results)} matching spectra")
                            
                            # Display results in a dataframe
                            st.dataframe(
                                results,
                                use_container_width=True,
                                height=300, 
                                hide_index=True
                            )
                            
                            # Download button for results
                            csv = results.to_csv(index=False)
                            st.download_button(
                                label="📥 Download Results as CSV",
                                data=csv,
                                file_name=f"massql_results_{uploaded_file.name}.csv",
                                mime="text/csv"
                            )
                        else:
                            st.warning("No matching spectra found")
                        
                        # Clean up temporary file
                        if temp_path.exists():
                            temp_path.unlink()
                            
                    except Exception as e:
                        error_msg = str(e)
                        
                        # Parse the error for better user experience
                        if "UnexpectedCharacters" in error_msg or "No terminal matches" in error_msg:
                            st.error("""❌ **Query Syntax Error**                          
                            Common issues:
                            - Missing `WHERE` or `FILTER` keyword after the function
                            - Incorrect condition syntax
                            - Invalid function names

                            **Correct format:** `QUERY scaninfo(MS2DATA) WHERE MS2PROD=226.1:TOLERANCEPPM=10`

                            Please check the example queries above for correct syntax.""")
                        
                        elif "process_query" in error_msg or "parse" in error_msg:
                            st.error(f"""❌ **File Processing Error**
                            
                            Could not process the file: `{uploaded_file.name}`

                            Possible reasons:
                            - File format not supported or corrupted
                            - Insufficient permissions
                            - File is currently open in another program""")
                                                
                        else:
                            st.error(f"❌ **Error processing file**: {uploaded_file.name}\n\nPlease verify your query syntax and file format.")
                        
                        # Clean up temporary file on error
                        if 'temp_path' in locals() and temp_path.exists():
                            temp_path.unlink()
                    
                    st.divider()
                    
        except Exception as e:
            st.error("❌ **Query Execution Error**")
            st.markdown("""
            An error occurred while executing your query.
            
            Please:
            1. Check your query syntax using the examples below
            2. Verify that your data files are valid and accessible
            3. Ensure you have the correct file format (mzML, mzXML, MGF, or MSP)
            """)

elif search_button and not selected_files:
    st.warning("⚠️ Please select at least one data file to search.")

elif search_button and not query:
    st.warning("⚠️ Please enter a MassQL query.")



