import streamlit as st
import pandas as pd
import json
import time
import os
import io
import zipfile
from openai import OpenAI
from datetime import datetime
from Bio import Entrez
import requests
import pdf2image
from PIL import Image
import pytesseract
from pytesseract import Output, TesseractError

# Streamlit UI Configuration
st.set_page_config(
    page_title="Clinical Research Analysis Tool",
    page_icon=":mag:",
    layout="wide",
    initial_sidebar_state="expanded",
)

# Helper functions for PDF processing
def convert_pdf_to_txt_file(pdf_file):
    """Convert PDF to a single text file"""
    import PyPDF2
    pdf_reader = PyPDF2.PdfReader(pdf_file)
    text = ""
    for page in pdf_reader.pages:
        text += page.extract_text() + "\n\n"
    return text, len(pdf_reader.pages)

def convert_pdf_to_txt_pages(pdf_file):
    """Convert PDF to multiple text pages"""
    import PyPDF2
    pdf_reader = PyPDF2.PdfReader(pdf_file)
    text_pages = []
    for page in pdf_reader.pages:
        text_pages.append(page.extract_text())
    return text_pages, len(pdf_reader.pages)

def save_pages(text_data):
    """Save text pages to a ZIP file"""
    zip_path = "pdf_to_txt.zip"
    with zipfile.ZipFile(zip_path, 'w') as zipf:
        for i, text in enumerate(text_data):
            zipf.writestr(f"page_{i+1}.txt", text)
    return zip_path

def displayPDF(file):
    """Display PDF in Streamlit"""
    # Create a temporary PDF file
    with open("temp.pdf", "wb") as f:
        f.write(file)
    
    # Display the PDF using an iframe
    pdf_display = f"""
        <iframe src="temp.pdf" width="700" height="1000" type="application/pdf"></iframe>
    """
    st.markdown(pdf_display, unsafe_allow_html=True)

def images_to_txt(pdf_file, lang='eng'):
    """Convert PDF to text using OCR"""
    images = pdf2image.convert_from_bytes(pdf_file)
    texts = []
    for image in images:
        text = pytesseract.image_to_string(image, lang=lang)
        texts.append(text)
    return texts, len(images)

# State Management
if 'analysis_results' not in st.session_state:
    st.session_state['analysis_results'] = []
if 'progress' not in st.session_state:
    st.session_state['progress'] = 0
if 'total_papers' not in st.session_state:
    st.session_state['total_papers'] = 0
if 'search_completed' not in st.session_state:
    st.session_state['search_completed'] = False
if 'openai_api_key' not in st.session_state:
    st.session_state['openai_api_key'] = None
if 'api_key_valid' not in st.session_state:
    st.session_state['api_key_valid'] = False
if 'df' not in st.session_state:
    st.session_state['df'] = pd.DataFrame()
if 'pdf_texts' not in st.session_state:
    st.session_state['pdf_texts'] = []
if 'pdf_analysis_completed' not in st.session_state:
    st.session_state['pdf_analysis_completed'] = False
if 'show_clear_confirmation' not in st.session_state:
    st.session_state['show_clear_confirmation'] = False
if 'show_new_search_dialog' not in st.session_state:
    st.session_state['show_new_search_dialog'] = False
if 'search_action' not in st.session_state:
    st.session_state['search_action'] = None

# Function to reset the app state
def reset_app_state():
    st.session_state['analysis_results'] = []
    st.session_state['progress'] = 0
    st.session_state['total_papers'] = 0
    st.session_state['search_completed'] = False
    st.session_state['pdf_analysis_completed'] = False
    st.session_state['df'] = pd.DataFrame()
    st.session_state['pdf_texts'] = []
    st.session_state['show_clear_confirmation'] = False
    st.session_state['show_new_search_dialog'] = False
    st.session_state['search_action'] = None

# Sidebar for Inputs
st.sidebar.header("Clinical Research Analysis Tool")

# Tab selection
tab_selection = st.sidebar.radio("Select Source", ["PubMed Search", "PDF Upload"])

# OpenAI API Key Input in Sidebar
with st.sidebar:
    # Only show the API key input if the key is not yet validated
    if not st.session_state['api_key_valid']:
        openai_api_key = st.text_input("OpenAI API Key", type="password")

        if openai_api_key:
            st.session_state['openai_api_key'] = openai_api_key

            # Validate API Key
            try:
                headers = {
                    'accept': 'application/json',
                    'x-baychatgpt-accesstoken': openai_api_key
                }
                response = requests.get('https://chat.int.bayer.com/api/v2/users/me', headers=headers)
                response.raise_for_status()
                st.session_state['api_key_valid'] = True
                st.success("API Key Validated! :white_check_mark:")
            except requests.exceptions.RequestException as e:
                st.error(f"API Key Invalid: {e}")
                st.session_state['api_key_valid'] = False

        elif st.session_state['openai_api_key'] is None:
            st.info("Please add your OpenAI API key to continue.", icon="🗝️")
    else:
        # If API key is already validated, show a success message instead of the input box
        st.success("API Key Authenticated ✓")
        if st.button("Change API Key"):
            st.session_state['api_key_valid'] = False
            st.session_state['openai_api_key'] = None
            st.rerun()

    # Disable the analysis buttons if the API key is not provided or invalid
    start_analysis_disabled = not (st.session_state['openai_api_key'] and st.session_state['api_key_valid'])
    
    # Add a clear table button in the sidebar
    st.sidebar.markdown("---")
    if st.sidebar.button("Clear Results Table"):
        st.session_state['show_clear_confirmation'] = True

# Show confirmation popup for clearing the table
if st.session_state['show_clear_confirmation']:
    with st.sidebar:
        st.warning("⚠️ Are you sure you want to clear all results?")
        col1, col2 = st.columns(2)
        if col1.button("Yes, Clear"):
            reset_app_state()
            st.rerun()
        if col2.button("Cancel"):
            st.session_state['show_clear_confirmation'] = False
            st.rerun()

# Configure NCBI Search
if "ncbi_email" in st.secrets:
    Entrez.email = st.secrets["ncbi_email"]
if "PM_Key" in st.secrets:
    Entrez.api_key = st.secrets["PM_Key"]

# Functions
def search_and_fetch_pubmed(query, max_results):
    """Search PubMed and fetch details in one function."""
    with Entrez.esearch(db="pubmed", term=query, retmax=max_results) as handle:
        record = Entrez.read(handle)

    id_list = record["IdList"]
    if not id_list:
        return []

    ids = ','.join(id_list)
    with Entrez.efetch(db="pubmed", id=ids, retmode="xml") as handle:
        records = Entrez.read(handle)

    return records

def add_to_excel(output_list, excel_file_path):
    """Add a list of dictionaries to an Excel sheet."""
    df = pd.DataFrame(output_list)
    df.to_excel(excel_file_path, index=False)
    return excel_file_path

def analyze_paper(paper, openai_api_key, is_pdf=False):
    """Call OpenAI to analyze the paper and return the results."""
    client = OpenAI(api_key=openai_api_key, base_url="https://chat.int.bayer.com/api/v1")

    system_prompt = """
    You are a bot speaking with another program that takes JSON formatted text as an input. Only return results in JSON format, with NO PREAMBLE.
    The user will input the results from a PubMed search or a full-text clinical trial PDF. Your job is to extract the exact information to return:
      'Title': The complete article title
      'PMID': The Pubmed ID of the article (if available, otherwise 'NA')
      'Full Text Link' : If available, the DOI URL, otherwise, NA
      'Subject of Study': The type of subject in the study. Human, Animal, In-Vitro, Other
      'Disease State': Disease state studied, if any, or if the study is done on a healthy population. leave blank if disease state or healthy patients is not mentioned explicitly.
      'Number of Subjects Studied': If human, the total study population. Otherwise, leave blank. This field needs to be an integer or empty.
      'Type of Study': Type of study done. Can be '3. RCT' for randomized controlled trial, '1. Meta-analysis','2. Systematic Review','4. Cohort Study', or '5. Other'. If it is '5. Other', append a short description
      'Study Design': Brief and succinct details about study design, if applicable
      'Intervention': Intervention(s) studied, if any. Intervention is the treatment applied to the group.
      'Intervention Dose': Go in detail here about the intervention's doses and treatment duration if available.
      'Intervention Dosage Form': A brief description of the dosage form - ie. oral, topical, intranasal, if available.
      'Control': Control or comarators, if any
      'Primary Endpoint': What the primary endpoint of the study was, if available. Include how it was measured too if available.
      'Primary Endpoint Result': The measurement for the primary endpoints
      'Secondary Endpoints' If available
      'Safety Endpoints' If available
      'Results Available': Yes or No
      'Primary Endpoint Met': Summarize from results whether or not the primary endpoint(s) was met: Yes or No or NA if results unavailable
      'Statistical Significance': alpha-level and p-value for primary endpoint(s), if available
      'Clinical Significance': Effect size, and Number needed to treat (NNT)/Number needed to harm (NNH), if available
      'Main Author': Last name, First initials
      'Other Authors': Last name, First initials; Last name First initials; ...
      'Journal Name': Full journal name
      'Date of Publication': YYYY-MM-DD
      'Error': Error description, if any. Otherwise, leave emtpy
    """

    if is_pdf:
        system_prompt += """
        Note: This is a full-text PDF of a clinical trial. Extract as much detail as possible from the full text.
        """

    conversation = client.chat.completions.create(
        model='gpt-4o-mini',  # Use appropriate model
        messages=[
            {
                "role": "system",
                "content": system_prompt
            },
            {
                "role": "user",
                "content": str(paper)
            }
        ]
    )

    cleaned_str = conversation.choices[0].message.content.replace("```json", "").replace("```", "").strip()
    return json.loads(cleaned_str)

# Function to perform PubMed search and analysis
def perform_pubmed_analysis(query, max_results, action="new"):
    with st.spinner(f"Searching PubMed for '{query}'..."):
        papers = search_and_fetch_pubmed(query, max_results)
        if 'PubmedArticle' in papers:
            st.session_state['total_papers'] = len(papers['PubmedArticle'])
            st.write(f"Found {st.session_state['total_papers']} papers.")
        else:
            st.error("No papers found. Try a different search query.")
            st.session_state['total_papers'] = 0
            return

    if st.session_state['total_papers'] > 0:
        progress_bar = st.progress(0)
        
        # Initialize or append to results based on action
        if action == "new":
            st.session_state['analysis_results'] = []
        
        for i, paper in enumerate(papers['PubmedArticle']):
            try:
                with st.spinner(f"Analyzing paper {i+1}/{st.session_state['total_papers']}..."):
                    result = analyze_paper(paper, st.session_state['openai_api_key'])
                    st.session_state['analysis_results'].append(result)
            except Exception as e:
                st.error(f"Error analyzing paper {i+1}: {e}")

            st.session_state['progress'] = (i + 1) / st.session_state['total_papers']
            progress_bar.progress(st.session_state['progress'])

        st.session_state['search_completed'] = True

# Function to perform PDF analysis
def perform_pdf_analysis(pdf_files, ocr_box=False, language_option="English", action="new"):
    languages = {
        'English': 'eng',
        'French': 'fra',
        'Arabic': 'ara',
        'Spanish': 'spa',
    }
    
    st.session_state['total_papers'] = len(pdf_files)
    progress_bar = st.progress(0)
    
    # Initialize or append to results based on action
    if action == "new":
        st.session_state['pdf_texts'] = []
        st.session_state['analysis_results'] = []
    
    for i, pdf_file in enumerate(pdf_files):
        try:
            with st.spinner(f"Processing file {i+1}/{len(pdf_files)}: {pdf_file.name}"):
                # Extract text from PDF
                file_extension = pdf_file.name.split(".")[-1]
                path = pdf_file.read()
                
                if file_extension == "pdf":
                    if ocr_box:
                        texts, _ = images_to_txt(path, languages[language_option])
                        text_content = "\n\n".join(texts)
                    else:
                        text_content, _ = convert_pdf_to_txt_file(io.BytesIO(path))
                else:  # Image files
                    pil_image = Image.open(io.BytesIO(path))
                    text_content = pytesseract.image_to_string(pil_image, lang=languages[language_option])
                
                # Store extracted text
                st.session_state['pdf_texts'].append({
                    'filename': pdf_file.name,
                    'content': text_content
                })
                
                # Analyze the PDF content
                with st.spinner(f"Analyzing content of {pdf_file.name}..."):
                    result = analyze_paper(text_content, st.session_state['openai_api_key'], is_pdf=True)
                    # Add filename to result
                    result['Filename'] = pdf_file.name
                    st.session_state['analysis_results'].append(result)
        
        except Exception as e:
            st.error(f"Error processing {pdf_file.name}: {e}")
        
        st.session_state['progress'] = (i + 1) / len(pdf_files)
        progress_bar.progress(st.session_state['progress'])
    
    st.session_state['pdf_analysis_completed'] = True
    st.session_state['search_completed'] = True

# Main UI based on selected tab
if tab_selection == "PubMed Search":
    st.title("PubMed Clinical Trial Analysis")
    
    # PubMed search parameters
    query = st.text_input("Search Query", "loratadine") + " AND (clinicaltrial[filter]"
    max_results = st.number_input("Max Results", min_value=1, max_value=400, value=20, step=1)
    
    # Configure file name
    query_word = query.split()[0]  # Extract the first word from the query
    filename = f"Results {query_word} {time.strftime('%y%m%d')}.xlsx"
    
    # Start PubMed analysis button
    if st.button("Start PubMed Analysis", disabled=start_analysis_disabled):
        # Check if there are existing results
        if st.session_state['analysis_results'] and len(st.session_state['analysis_results']) > 0:
            st.session_state['show_new_search_dialog'] = True
        else:
            # No existing results, proceed with new search
            perform_pubmed_analysis(query, max_results, action="new")

elif tab_selection == "PDF Upload":
    st.title("Clinical Trial PDF Analysis")
    
    # PDF upload parameters
    languages = {
        'English': 'eng',
        'French': 'fra',
        'Arabic': 'ara',
        'Spanish': 'spa',
    }
    
    # OCR options
    ocr_box = st.checkbox('Enable OCR (for scanned documents)')
    if ocr_box:
        language_option = st.selectbox('Select the document language', list(languages.keys()))
    else:
        language_option = "English"
    
    # Multiple file uploader
    pdf_files = st.file_uploader("Upload Clinical Trial PDF(s)", type=['pdf', 'png', 'jpg'], accept_multiple_files=True)
    
    if pdf_files:
        st.write(f"Uploaded {len(pdf_files)} file(s)")
        
        # Process PDFs button
        if st.button("Process and Analyze PDFs", disabled=start_analysis_disabled):
            # Check if there are existing results
            if st.session_state['analysis_results'] and len(st.session_state['analysis_results']) > 0:
                st.session_state['show_new_search_dialog'] = True
            else:
                # No existing results, proceed with new analysis
                perform_pdf_analysis(pdf_files, ocr_box, language_option, action="new")

# Show dialog for new search when results already exist
if st.session_state['show_new_search_dialog']:
    dialog_container = st.container()
    with dialog_container:
        st.warning("### Existing results found")
        st.write("Would you like to start a new search with a fresh table or add to the existing table?")
        
        col1, col2, col3 = st.columns([2, 2, 1])
        
        if col1.button("Start Fresh", key="fresh_button"):
            st.session_state['search_action'] = "new"
            st.session_state['show_new_search_dialog'] = False
            
            # Perform the appropriate action based on the selected tab
            if tab_selection == "PubMed Search":
                perform_pubmed_analysis(query, max_results, action="new")
            else:  # PDF Upload
                perform_pdf_analysis(pdf_files, ocr_box, language_option, action="new")
            
            # Force a rerun to refresh the UI and close the dialog
            st.rerun()
        
        if col2.button("Add to Existing", key="add_button"):
            st.session_state['search_action'] = "append"
            st.session_state['show_new_search_dialog'] = False
            
            # Perform the appropriate action based on the selected tab
            if tab_selection == "PubMed Search":
                perform_pubmed_analysis(query, max_results, action="append")
            else:  # PDF Upload
                perform_pdf_analysis(pdf_files, ocr_box, language_option, action="append")
            
            # Force a rerun to refresh the UI and close the dialog
            st.rerun()
            
        if col3.button("Cancel", key="cancel_button"):
            st.session_state['show_new_search_dialog'] = False
            st.rerun()

# Display Results and Download Button (common for both tabs)
if st.session_state['search_completed']:
    st.success("Analysis Complete!")

    if st.session_state['analysis_results']:
        df = pd.DataFrame(st.session_state['analysis_results'])

        # Convert lists to strings in problematic columns
        for col in df.columns:
            if df[col].dtype == 'object':
                try:
                    df[col] = df[col].astype(str)
                except:
                    print(f"can't convert {col}")
                    pass
        
        # Add a selection column to the dataframe
        df.insert(0, " ", False)
        
        # Store the dataframe in session state
        st.session_state['df'] = df
        
        # Display the interactive table with row selection
        st.write("### Select papers to download:")
        
        # Use the data editor with a checkbox column
        edited_df = st.data_editor(
            df,
            use_container_width=True,
            num_rows="fixed",
            hide_index=True,
            key="results_table"
        )
        
        # Get the selected rows
        selected_df = edited_df[edited_df[" "] == True] if " " in edited_df.columns else pd.DataFrame()
        num_selected = len(selected_df)
        st.write(f"Selected {num_selected} rows")
        
        # File naming based on source
        if tab_selection == "PubMed Search":
            filename = f"PubMed_{query_word}_{time.strftime('%y%m%d')}.xlsx"
        else:
            filename = f"PDF_Analysis_{time.strftime('%y%m%d')}.xlsx"
        
        # Download buttons
        col1, col2 = st.columns(2)
        
        # Download all results
        df_to_save = df.drop(columns=[" "])
        excel_file_all = add_to_excel(df_to_save.to_dict('records'), filename)
        with open(excel_file_all, "rb") as f:
            col1.download_button(
                label="Download All Results",
                data=f,
                file_name=filename,
                mime="application/vnd.ms-excel",
            )
        
        # Download selected rows
        if num_selected > 0:
            selected_filename = f"Selected_{time.strftime('%y%m%d')}.xlsx"
            selected_excel_file = f"Selected_{filename}"
            selected_df_to_save = selected_df.drop(columns=[" "])
            selected_df_to_save.to_excel(selected_excel_file, index=False)
            
            with open(selected_excel_file, "rb") as f:
                col2.download_button(
                    label=f"Download Selected ({num_selected}) Rows",
                    data=f,
                    file_name=selected_filename,
                    mime="application/vnd.ms-excel",
                )
        else:
            col2.write("Select rows to enable partial download")
            
        # If PDF analysis was done, offer text download
        if tab_selection == "PDF Upload" and st.session_state['pdf_analysis_completed']:
            st.write("### Download Extracted Text:")
            
            for i, pdf_text in enumerate(st.session_state['pdf_texts']):
                col1, col2 = st.columns([3, 1])
                col1.write(f"{i+1}. {pdf_text['filename']}")
                col2.download_button(
                    label="Download Text",
                    data=pdf_text['content'],
                    file_name=f"{pdf_text['filename']}.txt",
                    mime="text/plain",
                    key=f"download_text_{i}"
                )
    else:
        st.info("No results to display.")

# Footer with styling
st.markdown("""
<style>
.footer {
    position: fixed;
    left: 0;
    bottom: 0;
    width: 100%;
    background-color: rgba(240, 242, 246, 0.9);
    color: #262730;
    text-align: center;
    padding: 10px;
    font-size: 14px;
}
</style>
""", unsafe_allow_html=True)

# Hide Streamlit branding
hide = """
<style>
footer {visibility: hidden;}
.viewerBadge_container__1QSob {visibility: hidden;}
#MainMenu {visibility: hidden;}
</style>
"""
st.markdown(hide, unsafe_allow_html=True)