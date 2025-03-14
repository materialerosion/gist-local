# PubMed Analysis Tool

This Streamlit application allows you to search PubMed for relevant research articles and analyze them using OpenAI. It extracts key information from the abstracts and presents the results in a structured format, which can then be downloaded as an Excel file.

## Features

*   **PubMed Search:**  Search PubMed using a custom query with a built-in `AND (clinicaltrial[filter])` to target clinical trials.
*   **OpenAI Analysis:** Leverages OpenAI's language models to extract structured information from the retrieved abstracts.
*   **Structured Output:**  Displays extracted information in a table format, including:
    *   Title
    *   PMID
    *   Full Text Link
    *   Subject of Study
    *   Disease State
    *   Number of Subjects Studied
    *   Type of Study
    *   Study Design
    *   Intervention
    *   Intervention Dose
    *   Intervention Dosage Form
    *   Control
    *   Primary Endpoint
    *   Primary Endpoint Result
    *   Secondary Endpoints
    *   Safety Endpoints
    *   Results Available
    *   Primary Endpoint Met
    *   Statistical Significance
    *   Clinical Significance
    *   Main Author
    *   Other Authors
    *   Journal Name
    *   Date of Publication
    *   Error
*   **Downloadable Results:**  Download the analysis results as an Excel (`.xlsx`) file.
*   **API Key Validation:** Validates your OpenAI API key before starting the analysis.
*   **Progress Tracking:**  Provides a progress bar to indicate the status of the search and analysis.

## Prerequisites

*   **Python 3.7+**
*   **Streamlit:** `pip install streamlit`
*   **Pandas:** `pip install pandas`
*   **Biopython:** `pip install biopython`
*   **OpenAI Python Library:** `pip install openai`
*   **Requests:** `pip install requests`

## Setup

1.  **Clone the repository:**

    ```bash
    git clone <repository_url>
    cd <repository_directory>
    ```

2.  **Install the required packages:**

    ```bash
    pip install streamlit pandas biopython openai requests
    ```

3.  **Configure NCBI API Key:**

    *   Create a file named `.streamlit/secrets.toml` in your project directory.
    *   Add your NCBI email and API key to the file:

        ```toml
        ncbi_email = "your_email@example.com"
        PM_Key = "YOUR_NCBI_API_KEY"
        ```

    *   You can obtain an NCBI API key from [NCBI](https://www.ncbi.nlm.nih.gov/account/settings/). This is optional but highly recommended to avoid rate limiting.

4.  **Configure OpenAI API Key:**

    *   This tool is designed to use a custom OpenAI endpoint. You'll need an API key that works with `https://chat.int.bayer.com/api/v1`.
    *   The application will prompt you for your OpenAI API key in the sidebar.

## Usage

1.  **Run the Streamlit application:**

    ```bash
    streamlit run streamlit_app.py
    ```


2.  **Enter Search Parameters:**

    *   In the sidebar, enter your PubMed search query in the "Search Query" field.  The tool automatically appends `AND (clinicaltrial[filter])` to your query to specifically target clinical trials.
    *   Specify the maximum number of results to retrieve in the "Max Results" field.
    *   Enter your OpenAI API key in the "OpenAI API Key" field.

3.  **Start Analysis:**

    *   Click the "Start Analysis" button.  The button will be disabled until a valid OpenAI API key is provided.

4.  **View Results:**

    *   The application will display a progress bar while searching PubMed and analyzing the abstracts.
    *   Once the analysis is complete, the extracted information will be displayed in a table.

5.  **Download Results:**

    *   Click the "Download Excel File" button to download the results as an Excel file.

## Important Notes

*   **NCBI API Key:**  Using an NCBI API key is highly recommended to avoid rate limiting when searching PubMed.
*   **OpenAI API Key:**  A valid OpenAI API key is required to analyze the abstracts. Ensure your API key has access to the specified endpoint (`https://chat.int.bayer.com/api/v1`). The tool validates the API key before starting the analysis.
*   **Rate Limiting:**  Be mindful of the usage limits for both the NCBI and OpenAI APIs.
*   **Error Handling:**  The application includes basic error handling, but you may encounter errors due to network issues, API rate limits, or invalid input.
*   **JSON Formatting:** The OpenAI model is instructed to return results in JSON format. Any deviations from this format may cause errors.
*   **Modifying the Query:** The tool automatically adds `" AND (clinicaltrial[filter])"` to the end of your query to filter for clinical trials.  Keep this in mind when constructing your search query.

## Troubleshooting

*   **API Key Errors:**  Double-check that your NCBI and OpenAI API keys are correct and properly configured.
*   **Rate Limiting:**  If you encounter rate limiting errors, try reducing the number of results or using an NCBI API key.
*   **Network Errors:**  Check your internet connection.
*   **JSONDecodeError:** This error usually indicates that the OpenAI model did not return a valid JSON response.  This can be due to the query, the model, or temporary issues with the OpenAI API.
*   **No results to display:** Check your query is valid and returns results.

## Disclaimer

This tool is provided as-is, and no warranty is provided. Use at your own risk. The accuracy of the extracted information depends on the quality of the abstracts and the capabilities of the OpenAI model.
