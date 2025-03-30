# sangerqc: Quality Control Tool for Sanger Sequencing Data

## Overview
sangerqc is a web application built using Streamlit that provides a quality control tool for Sanger sequencing data. The application allows users to upload ABI files, set quality score thresholds, and generate modified FASTA files with low-quality bases replaced by a specified character. Users can download the processed results in a ZIP file.

## Features
- Upload multiple ABI files (.ab1) for processing.
- Set a quality score threshold to filter low-quality bases.
- Specify a replacement character for low-quality bases.
- Download processed FASTA files in a ZIP format.

## Usage

### Run without installation
You can use the application without installation by visiting the following link:  
[sangerqc on Streamlit](https://sangerqc.streamlit.app/)

### Run locally
1. Run the application:
   ```
   streamlit run app.py
   ```

2. Open your web browser and go to `http://localhost:8501`.

3. Follow the on-screen instructions to upload ABI files, set quality thresholds, and process the files.

### Running on Streamlit Community Cloud

## Contact
For any questions or support, please contact: [yamaguchi.keito.y0@gmail.com](mailto:yamaguchi.keito.y0@gmail.com)