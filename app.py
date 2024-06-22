import os
import streamlit as st
import zipfile
from io import BytesIO
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from tempfile import TemporaryDirectory

def abi_to_fasta(abi_file, output_file, quality_threshold=20):
    try:
        record = SeqIO.read(abi_file, "abi")
    except Exception as e:
        return f"Error reading {abi_file}: {e}"

    modified_seq = ""
    for letter, score in zip(record.seq, record.letter_annotations["phred_quality"]):
        if score < quality_threshold:
            modified_seq += "-"
        else:
            modified_seq += letter

    try:
        modified_record = SeqRecord(Seq(modified_seq), id=record.id, description="modified sequence")
        SeqIO.write(modified_record, output_file, "fasta")
    except Exception as e:
        return f"Error writing to {output_file}: {e}"

    return None

st.title("Sanger-QC: A quality control tool for sanger sequencing data.")

st.markdown("""
## How to Use
1. **Upload ABI files**:
   - Drag and drop or click the "Browse files" button and select the ABI files (.ab1) you want to process. You can upload multiple files at once.
2. **Set Quality Score Threshold**:
   - Use the "Quality Score Threshold" slider to set the quality score threshold. The default value is 20.
3. **Choose Replacement Character**:
   - Select whether to replace low-quality bases with `-` or `N`.
4. **Start Processing**:
   - Click the "Process" button to start processing all uploaded files. Sequences with quality scores below the threshold will be replaced with the chosen character in the resulting FASTA files.
5. **Download Results**:
   - Once processing is complete, a button to download a ZIP file containing the resulting FASTA files will appear. Click this button to download the ZIP file.
""")

uploaded_files = st.file_uploader("Upload ABI files", type=["ab1"], accept_multiple_files=True)
quality_threshold = st.slider("Quality Score Threshold", min_value=0, max_value=50, value=20)
replacement_char = st.radio("Replacement Character for Low-Quality Bases", ["ー", "N"])

if st.button("Process"):
    if not uploaded_files:
        st.error("Please upload at least one ABI file.")
    else:
        with TemporaryDirectory() as temp_dir:
            output_files = []
            error_messages = []
            for uploaded_file in uploaded_files:
                abi_path = os.path.join(temp_dir, uploaded_file.name)
                with open(abi_path, "wb") as f:
                    f.write(uploaded_file.getbuffer())

                output_file = os.path.join(temp_dir, os.path.splitext(uploaded_file.name)[0] + "_modified.fasta")
                error_message = abi_to_fasta(abi_path, output_file, quality_threshold, replacement_char)
                if error_message:
                    error_messages.append(error_message)
                else:
                    output_files.append(output_file)

            if error_messages:
                for error_message in error_messages:
                    st.error(error_message)
            else:
                # Create a zip file of the output files
                zip_buffer = BytesIO()
                with zipfile.ZipFile(zip_buffer, "w") as zip_file:
                    for output_file in output_files:
                        zip_file.write(output_file, os.path.basename(output_file))
                zip_buffer.seek(0)

                st.success("All files processed.")
                st.download_button(
                    label="Download All Results (ZIP)",
                    data=zip_buffer,
                    file_name="filtered_sequences.zip",
                    mime="application/zip"
                )