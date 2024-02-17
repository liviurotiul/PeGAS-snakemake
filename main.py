#==============================================================================================
# This is the entry point for the app
#==============================================================================================

from tqdm import tqdm
from genericpath import isfile
import os
import tempfile
import sys
import pyfiglet
import time

from app.utils import *
from app.extract_data import *
from app.paths import *
from app.report import *
from app.messages import *

import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=UserWarning)

import os
import tempfile

def process_sample_data(data, file_names_dict):
    """
    Process sample data by performing quality control, assembly, MLST, and Abricate analysis.

    Args:
        data (str): Name of the sample data.
        file_names_dict (dict): Dictionary mapping sample data to file names.
        RAW_DATA_PATH (str): Path to the raw data folder.
        DATABASE_PATH (str): Path to the database.

    Returns:
        str: Path to the added sample in the database.
    """

    # Work in a temporary folder
    with tempfile.TemporaryDirectory() as tempdirname:

        # Copy the files to the temp folder
        os.system(f"cp {RAW_DATA_PATH}/{file_names_dict[data][0]} {tempdirname}")
        os.system(f"cp {RAW_DATA_PATH}/{file_names_dict[data][1]} {tempdirname}")

        # WARNING: the names of the fastq files should be identical except for the R1/R2
        R1_path = os.path.join(tempdirname, os.listdir(tempdirname)[0])
        R2_path = os.path.join(tempdirname, os.listdir(tempdirname)[1])

        # Perform quality control before assembly
        fastqc_run_fastq(tempdirname, R1_path, R2_path)

        # Run assembly
        if shovill_run(tempdirname, R1_path, R2_path):
            return None  # Sometimes spades fails, so we return None to indicate an interrupted run

        contigs_path = os.path.join(tempdirname, "shovill", "contigs.fasta")

        # Run MLST
        mlst_run(tempdirname, contigs_path)

        # Run Abricate
        abricate_run(tempdirname, contigs_path)

        # Add to database
        path_to_sample = db_add(tempdirname, data)

    return path_to_sample

tqdm.write(pyfiglet.figlet_format("PeGAS"))
tqdm.write("")
tqdm.write("")

paths = [DATABASE_PATH, RAW_DATA_PATH, TABLES_PATH, PANGENOME_PATH, WEB_PATH, APP_PATH]
for p in paths:
    if not os.path.exists(p):
        tqdm.write(f"{p} was not found, please check if the folder structure has been altered or if any of\
            the {p} has been accidentely deleted, Exiting")
        raise FileNotFoundError

#==============================================================================================
# SETUP
#==============================================================================================

db_entries = find_all_processed_samples(DATABASE_PATH)

new_entries, file_names_dict = get_new_entries(RAW_DATA_PATH)

#==============================================================================================
# Begin the pipline
#==============================================================================================

tqdm.write("")
tqdm.write(BEGIN)

start = time.time()
progress_bar = tqdm(new_entries)

processed_samples = 0
past_processed_samples = len(db_entries)
import pdb;pdb.set_trace
for data in progress_bar:

    tqdm.write(f"Processing {data}")
    progress_bar.set_description(f"Processing {data}")

    if db_find_shovill(DATABASE_PATH, data) == "":

        # # Work in a temporary folder
        # with tempfile.TemporaryDirectory() as tempdirname:

        #     # Copy the files to the temp folder
        #     os.system(f"cp {RAW_DATA_PATH}/{file_names_dict[data][0]} {tempdirname}")
        #     os.system(f"cp {RAW_DATA_PATH}/{file_names_dict[data][1]} {tempdirname}")

        #     # WARNING: the names of the fastq files sould be identical except for the R1/R2
        #     R1_path = os.path.join(tempdirname, os.listdir(tempdirname)[0])
        #     R2_path = os.path.join(tempdirname, os.listdir(tempdirname)[1])

        #     # Perform quality control b4 assembly
        #     progress_bar.set_description(f"Processing {data} - FastQC")
        #     fastqc_run_fastq(tempdirname, R1_path, R2_path)

        #     # Run assembly
        #     progress_bar.set_description(f"Processing {data} - Shovill")
        #     if shovill_run(tempdirname, R1_path, R2_path):
        #         continue # Sometimes spades fails due to unknown reasons so I made this chekc to make sure the run is not interrupted

        #     contigs_path = os.path.join(tempdirname, "shovill", "contigs.fasta")

        #     # Run MLST
        #     progress_bar.set_description(f"Processing {data} - MLST")
        #     mlst_run(tempdirname, contigs_path)

        #     # Run abricate
        #     progress_bar.set_description(f"Processing {data} - Abricate")
        #     abricate_run(tempdirname, contigs_path)
            
        #     # Add to db
        #     path_to_sample = db_add(tempdirname, data)
        process_sample_data(data, file_names_dict)
        processed_samples += 1

    # Run Prokka
    if not db_find_prokka(DATABASE_PATH, data):

        # tqdm.write(f"prokka folder not found for {data}")

        progress_bar.set_description(f"Processing {data} - Prokka")
        prokka_run(f"{DATABASE_PATH}/{db_find_entry_folder(data)}")

tqdm.write(f"Number of processed samples this run: {processed_samples}")
tqdm.write(f"Number of samples already found in the database: {past_processed_samples}")
tqdm.write(f"Number of samples that were not processed: {len(progress_bar) - processed_samples}")

# Pangenome analysis
collect_gffs()

roary_run(PANGENOME_PATH)

# Write report
tqdm.write(GEN_REPORT)
tables = generate_tables()

generate_report(tables)

end = time.time()
tqdm.write(f"Time elapsed: {end - start}")