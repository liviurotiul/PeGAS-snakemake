import itertools
import pandas as pd
from tqdm import tqdm
from scripts.utils import species_dict
from scripts.utils import forbidden_files

raw_data_path = config["raw_data"]
# list_of_samples = config["sample_list"]
try:
    list_of_samples = config["sample_list"]
except KeyError:
    list_of_samples = None

if not list_of_samples:

    list_of_samples = os.listdir(raw_data_path)
    R1 = [f.replace("_L001_R1_001.fastq.gz", "") for f in list_of_samples if f.endswith("_L001_R1_001.fastq.gz")]
    R2 = [f.replace("_L001_R2_001.fastq.gz", "") for f in list_of_samples if f.endswith("_L001_R2_001.fastq.gz")]
    list_of_samples = list(set(R1+R2))

else:
    # Read list of samples; it is a text file with one sample name per line
    with open(list_of_samples, 'r') as f:
        list_of_samples = f.read().splitlines()

# Make the list of files by adding the R1 and R2 to the sample name
R1_files = [f"{sample}_L001_R1_001.fastq.gz" for sample in list_of_samples]
R2_files = [f"{sample}_L001_R2_001.fastq.gz" for sample in list_of_samples]
list_of_files = R1_files + R2_files

if not os.path.exists("raw_data"):
    os.system("mkdir raw_data")

# Copy the contents of raw_data_path to the raw_data folder in the working directory
# Skipping the ones that already exist
# Only copy the files that are in the list_of_files
for file in list_of_files:
    if file not in os.listdir("raw_data"):
        os.system(f"cp {raw_data_path}/{file} raw_data")

# Delete files that are not in the list of samples
for file in os.listdir("raw_data"):
    if file not in list_of_files:
        os.system(f"rm -rf raw_data/{file}")

if "results" in os.listdir():
    # Delete files that are not in the list of samples from the results folder
    for folder in os.listdir("results"):
        if folder not in list_of_samples:
            os.system(f"rm -rf results/{folder}")
    
if "pangenome" in os.listdir():
    # Delete files that are not in the list of samples from the pangenome folder 
    for species in os.listdir("pangenome"):
        for file in os.listdir(f"pangenome/{species}"):

            # If the file is not in the list of samples delete it
            if file.replace(".gff", "") not in list_of_samples:
                os.system(f"rm pangenome/{species}/{file}")

                # Should also delete the pangenomic analysis since it is not valid anymore
                if "output" in os.listdir(f"pangenome/{species}"):
                    os.system(f"rm -rf pangenome/{species}/output")

                # If the species folder is left with only one file delete it
                if len(os.listdir(f"pangenome/{species}")) < 2:
                    os.system(f"rm -rf pangenome/{species}")

file_dict = {}
file_set = set()

files = os.listdir("raw_data")
files = [f for f in files if f.endswith('.fastq.gz')]

# Exclude files that are copies of other files
files = [f for f in files if not any([f"({i})" in f for i in range(10)])]

# Process each file
for f in tqdm(files, desc="Processing files"):

    entry_name = f.replace('_R1', '').replace('_R2', '').replace('.fastq.gz', '')
    file_dict.setdefault(entry_name, []).append(f)
    file_set.add(entry_name)

for entry in file_set:
    if len(file_dict[entry]) != 2:
        raise ValueError(f"Entry {entry} does not have a pair of files")

fastq_file_names = [sorted([R1, R2]) for R1, R2 in file_dict.values()]

sample_names = [file[0].replace("_L001_R1_001.fastq.gz", "") for file in fastq_file_names]

rule all:
    input:
        # define fastqc output like data/results/fastqc/sample_name/R*
        fastqc=expand(
            "results/{sample}/fastqc/{sample}_L001_R{read}_001_fastqc.html",
            sample=sample_names,
            read=[1,2]
        ),
        shovill=expand(
            "results/{sample}/shovill/contigs.fa",
            sample=sample_names
        ),
        abricate_ncvi=expand(
            "results/{sample}/abricate_ncbi.tsv",
            sample=sample_names
        ),
        abricate_plasmidfinder=expand(
            "results/{sample}/abricate_plasmidfinder.tsv",
            sample=sample_names
        ),
        abricate_vfdb=expand(
            "results/{sample}/abricate_vfdb.tsv",
            sample=sample_names
        ),
        mlst=expand(
            "results/{sample}/mlst.tsv",
            sample=sample_names
        ),
        dataframe="dataframe/results.csv",
        prokka=expand(
            "results/{sample}/prokka/{sample}.gff",
            sample=sample_names
        ),
        pangenome_flag="flags/.pangenome",
        Box_Contig_Length="report/QC_Box_Plot.html",
        DF_Full_Table="report/Gene_Table.html",
        DF_Reads_Table="report/QC_Table.html",
        Heatmap_Pangenomic="report/Pangenomic_Heatmap.html",
        Heatmap_Plasmids_Full_Figure_Coverage="report/Plasmid_Gene_Heatmap.html",
        Heatmap_Resistance_Full_Figure="report/Resistance_Heatmap.html",
        Heatmap_Virulence_Full_Figure_Coverage="report/Virulence_Heatmap.html",
        Network_Samples_Figure="report/Network_Chart.html",
        Pangenome_Pie_Chart="report/Pangenome_Pie_Chart.html",
        Scatter_Contig_Length="report/QC_Scatter_Plot.html",
        Subtype_HTML_String="report/Subtype_Pie_Charts.html",
        Sunburst_Figure="report/Sunburst_Chart.html",
        flag="flags/.report"

rule fastqc:
    input:
        R1="raw_data/{sample}_L001_R1_001.fastq.gz",
        R2="raw_data/{sample}_L001_R2_001.fastq.gz"
    output:
        R1="results/{sample}/fastqc/{sample}_L001_R1_001_fastqc.html",
        R2="results/{sample}/fastqc/{sample}_L001_R2_001_fastqc.html"
    shell:
        "fastqc {input.R1} {input.R2} -o results/{wildcards.sample}/fastqc"

rule shovill:
    input:
        R1="raw_data/{sample}_L001_R1_001.fastq.gz",
        R2="raw_data/{sample}_L001_R2_001.fastq.gz"
    output:
        assembly="results/{sample}/shovill/contigs.fa"
    shell:
        "shovill --trim --outdir results/{wildcards.sample}/shovill --R1 {input.R1} --R2 {input.R2} --force --cpus 4"

rule abricate_ncbi:
    input:
        assembly="results/{sample}/shovill/contigs.fa"
    output:
        abricate="results/{sample}/abricate_ncbi.tsv"
    shell:
        "abricate --db ncbi {input.assembly} > {output.abricate}"

rule abricate_plasmidfinder:
    input:
        assembly="results/{sample}/shovill/contigs.fa"
    output:
        abricate="results/{sample}/abricate_plasmidfinder.tsv"
    shell:
        "abricate --db plasmidfinder {input.assembly} > {output.abricate}"

rule abricate_vfdb:
    input:
        assembly="results/{sample}/shovill/contigs.fa"
    output:
        abricate="results/{sample}/abricate_vfdb.tsv"
    shell:
        "abricate --db vfdb {input.assembly} > {output.abricate}"

rule mlst:
    input:
        assembly="results/{sample}/shovill/contigs.fa"
    output:
        mlst="results/{sample}/mlst.tsv"
    shell:
        "mlst {input.assembly} > {output.mlst}"

rule prokka:
    input:
        assembly="results/{sample}/shovill/contigs.fa"
    output:
        prokka="results/{sample}/prokka/{sample}.gff"
    conda:
        "prokka_env.yml"
    shell:
        """
            conda run -n prokka_env prokka --centre X \
            --compliant {input.assembly} \
            --outdir results/{wildcards.sample}/prokka/ \
            --force --cpus 4 \
            --prefix {wildcards.sample}
        """

checkpoint build_dataframe:
    input:
        abricate_ncbi=expand("results/{sample}/abricate_ncbi.tsv", sample=sample_names),
        abricate_plasmidfinder=expand("results/{sample}/abricate_plasmidfinder.tsv", sample=sample_names),
        abricate_vfdb=expand("results/{sample}/abricate_vfdb.tsv", sample=sample_names),
        mlst=expand("results/{sample}/mlst.tsv", sample=sample_names),
        shovill=expand("results/{sample}/shovill/contigs.fa", sample=sample_names)
    output:
        dataframe="dataframe/results.csv"
    shell:
        "python scripts/build_dataframe.py"

# Helper function to get the species based on the sample
def get_species_list_for_roary():
    df = pd.read_csv("dataframe/results.csv")

    # Group by species and samples
    df = df.groupby(["SPECIES", "SAMPLE"]).size().reset_index()
    df = df[['SPECIES', 'SAMPLE']]

    # Group by species and count the number of samples
    df_grouped = df.groupby("SPECIES").count().reset_index()

    # Get the species that have more than one sample
    eligeble_species = df_grouped.loc[df_grouped["SAMPLE"] > 1, "SPECIES"].tolist()

    df.loc[df["SPECIES"].isin(eligeble_species), "SAMPLE"].tolist()

    # Return only the species that have more than one sample
    return df.loc[df["SPECIES"].isin(eligeble_species), "SPECIES"].tolist()

def get_species_list_for_roary_unique():
    df = pd.read_csv("dataframe/results.csv")

    # Group by species and samples
    df = df.groupby(["SPECIES", "SAMPLE"]).size().reset_index()
    df = df[['SPECIES', 'SAMPLE']]

    df_grouped = df.groupby("SPECIES").count().reset_index()

    eligeble_species = df_grouped.loc[df_grouped["SAMPLE"] > 1, "SPECIES"].tolist()

    return list(set(df.loc[df["SPECIES"].isin(eligeble_species), "SPECIES"].tolist()))

rule pangenome:
    input:
        dataframe="dataframe/results.csv"
    output:
        "flags/.pangenome"
    run:
        # Create folders in output directory for each species based on the dataframe
        df = pd.read_csv(input.dataframe)

        # Get eligeble species for pangenome analysis
        eligeble_species = get_species_list_for_roary_unique()

        if not os.path.exists(f"pangenome"):
            os.system(f"mkdir pangenome")

        # For every species get the samples associated with that species and copy theri
        # gffs from prokka in their corresponding folder
        for species in eligeble_species:

            # Get the samples associated with the species
            samples = df.loc[df["SPECIES"] == species, "SAMPLE"].tolist()

            if not os.path.exists(f"pangenome/{species}"):
                os.system(f"mkdir pangenome/{species}")

            # If the folder already exists delete it
            if os.path.exists(f"pangenome/{species}/output"):
                os.system(f"rm -rf pangenome/{species}/output")

            # Copy the samples to the corresponding folder
            for sample in samples:
                if f"{sample}.gff" not in os.listdir(f"pangenome/{species}"):
                    os.system(
                        f"cp results/{sample}/prokka/{sample}.gff pangenome/{species}/{sample}.gff"
                    )
        
        # Run roary for each species
        for species in eligeble_species:
            os.system(f"roary pangenome/{species}/*.gff -f pangenome/{species}/output -e -n -p 16")

        os.system("touch flags/.pangenome")

rule build_report:
    input:
        pangenome="flags/.pangenome",
    output:
        Box_Contig_Length="report/QC_Box_Plot.html",
        DF_Full_Table="report/Gene_Table.html",
        DF_Reads_Table="report/QC_Table.html",
        Heatmap_Pangenomic="report/Pangenomic_Heatmap.html",
        Heatmap_Plasmids_Full_Figure_Coverage="report/Plasmid_Gene_Heatmap.html",
        Heatmap_Resistance_Full_Figure="report/Resistance_Heatmap.html",
        Heatmap_Virulence_Full_Figure_Coverage="report/Virulence_Heatmap.html",
        Network_Samples_Figure="report/Network_Chart.html",
        Pangenome_Pie_Chart="report/Pangenome_Pie_Chart.html",
        Scatter_Contig_Length="report/QC_Scatter_Plot.html",
        Subtype_HTML_String="report/Subtype_Pie_Charts.html",
        Sunburst_Figure="report/Sunburst_Chart.html",
        flag="flags/.report"
    shell:
        "python scripts/build_report.py"
