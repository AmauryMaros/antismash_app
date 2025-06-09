import streamlit as st

# Apply custom CSS for text justification
st.markdown("""
    <style>
    .justified-text {
        text-align: justify;
    }
    </style>
    """, unsafe_allow_html=True)

def page():
    st.title("AntiSMASH (Antibiotics & Secondary Metabolite Analysis Shell)")
    
    st.markdown("""
    <div class="justified-text">
    AntiSMASH is a bioinformatics tool used to identify and analyze biosynthetic gene clusters (BGCs) in bacterial and fungal genomes. 
                
    These clusters are groups of genes that work together to produce secondary metabolites, which are compounds not directly involved in growth or reproduction but often serve important ecological functions, such as antibiotics or pigments.


    </div>
    """, unsafe_allow_html=True)

    st.subheader("Workflow", divider='grey')

    st.markdown("""
    <div class="justified-text">
    
    1. **Input**: Typically, the input is a genome sequence in FASTA format or a GenBank file.
    
    2. **Running AntiSMASH**: AntiSMASH can be run locally or via antiSMASH web version.
    
    </div>
    """, unsafe_allow_html=True)

    # Display the SLURM script code block
    st.code("""
    #!/bin/bash

    #SBATCH --job-name=antismash
    #SBATCH --output=/local/scratch/amaros/antismash/01_run_antiSMASH/logs/%x_%A_%a.out
    #SBATCH --error=/local/scratch/amaros/antismash/01_run_antiSMASH/logs/%x_%A_%a.err
    #SBATCH --array=1-1500
    #SBATCH --mem-per-cpu=40G
    #SBATCH --cpus-per-task=4

    # Conda activation
    source /usr/local/packages/miniconda3/etc/profile.d/conda.sh
    conda activate /usr/local/packages/miniconda3/envs/antismash

    # Go to working directory
    cd /local/scratch/amaros/antismash/01_run_antiSMASH/

    infile=$(sed -n "${SLURM_ARRAY_TASK_ID}p" MAGs_list/mag_0001_1500.txt)

    # Run AntiSMASH
    /usr/local/packages/miniconda3/envs/antismash/bin/antismash \\
        --cpus 4 \\
        /local/projects-t3/LSVF/VIRGO2/final_bins/${infile}.fasta \\
        --output-dir /local/scratch/amaros/antismash/results_antismash/${infile}/ \\
        --genefinding-tool prodigal \\
        --cb-general --cb-subclusters --cb-knownclusters --fullhmmer \\
        --clusterhmmer --tigrfam --asf --cc-mibig --pfam2go --rre --smcog-trees
    """, language="bash")

    st.markdown("""
    <div class="justified-text">
                
    3. **Post-processing**: AntiSMASH outputs can be processed with custom scripts to extract JSON data for further analysis.
    </div>
    """, unsafe_allow_html=True)

    st.code("""
[amaros@ravellab 02_pipeline]$ ll
total 248
-rwxr-xr-x. 1 amaros igs  2732 Nov 25 23:11 01_process_antismash_output.py
-rwxr-xr-x. 1 amaros igs  1602 Nov 25 23:11 02_process_antismash_json.py
-rwxr-xr-x. 1 amaros igs  2162 Nov 26 11:28 03_get_ctg_coordinates.py
-rwxr-xr-x. 1 amaros igs  2108 Nov 26 11:50 04_get_ctg_sequences.py
-rwxr-xr-x. 1 amaros igs   875 Nov 26 11:43 05_blast_ctg.py
-rwxr-xr-x. 1 amaros igs   304 Nov 26 11:43 blast_ctg.sh
-rwxr-xr-x. 1 amaros igs 19800 Nov 25 23:11 functions.py
-rwxr-xr-x. 1 amaros igs   808 Nov 25 23:11 get_json_status.py
-rwxr-xr-x. 1 amaros igs   597 Nov 25 23:11 parse_jsons.sh
drwxr-xr-x. 2 amaros igs    85 Nov 25 23:11 __pycache__
    """, language="bash")


    st.subheader("antiSMASH documentation", divider='grey')
    st.markdown("<a href='https://docs.antismash.secondarymetabolites.org/' target='_blank'>antiSMASH documentation website</a>",unsafe_allow_html=True)
    st.markdown("<a href='https://docs.antismash.secondarymetabolites.org/glossary/' target='_blank'>antiSMASH Glossary</a>",unsafe_allow_html=True)
    st.markdown("<a href='https://antismash-db.secondarymetabolites.org/' target='_blank'>antiSMASH database</a>",unsafe_allow_html=True)






# # Run the page function
# page()


    # #### Key Features and Capabilities

    # 1. **Biosynthetic Gene Cluster Detection**: AntiSMASH detects BGCs by identifying and annotating genes known to be associated with secondary metabolite synthesis pathways. This includes pathways for polyketides, nonribosomal peptides, terpenes, and more.

    # 2. **Annotation and Classification**: After detecting BGCs, AntiSMASH classifies them based on type and annotates functional genes, providing users with insights into potential chemical products and biosynthetic capabilities.

    # 3. **Comparative Analysis**: AntiSMASH can compare BGCs across multiple genomes, helping researchers understand evolutionary relationships and predict possible new compounds.

    # 4. **Genome Mining**: AntiSMASH integrates with databases like MIBiG (Minimum Information about a Biosynthetic Gene Cluster) to provide data on known BGCs, allowing researchers to identify similar or novel clusters that may produce uncharacterized compounds.

    # 5. **Output**: AntiSMASH provides both JSON and visual outputs, allowing users to view cluster locations, gene annotations, and predicted products in a user-friendly format. The output can be processed for downstream analyses or visualized as gene cluster diagrams.