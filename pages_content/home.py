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
    st.title("Biosynthetic Gene Clusters in the vaginal microbiome")
    
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
    #SBATCH --output=path_to_antismash/01_run_antiSMASH/logs/%x_%A_%a.out
    #SBATCH --error=path_to_antismash/01_run_antiSMASH/logs/%x_%A_%a.err
    #SBATCH --array=1-1500
    #SBATCH --mem-per-cpu=40G
    #SBATCH --cpus-per-task=4

    # Conda activation
    source /usr/local/packages/miniconda3/etc/profile.d/conda.sh
    conda activate /usr/local/packages/miniconda3/envs/antismash

    # Go to working directory
    cd path_to_antismash/01_run_antiSMASH/

    infile=$(sed -n "${SLURM_ARRAY_TASK_ID}p" MAGs_list/mag_0001_1500.txt)

    # Run AntiSMASH
    /usr/local/packages/miniconda3/envs/antismash/bin/antismash \\
        --cpus 4 \\
        path_to_MAGS/${infile}.fasta \\
        --output-dir path_to_antismash/results_antismash/${infile}/ \\
        --genefinding-tool prodigal \\
        --cb-general --cb-subclusters --cb-knownclusters --fullhmmer \\
        --clusterhmmer --tigrfam --asf --cc-mibig --pfam2go --rre --smcog-trees
    """, language="bash")

    st.markdown("""
    <div class="justified-text">
                
    3. **Post-processing**: AntiSMASH outputs can be processed with custom scripts to extract JSON data for further analysis.
    </div>
    """, unsafe_allow_html=True)


    st.subheader("antiSMASH documentation", divider='grey')
    st.markdown("<a href='https://docs.antismash.secondarymetabolites.org/' target='_blank'>antiSMASH documentation website</a>",unsafe_allow_html=True)
    st.markdown("<a href='https://docs.antismash.secondarymetabolites.org/glossary/' target='_blank'>antiSMASH Glossary</a>",unsafe_allow_html=True)
    st.markdown("<a href='https://antismash-db.secondarymetabolites.org/' target='_blank'>antiSMASH database</a>",unsafe_allow_html=True)
