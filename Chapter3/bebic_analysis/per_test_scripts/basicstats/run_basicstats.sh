#!/bin/bash -e

name_of_test="test11"
modBase="A"
motif="A"
bebic="BEBIC_B"

# Sample names
sample_names=(
    "B3_Control"
    "B3_Test"
    "B4_Control"
    "B4_Test"
    "B5_Control"
    "B5_Test"
    "B6_Control"
    "B6_Test"
)

# working directory
wd=$(pwd)

ad="/external/analyses/lucy/bebic_pipeline"  # analysis dir
path2_intermediate_data="${ad}/${name_of_test}/intermediate_data_files"
path2pileups="${ad}/${name_of_test}/pileups"
path2_plots="${wd}/plots"
# path2pileups="${ad}/${name_of_test}/pileups"

mkdir -p "${path2_intermediate_data}"
mkdir -p "${path2_plots}"

# Activate the environment
source /opt/miniforge3/etc/profile.d/conda.sh
source /opt/miniforge3/etc/profile.d/mamba.sh
export PYTHONNOUSERSITE=1 # don't add python user site library to path
mamba activate gen_venv

# Chromosomes to process
chromosomes=(
    "chr1"
    "chr2"
    "chr3"
    "chr4"
    "chr5"
    "chr6"
    "chr7"
    "chr8"
    "chr9"
    "chr10"
    "chr11"
    "chr12"
    "chr13"
    "chr14"
    "chr15"
    "chr16"
    "chr17"
    "chr18"
    "chr19"
    "chr20"
    "chr21"
    "chr22"
    "chrX"
    "chrY"
    "chrM"
    # "chrEBV"
    # "chrL"
    # "NC_008769.1"
    # "unmapped"
)

for sample_name in "${sample_names[@]}"; do
    echo "Running sample ${sample_name}" 
    for chrom in "${chromosomes[@]}"; do
        echo "  Running chromosome ${chrom}" 
        Rscript scripts/basicstats/basicstats.getchromdata.R --sample_name "${sample_name}" --chrom "${chrom}" --modBase "${modBase}" --motif "${motif}" --name_of_test "${name_of_test}" --path2pileups "${path2pileups}" --path2_intermediate_data "${path2_intermediate_data}"
    done
done

echo ""
echo "Running plot_all_sample_data.R" 
Rscript scripts/basicstats/basicstats.plot_all_sample_data.R --bebic "${bebic}" --name_of_test "${name_of_test}" --modBase "${modBase}" --motif "${motif}" --path2_intermediate_data "${path2_intermediate_data}" --bw_val 1.2 

echo "Running plot_MDS.R" 
Rscript scripts/basicstats/basic_stats.plot_MDS.R --bebic "${bebic}" --name_of_test "${name_of_test}" --modBase "${modBase}" --motif "${motif}" --path2_intermediate_data "${path2_intermediate_data}" --path2_plots "${path2_plots}"
