#!/bin/bash -e

# Activate the environment
source /home/dejlu879/miniforge3/etc/profile.d/conda.sh
source /home/dejlu879/miniforge3/etc/profile.d/mamba.sh
export PYTHONNOUSERSITE=1 # don't add python user site library to path
mamba activate bebic_venv
# modkit version 0.4.1

# Options to pick:
name_of_test="test11"
# mod_type="trad"

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
    # "C3_Control"
    # "C3_Test"
    # "C5_Control"
    # "C5_Test"
    # "C6_Control"
    # "C6_Test"
    # "B3B5_Baseline"
    # "B4B6_Baseline"
)

ad="/external/analyses/lucy/bebic_pipeline"  # analysis dir

# path2_intermediate_data <- paste0(ad, "/", name_of_test, "/intermediate_data_files")
# pileup directory
pileup_dir="${ad}/${name_of_test}/pileups"
mkdir -p $pileup_dir

files_to_skip_dmr="${pileup_dir}/non-existent_pileups.txt"
touch $files_to_skip_dmr

# References per chromosome
# path2refs="/home/dejlu879/refs/T2T_chromsonly_ebv_lambda_bcg"
reference="/home/dejlu879/refs/T2T/T2T_chromsonly.ebv.lambda.bcg.fa"
bed_cytosines="/home/dejlu879/refs/T2T/T2T_chromsonly_ebv_lambda_bcg.noAGG.bed"

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

echo "============== Beginning pileup =============="

for chr in "${chromosomes[@]}"; do
    # reference="${path2refs}/${chr}.fa"
    if [ ! -e "${reference}.fai" ]; then
        samtools faidx "${reference}"
    fi

    for sample in "${sample_names[@]}"; do
        echo "-----| ${chr}, ${sample} open pileup started ... |-----"

        pileup_bed="${pileup_dir}/${sample}.${chr}.${name_of_test}.pileup.bed"
        log_file="${pileup_dir}/${sample}.${chr}.${name_of_test}.pileup.log"
        
        if [ -e "${pileup_bed}.gz" ]; then
            echo "${pileup_bed}.gz exists, skipping this one..."
            continue
        fi

        if [[ "${sample}" == B* ]]; then
            dir_name="20240416-BEBIC_B"
        elif [[ "${sample}" == C* ]]; then
            dir_name="20240628-BEBIC_C"
        fi

        path2reads_bam="/external/analyses/${dir_name}/t2t_bams_by_sample/${sample}/split_by_chr/${sample}.t2t.chr_ebv_l_bcg.${chr}.bam"

        if [ ! -e ${path2reads_bam} ]; then
            echo "${path2reads_bam} does not exist, adding to file and skipping."
            echo "${path2reads_bam}" >> $files_to_skip_dmr
            continue
        fi
        
        modkit pileup ${path2reads_bam} ${pileup_bed} \
        --log-filepath $log_file \
        --threads 25 \
        --motif A 0 \
        --combine-mods \
        --include-bed ${bed_cytosines} \
        --ref "${reference}"

        # had to remove header so that indexing will work
        # --with-header \

        echo "${chr}, ${sample} pileup done"

        bgzip ${pileup_bed}
        tabix -p bed ${pileup_bed}.gz
        echo "------------------------------"
        sleep 5

    done
    sleep 1
done

echo "========== All completed yay! =========="

# modkit version: mod_kit 0.4.1
