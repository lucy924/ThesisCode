#!/bin/bash -e

# Activate the environment
source /home/dejlu879/miniforge3/etc/profile.d/conda.sh
source /home/dejlu879/miniforge3/etc/profile.d/mamba.sh
export PYTHONNOUSERSITE=1 # don't add python user site library to path
mamba activate gen_venv_mamba

# analyses directory
ad="/external/analyses/20241104-BEBIC_Karyotyping"
cd $ad

path2alignedbam="/external/analyses/20240416-BEBIC_B/BT0/demuxed/B3B5_T0.bam"
path2unalignedbam="/external/analyses/20240416-BEBIC_B/BT0/unaligned_bams/B3B5_T0.unaligned.bam"
path2ref="/home/dejlu879/refs/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"

if [ -f "$path2unalignedbam" ]; then
    echo "$path2unalignedbam exists."
    # Perform your action here
else
    echo "$path2unalignedbam does not exist. Performing unalignment."
    samtools reset -o $path2unalignedbam $path2alignedbam
fi

nextflow run epi2me-labs/wf-human-variation \
	--bam $path2unalignedbam \
	--ref $path2ref \
	--sample_name 'B3B5_T0' \
	--cnv \
	--sv \
	--phased \
    --override_basecaller_cfg dna_r10.4.1_e8.2_400bps_sup@v4.3.0 \
    --out_dir $ad \
    --bam_min_coverage 0.2 \
    -resume \
	-profile standard
