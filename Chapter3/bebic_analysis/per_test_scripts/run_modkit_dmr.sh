#!/bin/bash -e

# mamba activate bebic_venv
# modkit version 0.4.1

############## Set up environment ########################################

source /opt/miniforge3/etc/profile.d/conda.sh
source /opt/miniforge3/etc/profile.d/mamba.sh
export PYTHONNOUSERSITE=1 # don't add python user site library to path
mamba activate gen_venv

############# USER CHANGES AND OPTIONS ###################################
name_of_test="test11"
# description of test: test7 is the AGG motif test for A modified positions in motif AGG, read depth threshold of 20 and effect size of 10%, with Control Vs Test for BEBIC_C.

bebic="BEBIC_B"
# bebic="BEBIC_C"

# Experiment names
exp_names=(
    "B3"
    "B4"
    "B5"
    "B6"
)
# exp_names=(
#     "C3"
#     "C5"
#     "C6"
# )

threads=12
read_depth_threshold=20
min_effect_size=0.1  # 0.1 = 10% difference

# working directory
wd=$(pwd)
# $(...): This syntax captures the output of the command inside the parentheses and assigns it to the variable current_directory.

# analyses directory
ad="/external/analyses/lucy/bebic_pipeline"

# References per chromosome
path2refs="/home/dejlu879/refs/T2T/" # all aligned to T2T_chromsonly.ebv.lambda.bcg.fa"

# Pileup paths
path2pileups="${ad}/${name_of_test}/pileups"
files_to_skip_dmr="${path2pileups}/non-existent_pileups.txt"
touch $files_to_skip_dmr

path2results="${ad}/${name_of_test}/modkit"
mkdir -p "${path2results}"

################# Stays the same #######################################

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
)

# ----------------- Checking for non-existing files ----------------- #

# Function to read the list of nonexistent files into an array
update_nonexistent_files() {
  mapfile -t nonexistent_files < "${files_to_skip_dmr}"
}

# Function to check if a file is in the list of nonexistent files
is_nonexistent() {
  local file=$1
  for non_file in "${nonexistent_files[@]}"; do
    if [[ "${file}" == "${non_file}" ]]; then
      return 0  # File is in the list of nonexistent files
    fi
  done
  return 1  # File is not in the list
}

# ----------------- Update exp lists ----------------- #
# Update the list of broken experiment-chromosome pairs
update_nonexistent_files

################# Start Analysis #######################################

echo "============== Beginning DMR on ${bebic} =============="
# mod_type="CpG"
path2all_test_results="${path2results}/all_test"

for chr in "${chromosomes[@]}"; do
    reference="${path2refs}/T2T_chromsonly.ebv.lambda.bcg.fa"  # this was ref used for pileup
    echo "-------------- DMR ${chr} starting --------------"

    if [[ "${bebic}" == "BEBIC_B" ]]; then
      control_bed_B3="${path2pileups}/B3_Control.${chr}.${name_of_test}.pileup.bed.gz"
      test_bed_B3="${path2pileups}/B3_Test.${chr}.${name_of_test}.pileup.bed.gz"
      control_bed_B4="${path2pileups}/B4_Control.${chr}.${name_of_test}.pileup.bed.gz"
      test_bed_B4="${path2pileups}/B4_Test.${chr}.${name_of_test}.pileup.bed.gz"
      control_bed_B5="${path2pileups}/B5_Control.${chr}.${name_of_test}.pileup.bed.gz"
      test_bed_B5="${path2pileups}/B5_Test.${chr}.${name_of_test}.pileup.bed.gz"
      control_bed_B6="${path2pileups}/B6_Control.${chr}.${name_of_test}.pileup.bed.gz"
      test_bed_B6="${path2pileups}/B6_Test.${chr}.${name_of_test}.pileup.bed.gz"
    elif [[ "${bebic}" == "BEBIC_C" ]]; then
      control_bed_C3="${path2pileups}/C3_Control.${chr}.${name_of_test}.pileup.bed.gz"
      test_bed_C3="${path2pileups}/C3_Test.${chr}.${name_of_test}.pileup.bed.gz"
      control_bed_C5="${path2pileups}/C5_Control.${chr}.${name_of_test}.pileup.bed.gz"
      test_bed_C5="${path2pileups}/C5_Test.${chr}.${name_of_test}.pileup.bed.gz"
      control_bed_C6="${path2pileups}/C6_Control.${chr}.${name_of_test}.pileup.bed.gz"
      test_bed_C6="${path2pileups}/C6_Test.${chr}.${name_of_test}.pileup.bed.gz"
    fi

    dmr_result="${path2all_test_results}/raw_DMPs.${bebic}.${chr}.bed"
    dmr_log="${dmr_result}.log"
    dmr_segment="${path2all_test_results}/raw_DMRs.${bebic}.${chr}.bed"

    if [ -e "$dmr_result" ]; then
        echo "Analysis for this done already."
        continue
    fi

    if [ "${chr}" == "chrEBV" ]; then
        echo "No valid pairs for ${bebic} on chrEBV."
        continue
    fi

    # Check each chromosome pileup
    update_nonexistent_files

    for exp in "${exp_names[@]}"; do
      echo "Checking Exp ${exp}"

      # Check if either file is in the list of nonexistent files
      if [[ "${exp}" == B* ]]; then
          dir_name="20240416-BEBIC_B"
      elif [[ "${exp}" == C* ]]; then
          dir_name="20240628-BEBIC_C"
      fi

      # B3_Control.t2t.chr_ebv_l_bcg.chrEBV.bam
      control_bam="/external/analyses/${dir_name}/t2t_bams_by_sample/${exp}_Control/split_by_chr/${exp}_Control.t2t.chr_ebv_l_bcg.${chr}.bam"
      test_bam="/external/analyses/${dir_name}/t2t_bams_by_sample/${exp}_Test/split_by_chr/${exp}_Test.t2t.chr_ebv_l_bcg.${chr}.bam"

      if is_nonexistent "${control_bam}"; then
          echo "$(basename "${control_bam}") is in the list of nonexistent files"
          # continue 2  # Skip to the next experiment
      fi
      if is_nonexistent "${test_bam}"; then
          echo "$(basename "${test_bam}") is in the list of nonexistent files"
          # continue 2  # Skip to the next experiment
      fi
    done

    echo "    - performing DMR on ${chr} pileups, region discovery"
    # if [ "${chr}" == "chrY" ]; then
    #   echo "Ignore chrY, not enough datapoints"
    #   continue
    # fi

    # if [ "${chr}" == "chrM" ]; then
    #   if [ "${bebic}" == "BEBIC_B" ]; then
    #     echo "Ignore chrM, not enough datapoints??"
    #     continue
    #   fi
    # fi

    if [[ "${bebic}" == "BEBIC_B" ]]; then
      modkit dmr pair \
      -a "${control_bed_B3}" \
      -a "${control_bed_B4}" \
      -a "${control_bed_B5}" \
      -a "${control_bed_B6}" \
      -b "${test_bed_B3}" \
      -b "${test_bed_B4}" \
      -b "${test_bed_B5}" \
      -b "${test_bed_B6}" \
      -o "${dmr_result}" \
      --header \
      --ref "${reference}" \
      --segment "${dmr_segment}" \
      --base A \
      --threads "${threads}" \
      --log-filepath "${dmr_log}" \
      -f \
      --min-valid-coverage "${read_depth_threshold}" \
      --delta "${min_effect_size}"

    elif [[ "${bebic}" == "BEBIC_C" ]]; then
      modkit dmr pair \
      -a "${control_bed_C3}" \
      -a "${control_bed_C5}" \
      -a "${control_bed_C6}" \
      -b "${test_bed_C3}" \
      -b "${test_bed_C5}" \
      -b "${test_bed_C6}" \
      -o "${dmr_result}" \
      --header \
      --ref "${reference}" \
      --segment "${dmr_segment}" \
      --base A \
      --threads "${threads}" \
      --log-filepath "${dmr_log}" \
      -f \
      --min-valid-coverage "${read_depth_threshold}" \
      --delta "${min_effect_size}"
    
    fi
    # -f = force overwrite of output file if it already exists
    echo "  DMR for ${chr} completed."

done

echo "DMR completed on ${bebic}."
echo " "

# # Now do tests per experiment pair
# echo "============== Beginning DMR on ${bebic} pair by pair =============="
# # mod_type="CpG"
# path2paired_results="${path2results}/paired"

# for chr in "${chromosomes[@]}"; do
#   reference="${path2refs}/T2T_chromsonly.ebv.lambda.bcg.fa"  # this was ref used for pileup

#   for exp in "${exp_names[@]}"; do
#     echo "-------------- Open DMR ${chr}, ${exp} starting --------------"
#     control_bed="${path2pileups}/${exp}_Control.${chr}.${name_of_test}.pileup.bed.gz"
#     test_bed="${path2pileups}/${exp}_Test.${chr}.${name_of_test}.pileup.bed.gz"

#     dmr_result="${path2paired_results}/raw_DMPs.${exp}.${chr}.bed"
#     dmr_log="${dmr_result}.log"
#     dmr_segment="${path2paired_results}/raw_DMRs.${exp}.${chr}.bed"
    
#     if [ -e "$dmr_result" ]; then
#         echo "Analysis for this done already."
#         continue
#     fi

#     if [ "${chr}" == "chrY" ]; then
#         if [ "${exp}" == "B3" ]; then
#             echo "Not enough datapoints :("
#             continue
#         fi
#         if [ "${exp}" == "B4" ]; then
#             echo "Not enough datapoints :("
#             continue
#         fi
#         if [ "${exp}" == "B5" ]; then
#             echo "Not enough datapoints :("
#             continue
#         fi
#         if [ "${exp}" == "B6" ]; then
#             echo "Not enough datapoints :("
#             continue
#         fi
#     fi

#     # if [ "${chr}" == "chr14" ]; then
#     #     if [ "${exp}" == "B6" ]; then
#     #         echo "Not enough datapoints :("
#     #         continue
#     #     fi
#     # fi

#     # Check each chromosome pileup
#     update_nonexistent_files

#     # Check if either file is in the list of nonexistent files
#     if [[ "${exp}" == B* ]]; then
#         dir_name="20240416-BEBIC_B"
#     elif [[ "${exp}" == C* ]]; then
#         dir_name="20240628-BEBIC_C"
#     fi

#     # B3_Control.t2t.chr_ebv_l_bcg.chrEBV.bam
#     control_bam="/external/analyses/${dir_name}/t2t_bams_by_sample/${exp}_Control/split_by_chr/${exp}_Control.t2t.chr_ebv_l_bcg.${chr}.bam"
#     test_bam="/external/analyses/${dir_name}/t2t_bams_by_sample/${exp}_Test/split_by_chr/${exp}_Test.t2t.chr_ebv_l_bcg.${chr}.bam"

#     if is_nonexistent "${control_bam}"; then
#         echo "$(basename "${control_bam}") is in the list of nonexistent files"
#         # continue 2  # Skip to the next experiment
#     fi
#     if is_nonexistent "${test_bam}"; then
#         echo "$(basename "${test_bam}") is in the list of nonexistent files"
#         # continue 2  # Skip to the next experiment
#     fi

#     echo "    - performing DMR on open pileups, region discovery"
#     # if [ "${chr}" == "chrY" ]; then
#     #   echo "Ignore chrY, not enough datapoints"
#     #   continue
#     # fi

#     modkit dmr pair \
#     -a "${control_bed}" \
#     -b "${test_bed}" \
#     -o "${dmr_result}" \
#     --header \
#     --ref "${reference}" \
#     --segment "${dmr_segment}" \
#     --base A \
#     --threads "${threads}" \
#     --log-filepath "${dmr_log}" \
#     -f \
#     --min-valid-coverage "${read_depth_threshold}" \
#     --delta "${min_effect_size}"

#   done
#   # -f = force overwrite of output file if it already exists
#   echo "  DMR for ${chr}, ${exp} completed."

# done
# echo "DMR completed on ${bebic} pair by pair."

echo "========== All completed yay! =========="
