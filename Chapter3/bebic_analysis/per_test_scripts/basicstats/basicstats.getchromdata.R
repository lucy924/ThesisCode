#!/usr/bin/env Rscript

# Load necessary library for argument parsing (optional)
library(optparse)
# library(bedtorch)
library("data.table")

# Define command line options
option_list <- list(
    make_option(c("--sample_name"),
        type = "character", default = NULL, help = "sample name eg B3_Control"
    ),
    make_option(c("--chrom"),
        type = "character", default = NULL, help = "Chromosome eg chr1"
    ),
    make_option(c("--modBase"),
        type = "character", default = NULL, help = "base modification looking at, ie A, C"
    ),
    make_option(c("--motif"),
        type = "character", default = NULL, help = "motif, ie A, C, CpG, AGG"
    ),
    make_option(c("--name_of_test"),
        type = "character", default = NULL, help = "test1"
    ),
    make_option(c("--path2pileups"),
        type = "character", default = NULL,
        help = "path to pileup files to use, e.g. /external/analyses/lucy/bebic_pipeline/test2/pileups"
    ),
    make_option(c("--path2_intermediate_data"),
        type = "character", default = NULL, help = "path to intermediate data fp"
    ),
    make_option(c("--cov_threshold"),
        type = "integer", default = 20, help = "coverage threshold to keep. default 20"
    )
)

# Parse the command line arguments
opt_parser <- OptionParser(option_list = option_list)
options <- parse_args(opt_parser)

sample_name <- options$sample_name
chrom <- options$chrom
modBase <- options$modBase
motif <- options$motif
name_of_test <- options$name_of_test
path2pileups <- options$path2pileups
path2_intermediate_data <- options$path2_intermediate_data
cov_threshold <- options$cov_threshold
# print("cov_threshold")
# print(cov_threshold)

# save file name
save_chrom_info_fp <- paste0(
    path2_intermediate_data, "/",
    sample_name, ".", chrom, ".", name_of_test, ".meth_data.Rdata"
)

# for (chrom in chr_names) {
log_level <- "    "
fp <- paste0(
    path2pileups, "/",
    sample_name, ".", chrom, ".", name_of_test, ".pileup.bed.gz"
    # e.g. /external/analyses/lucy/bebic_pipeline/test2/pileups/B3_Control.chr1.test2.pileup.bed.gz
    # sample_name, ".", chrom, ".open_pileup.betterCmods.bed.gz" # testing
)

# Check if file exists
if (!file.exists(fp)) {
    print(paste0(log_level, fp, " file does not exist, moving on"))
    flush.console() # Ensure the output is printed immediately
    quit(status = 1)
}

#  Read in modBase counts
modBase_df <- fread(
    file = fp, head = FALSE
)
columns = c("chrom", "start", "end", "name", "score", "strand", "thickStart", "thickEnd", "color", "valid_coverage", "percent_modified", "count_modified", "count_canonical", "count_other_mod", "count_delete", "count_fail", "count_diff", "count_nocall")
setnames(modBase_df, columns) # when run pileup without header

# remove bases below threshold
modBase_df <- modBase_df[modBase_df$valid_coverage >= cov_threshold, ]

# remove modified_base_code a from name column, it's a wrong call using the cpg
# ensure everything is a capital
modBase_df[, name := toupper(name)]
if (modBase == "A") {
    modBase_df <- subset(modBase_df, name == "A")
} else {
    modBase_df <- subset(modBase_df, name != "A")
}

# Get meth, unmeth, and beta_val
meth <- modBase_df$count_modified
unmeth <- modBase_df$count_canonical

# Get percent_mods
modBase_percent_modified <- modBase_df$percent_modified

# get associated location info
feature_locs <- paste0(modBase_df$chrom, ":", modBase_df$start, "-", modBase_df$end, ":", modBase_df$strand)

# save data
save(
    meth, unmeth, modBase_percent_modified,
    feature_locs,
    file = save_chrom_info_fp
)
