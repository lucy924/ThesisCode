suppressMessages(library("GenomicRanges"))
suppressMessages(library("tidyr"))
suppressMessages(library("data.table"))
suppressMessages(library("dplyr"))
suppressMessages(library("minfi"))

################# Pick variables #################
bebic <- "BEBIC_B"
# bebic <- "BEBIC_C"
exp_list <- c("B3", "B4", "B5", "B6")
# exp_list <- c("C3", "C5", "C6")
# modBase = "C"
modBase <- "A"
motif <- "A"
name_of_test <- "test11"

###################################################

################# Pipeline paths #################

wd <- getwd()
print(paste0("working directory: ", wd))

ad <- "/external/analyses/lucy/bebic_pipeline"
path2_intermediate_data <- paste0(ad, "/", name_of_test, "/intermediate_data_files")
path2_plots <- paste0(getwd(), "/plots/")
output_dir <- paste0(path2_intermediate_data, "/")

print(paste0("output directory for intermediate files: ", output_dir))
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

print(paste0("output directory for plots: ", path2_plots))
if (!dir.exists(path2_plots)) {
    dir.create(path2_plots)
}

###################################################

chr_names <- c(
    "chr1", "chr2", "chr3", "chr4", "chr5",
    "chr6", "chr7", "chr8", "chr9", "chr10",
    "chr11", "chr12", "chr13", "chr14", "chr15",
    "chr16", "chr17", "chr18", "chr19", "chr20",
    "chr21", "chr22", "chrX", "chrY", "chrM"
)

sample_names_list <- c()
for (exp in exp_list) {
    for (var in c("Control", "Test")) {
        sample_name <- (paste0(exp, "_", var))
        sample_names_list <- c(sample_names_list, sample_name)
    }
}

###################################################

filename_mset <- paste0(
    path2_intermediate_data, "/",
    name_of_test, ".gMset.Rdata"
)
# stuff in this file: gMset, new_colnames
load(file = filename_mset)

###################################################

beta_matrix <- getBeta(gMset)
colnames(beta_matrix) <- new_colnames # add colnames back

sample_names <- colnames(beta_matrix) # Get sample names from your meth matrix
if (bebic == "BEBIC_B") {
    pheno <- gsub("B[0-9]+_", "", sample_names) # Strip away the prefix to get "Control", "Test" in correct order
} else if (bebic == "BEBIC_C") {
    pheno <- gsub("C[0-9]+_", "", sample_names) # Strip away the prefix to get "Control", "Test" in correct order
}

# Do DMP
dmp_results <- dmpFinder(
    dat = beta_matrix,
    pheno = pheno,
    type = "categorical",
    #   qCutoff = 1,
    shrinkVar = TRUE # variance shrinkage for when small num of samples. will generate a warning about contrasts
)

# Extract genomic information based on the indices in dmp_results
indices <- as.numeric(rownames(dmp_results)) # Get the row indices

# Extract genomic locations
# Create the GRanges object
# sample_gr <- GRanges(
#     seqnames = seqnames(Mset),
#     ranges = IRanges(start = start(Mset), end = end(Mset)),
#     strand = strand(Mset)
# )

# Step 2: Map Indices to Genomic Locations
genomic_info <- gMset[indices]
genomic_locations <- data.frame(
    chromosome = seqnames(genomic_info),
    start = start(genomic_info),
    end = end(genomic_info),
    strand = strand(genomic_info)
)

# Step 3: Combine with DMP Results
dmp_results_with_genomic <- cbind(genomic_locations, dmp_results)
# cbind = column bind

# Filter to p-val 0.05
dmp_results_with_genomic_pval <- dmp_results_with_genomic[dmp_results_with_genomic$pval < 0.05, ]
# Filter to effect size >= 10%
dmp_results_with_genomic_pval <- dmp_results_with_genomic_pval[(dmp_results_with_genomic_pval$intercept <= -0.1 | dmp_results_with_genomic_pval$intercept >= 0.1), ]

print("saving result...")
###################################
# output to rdata
save(dmp_results_with_genomic_pval, file = paste0(path2_intermediate_data, "/", bebic, ".", motif, ".minfi_dmp.Rdata"))

print("done.")
