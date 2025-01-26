suppressMessages(library("GenomicRanges"))
suppressMessages(library("tidyr"))
suppressMessages(library("data.table"))
suppressMessages(library("dplyr"))
suppressMessages(library("minfi"))

################# Pick variables #################
name_of_test <- "test11"
bebic <- "BEBIC_B"
# bebic = "BEBIC_C"
exp_list = c("B3", "B4", "B5", "B6")
# exp_list = c("C3", "C5", "C6")
# modBase = "C"
modBase <- "A"
motif <- "A"

###################################################

################# Pipeline paths #################

wd <- getwd()
print(paste0("working directory: ", wd))
ad <- "/external/analyses/lucy/bebic_pipeline"
path2_intermediate_data <- paste0(ad, "/", name_of_test, "/intermediate_data_files")
output_dir <- path2_intermediate_data

print(paste0("output directory for intermediate files: ", output_dir))
if (!dir.exists(output_dir)) {
    dir.create(output_dir)
}

# output_csv <- paste0(output_dir, "/", bebic, ".", motif, "minfi_dmp_results.csv")
output_csv <- paste0(output_dir, "/sig_DMPs.", bebic, ".minfi.csv")

###################################################

################# Set up stuff #################
sample_names_list <- c()
for (exp in exp_list) {
    for (var in c("Control", "Test")) {
        sample_name = (paste0(exp, "_", var))
        sample_names_list = c(sample_names_list, sample_name)
    }
}

chr_names = c(
    "chr1", "chr2", "chr3", "chr4", "chr5",
    "chr6", "chr7", "chr8", "chr9", "chr10",
    "chr11", "chr12", "chr13", "chr14", "chr15",
    "chr16", "chr17", "chr18", "chr19", "chr20",
    "chr21", "chr22", "chrX", "chrY", "chrM"
)
colour_mapping <- c(
    "B3_Control" = "#b9647a",
    "B3_Test" = "#E3C1CA",
    "B4_Control" = "#c843c9",
    "B4_Test" = "#E9B5E9",
    "B5_Control" = "#d34143",
    "B5_Test" = "#EDB3B4",
    "B6_Control" = "#c059a3",
    "B6_Test" = "#E6BCDA",
    "C3_Control" = "#3b6ae1",
    "C3_Test" = "#BEE4E1",
    "C5_Control" = "#5cbab4",
    "C5_Test" = "#B1C4F3",
    "C6_Control" = "#6c8db3",
    "C6_Test" = "#C4D1E0"
)

###################################################

# Load in data

filename <- paste0(output_dir, "/", bebic, ".", motif, ".minfi_dmp.Rdata")
load(file = filename)
dmp_results <- dmp_results_with_genomic_pval
remove(dmp_results_with_genomic_pval)
gc() # garbage collection

# -------------------- #

# dataframe I'm aiming for:
# 	chrom	start	end	score	base	a_methylated	a_hydroxymethylated	a_total	b_methylated	b_hydroxymethylated	b_total	a_pct_modified	b_pct_modified	map_pvalue	effect_size	balanced_map_pvalue	balanced_effect_size	num_a_samples	num_b_samples	name
# 0	chr1	121695871	121695872	4.099420	C	15	1	141	34	0	146	0.113475	0.232877	0.327975	-0.119402	0.015769	-0.151005	3	4	1.1
# 1	chr1	124153984	124153985	4.901008	A	3	0	88	20	0	119	0.034091	0.168067	0.002442	-0.133976	0.005758	-0.119363	3	3	1.2
# 2	chr1	124177350	124177351	2.108251	C	42	2	140	28	1	149	0.314286	0.194631	0.033023	0.119655	0.034888	0.109184	4	3	1.3
# 3	chr1	124203828	124203829	2.788887	C	25	0	102	17	0	146	0.245098	0.116438	0.084494	0.128660	0.001078	0.177288	3	3	1.4
# 4	chr1	124204400	124204401	2.835215	C	49	1	105	104	1	165	0.476190	0.636364	0.041863	-0.160173	0.006696	-0.191870	3	4	1.5

# Important columns
# chrom	start	score	base    balanced_effect_size    name
# Note that I don't add the name until the dataframe has been merged C and A together.
# Which is good cos I did that code in python already

# -------------------- #
# first add base column
dmp_results$base <- modBase

# -------------------- #
# next remove columns end, strand, qval
dmp_results$end <- NULL
dmp_results$strand <- NULL
dmp_results$qval <- NULL
# next rename chromosome, intercept and f columns
colnames(dmp_results)[colnames(dmp_results) == "chromosome"] <- "chrom"
colnames(dmp_results)[colnames(dmp_results) == "intercept"] <- "balanced_effect_size"
colnames(dmp_results)[colnames(dmp_results) == "f"] <- "score"

dmp_results <- dmp_results %>% filter(pval <= 0.01)
# dmp_results_allBase <- rbind(dmp_results_Cmod, dmp_results_Amod)
dmp_results <- dmp_results[order(dmp_results$pval), ]

# Write to csv
write.csv(dmp_results, file = output_csv, row.names = FALSE)
