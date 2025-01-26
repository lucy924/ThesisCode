# #!/usr/bin/env Rscript

# Load necessary library for argument parsing (optional)
suppressMessages(library(optparse))
suppressMessages(library("data.table"))
suppressMessages(library(dplyr))
suppressMessages(library(minfi))

# Define command line options
option_list <- list(
    make_option(c("--bebic"),
        type = "character", default = NULL, help = "bebic, eg BEBIC_B"
    ),
    make_option(c("--modBase"),
        type = "character", default = NULL, help = "base modification looking at, eg A, C"
    ),
    make_option(c("--motif"),
        type = "character", default = NULL, help = "motif looking at, eg A, C, CpG, AGG"
    ),
    make_option(c("--name_of_test"),
        type = "character", default = NULL, help = "test1"
    ),
    make_option(c("--path2_intermediate_data"),
        type = "character", default = NULL, help = "path to intermediate data fp"
    ),
    make_option(c("--path2_plots"),
        type = "character", default = NULL, help = "path to plots fp"
    )
)
# Parse the command line arguments
opt_parser <- OptionParser(option_list = option_list)
options <- parse_args(opt_parser)

bebic <- options$bebic
modBase <- options$modBase
motif <- options$motif
name_of_test <- options$name_of_test
path2_intermediate_data <- options$path2_intermediate_data
path2_plots <- options$path2_plots
# bebic <- "BEBIC_B"
# modBase <- "CpG"
# name_of_test <- "test2"
# path2_intermediate_data <- "/external/analyses/lucy/bebic_pipeline/test2/intermediate_data_files"

chr_names <- c(
    "chr1", "chr2", "chr3", "chr4", "chr5",
    "chr6", "chr7", "chr8", "chr9", "chr10",
    "chr11", "chr12", "chr13", "chr14", "chr15",
    "chr16", "chr17", "chr18", "chr19", "chr20",
    "chr21", "chr22", "chrX", "chrY", "chrM"
)

# Define custom color mapping for each sample
colour_mapping <- c(
    "B3_Control" = "#b9647a", # Add more samples with colors as needed
    "B3_Test" = "#E3C1CA", # Example for additional samples
    "B4_Control" = "#c843c9", # Add more samples with colors as needed
    "B4_Test" = "#E9B5E9", # Example for additional samples
    "B5_Control" = "#d34143", # Add more samples with colors as needed
    "B5_Test" = "#EDB3B4", # Example for additional samples
    "B6_Control" = "#c059a3", # Add more samples with colors as needed
    "B6_Test" = "#E6BCDA" # Example for additional samples
)

if (bebic == "BEBIC_B") {
    exp_list <- c("B3", "B4", "B5", "B6")
} else if (bebic == "BEBIC_C") {
    exp_list <- c("C3", "C5", "C6")
}

sample_names_list <- c()
for (exp in exp_list) {
    for (var in c("Control", "Test")) {
        sample_name = (paste0(exp, "_", var))
        sample_names_list = c(sample_names_list, sample_name)
    }
}

# Load info needed for MDS
for (sample_name in sample_names_list) {
    print(paste0("MDS, loaded ", sample_name))
    
    mvals_dt_name <- paste0(
        path2_intermediate_data, "/",
        sample_name, ".", name_of_test, ".mvals_dt.Rdata"
    )

    if (file.exists(mvals_dt_name)) {
        next
    }

    # sample_WG_modBase_Mval <- list() # whole genome data
    sample_WG_modBase_meth <- list() # whole genome data
    sample_WG_modBase_unmeth <- list() # whole genome data
    sample_WG_modBase_featurelocs <- list() # whole genome data
    # sample_WG_modBase_strand <- list() # whole genome data

    for (chrom in chr_names) {
        save_chrom_info_fp <- paste0(
            path2_intermediate_data, "/",
            sample_name, ".", chrom, ".", name_of_test, ".meth_data.Rdata"
        )
        load(file = save_chrom_info_fp)
        # sample_WG_modBase_Mval[[chrom]] <- log(meth / unmeth)
        sample_WG_modBase_meth[[chrom]] <- meth
        sample_WG_modBase_unmeth[[chrom]] <- unmeth
        sample_WG_modBase_featurelocs[[chrom]] <- feature_locs
        # sample_WG_modBase_strand[[chrom]] <- modBase_strand
    }
    remove(
        meth, unmeth, modBase_percent_modified, feature_locs
    )
    suppressMessages(gc())

    # flatten out list - for RAM management
    sample_Mvals_dt <- data.table(
        feature_locs = unlist(sample_WG_modBase_featurelocs),
        # strand = unlist(sample_WG_modBase_strand),
        # Mvals = unlist(sample_WG_modBase_Mval),
        meth = unlist(sample_WG_modBase_meth),
        unmeth = unlist(sample_WG_modBase_unmeth),
        sample = rep(sample_name, length(sample_WG_modBase_featurelocs))
    )
    remove(sample_WG_modBase_featurelocs, sample_WG_modBase_meth, sample_WG_modBase_unmeth) #, sample_WG_modBase_strand)
    suppressMessages(gc())
    # sample_Mvals_dt <- sample_Mvals_dt %>%
    #     mutate(feature_locs = paste(sample_WG_modBase_featurelocs, sample_WG_modBase_strand, sep = "_"))
    # sample_Mvals_dt$strand = NULL

    # save Mvals dt
    save(sample_Mvals_dt, file = mvals_dt_name)
    remove(sample_Mvals_dt)
    suppressMessages(gc())
}

print(paste0("Construct table of all samples Mvals"))
all_samples_Mvals <- rbindlist(lapply(
    sample_names_list, function(sample_name) {
        load(paste0(
            path2_intermediate_data, "/",
            sample_name, ".", name_of_test, ".mvals_dt.Rdata"
        ))
        return(sample_Mvals_dt)
    }
))

# print(paste0("Clean up feature_locs"))
# all_samples_Mvals$feature_locs <- paste0(all_samples_Mvals$feature_locs, ":", all_samples_Mvals$strand)
# all_samples_Mvals$feature_locs <- gsub(":", "_", all_samples_Mvals$feature_locs)
# all_samples_Mvals$strand <- NULL

print(paste0("note un_meth means unmeth and meth data"))
named_list_of_un_meth_data <- list()
for (exp in exp_list) {
    cont_name <- paste0(exp, "_Control")
    test_name <- paste0(exp, "_Test")
    control <- all_samples_Mvals[all_samples_Mvals$sample == cont_name, ]
    test <- all_samples_Mvals[all_samples_Mvals$sample == test_name, ]
    exp_samples <- merge(control, test, by = "feature_locs", all.x=FALSE, all.y=FALSE)

    # tidy exp_samples from this:
    # feature_locs meth.x unmeth.x   sample.x meth.y unmeth.y sample.y
    # to this:
    # feature_locs B3_Control.meth B3_Control.unmeth B3_Test.meth B3_Test.unmeth
    exp_samples$sample.x <- NULL
    exp_samples$sample.y <- NULL
    colnames(exp_samples) <- c(
        "feature_locs",
        paste0(cont_name, ".meth"), paste0(cont_name, ".unmeth"),
        paste0(test_name, ".meth"), paste0(test_name, ".unmeth")
    )
    # save exp_samples externally
    named_list_of_un_meth_data[[exp]] <- exp_samples
}

# all_samples <- merge(
#     named_list_of_un_meth_data[[c("B3")]], named_list_of_un_meth_data[[c("B4")]],
#     by = "feature_locs"
# )

# note un_meth means unmeth and meth data
print(paste0("making all_samples_un_meth"))
all_samples_un_meth = list()
for (exp in exp_list) {
    if (length(all_samples_un_meth) == 0) {
        all_samples_un_meth <- copy(named_list_of_un_meth_data[[exp]]) # copy otherwise it gets overridden
    } else {
        all_samples_un_meth <- merge(
            all_samples_un_meth, named_list_of_un_meth_data[[exp]],
            by = "feature_locs"
        )
        # named_list_of_un_meth_data[exp]
    }
}

# ======================================== #
# Code that went into the for loop above
# all_samples_Mvals_control<-all_samples_Mvals[all_samples_Mvals$sample=="B3_Control",]
# all_samples_Mvals_test<-all_samples_Mvals[all_samples_Mvals$sample=="B3_Test",]
# #lining up all positions accross the samples
# # all_samples_Mvals<-as.data.frame(all_samples_Mvals)

# Final_samples<-merge(all_samples_Mvals_test, all_samples_Mvals_control, by="feature_locs")
# colnames(Final_samples) <- c(
#     "feature_locs", "Strand.test", "meth_test", "unmeth_test", "B3_Test",
#     "strand_control","meth_control","unmeth_control","B3_Control"
# )

# ======================================== #

unmeth_data <- as.data.frame(all_samples_un_meth)[
    ,
    colnames(all_samples_un_meth) %in% grep("\\.unmeth", colnames(all_samples_un_meth),
        value = TRUE
    )
]

meth_data <- as.data.frame(all_samples_un_meth)[
    ,
    colnames(all_samples_un_meth) %in% grep("\\.meth", colnames(all_samples_un_meth),
        value = TRUE
    )
]

# Get feature_locs
feature_locs = all_samples_un_meth$feature_locs
feature_locs <- gsub("_", ":", feature_locs)
feature_locs <- gsub("\\.", "*", feature_locs)
new_colnames <- sub("\\.unmeth$", "", colnames(unmeth_data)) # save for later

ifelse(sub(
    "\\.unmeth$", "", colnames(unmeth_data)
) == sub(
    "\\.meth$", "", colnames(meth_data)
), "proceed with next code", "something screwed up stop and check")

# remove colnames for the Mset to work
colnames(meth_data) <- NULL
colnames(unmeth_data) <- NULL

gMset <- GenomicMethylSet(
    gr = as(feature_locs, "GRanges"),
    Meth = as.matrix(meth_data),
    Unmeth = as.matrix(unmeth_data)
)
# Mset <- MethylSet(
#     Meth = as.matrix(meth_data),
#     Unmeth = as.matrix(unmeth_data)
# )

filename_mset = paste0(
    path2_intermediate_data, "/",
    name_of_test, ".gMset.Rdata"
)
save(gMset, new_colnames, file=filename_mset)

Mvalue <- getM(gMset)
colnames(Mvalue) <- new_colnames # add colnames back

# Create the targets table
targets <- data.table(
    Targets = new_colnames,
    Treatment = ifelse(grepl("Control", new_colnames), 0, 1),
    Replicate = as.integer(factor(sub("_(.*)", "", new_colnames)))
)

# Aaron's code
PCA_plot <- function(Mvalue, df) {
    require(minfi)
    require(ggplot2)
    require(limma)
    mycols <- c("salmon", "skyblue", "mediumpurple", "orange", "blue", "grey", "pink")

    MDS_Treatment_Slide <- as.data.frame(plotMDS(Mvalue, top = 1000, gene.selection = "pairwise", dim = c(1, 2)))
    ggplot(MDS_Treatment_Slide, aes(x = x, y = y, color = as.factor(df$Treatment), shape = as.factor(df$Replicate))) +
        # scale_color_manual(values = c("#3B6AE1", "#C843C9")) +
        scale_color_manual(values = c("0" = "#0072B2", "1" = "#E69F00"), labels = c("0" = "Control", "1" = "Test")) +
        theme_bw() +
        geom_point(size = 7) +
        xlab(paste("PC1")) +
        ylab(paste("PC2")) +
        coord_fixed(ratio = 1, xlim = c(-2.5, 3.5), ylim = c(-2.5, 3.5)) +  # Set the same limits for both axes
        scale_x_continuous(breaks = seq(-2, 3, by = 1)) +  # Add ticks at every 1 unit on the x-axis
        scale_y_continuous(breaks = seq(-2, 3, by = 1)) +  # Add ticks at every 1 unit on the y-axis
        theme(
            axis.line = element_line(color = "black"),
            plot.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            axis.text.x = element_text(size = 20),
            axis.text.y = element_text(size = 20)
        ) +
        labs(color = "Treatment", shape = "Replicate")
}

# all_samples_Mvals_wide_complete$strand = NULL
PCA_plot(Mvalue, targets)

ggsave(
    filename = paste0(path2_plots, "/MDS_plot.", name_of_test, ".png"), dpi = 800, scale = 1.5,
    limitsize = FALSE
)

# Lucy's bad code
# # Transform to wide as a datatable
# all_samples_Mvals = as.data.table(all_samples_Mvals)
# all_samples_Mvals_wide = dcast(all_samples_Mvals, feature_locs + strand ~ sample, value.var = 'Mvals')

# # remove features with no vals in any of the rows
# all_samples_Mvals_wide = as.data.frame(all_samples_Mvals_wide)

# all_samples_Mvals_wide_complete <- all_samples_Mvals_wide[complete.cases(
#     all_samples_Mvals_wide[, !names(all_samples_Mvals_wide) %in% c("feature_locs")]
# ), ]

# all_samples_Mvals_wide_complete = as.data.table(all_samples_Mvals_wide_complete)
# # Note that there are infintite vals in this dataframe at this time, this means there was 0 methylation at that point
# all_samples_Mvals_wide_complete$strand = NULL

# all_samples_Mvals_wide_complete <- as.data.frame(all_samples_Mvals_wide_complete)

# rownames(all_samples_Mvals_wide_complete) <- all_samples_Mvals_wide_complete[, 1]

# # Extract the headers
# headers <- names(all_samples_Mvals_wide_complete)[3:ncol(all_samples_Mvals_wide_complete)]

# #
