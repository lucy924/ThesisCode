# #!/usr/bin/env Rscript

# Load necessary library for argument parsing (optional)
library(optparse)
library("data.table")
library(ggplot2)

compute_density <- function(data2densify, sample_name, bw_val = FALSE) {
    # computes density per sample
    # data2densify is either the percent_mods list or the betavals list
    if (bw_val) {
        dens <- density(data2densify, bw = bw_val)
    } else {
        dens <- density(data2densify)
    }
    data.table(x = dens$x, y = dens$y, sample = sample_name)
}


# Define command line options
option_list <- list(
    make_option(c("--bebic"),
        type = "character", default = NULL, help = "bebic, eg BEBIC_B"
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
    make_option(c("--path2_intermediate_data"),
        type = "character", default = NULL, help = "path to intermediate data fp"
    ),
    make_option(c("--bw_val"),
        type = "integer", default = NULL, help = "smoothing parameter, default None/auto"
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
bw_val <- options$bw_val

path2plots <- paste0(getwd(), "/plots/")

# ad <- "/external/analyses/lucy/bebic_pipeline"
# path2_intermediate_data = paste0(ad, "/", name_of_test, "/", "intermediate_data_files")

chr_names <- c(
    "chr1", "chr2", "chr3", "chr4", "chr5",
    "chr6", "chr7", "chr8", "chr9", "chr10",
    "chr11", "chr12", "chr13", "chr14", "chr15",
    "chr16", "chr17", "chr18", "chr19", "chr20",
    "chr21", "chr22", "chrX", "chrY", "chrM"
)

# Define custom color mapping for each sample
# colour_mapping <- c(
#     "B3_Control" = "#b9647a", # Add more samples with colors as needed
#     "B3_Test" = "#E3C1CA", # Example for additional samples
#     "B4_Control" = "#c843c9", # Add more samples with colors as needed
#     "B4_Test" = "#E9B5E9", # Example for additional samples
#     "B5_Control" = "#d34143", # Add more samples with colors as needed
#     "B5_Test" = "#EDB3B4", # Example for additional samples
#     "B6_Control" = "#c059a3", # Add more samples with colors as needed
#     "B6_Test" = "#E6BCDA", # Example for additional samples
#     "C3_Control" = "#3b6ae1",
#     "C3_Test" = "#BEE4E1",
#     "C5_Control" = "#5cbab4",
#     "C5_Test" = "#B1C4F3",
#     "C6_Control" = "#6c8db3",
#     "C6_Test" = "#C4D1E0"
# )

colour_mapping <- c(
    "B3_Control" = "#e67f27", # Add more samples with colors as needed
    "B3_Test" = "#e67f2780", # Example for additional samples
    "B4_Control" = "#ae0056", # Add more samples with colors as needed
    "B4_Test" = "#ae005680", # Example for additional samples
    "B5_Control" = "#00740f", # Add more samples with colors as needed
    "B5_Test" = "#00740f80", # Example for additional samples
    "B6_Control" = "#68afff", # Add more samples with colors as needed
    "B6_Test" = "#68afff80", # Example for additional samples
    "C3_Control" = "#e67f27",
    "C3_Test" = "#e67f2780",
    "C5_Control" = "#00740f",
    "C5_Test" = "#00740f80",
    "C6_Control" = "#68afff",
    "C6_Test" = "#68afff80"
)

# ["#00740f",  # green D3
# "#ae0056",  # pink D2
# "#68afff",  # blue D4
# "#e67f27"]  # orange D1

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

# ============================================ #
# Me working stuff out

# sample_name <- "B3_Control"


# save_chrom_info_fp <- paste0(
#     path2_intermediate_data, "/",
#     sample_name, ".", chrom, ".", name_of_test, ".meth_data.Rdata"
# )

# ================================ #
# How to get from this:
# save(
#     meth, unmeth, beta_vals_list, modBase_percent_modified, modBase_seqnames,
#     modBase_start, modBase_end, feature_locs, modBase_strand,
#     file = save_chrom_info_fp
# )

# To this:
# r$> head(all_percent_density_data)
#            x            y     sample
#        <num>        <num>     <char>
# 1: -5.750146 0.0004625781 B3_Control
# 2: -5.531946 0.0006492134 B3_Control
# 3: -5.313746 0.0008975092 B3_Control
# 4: -5.095546 0.0012225948 B3_Control
# 5: -4.877345 0.0016416162 B3_Control
# 6: -4.659145 0.0021735396 B3_Control

# ?????????????

# compute density data, this saves RAM by doing it separately then plotting it after
# origin
compute_percent_density <- function(dt, sample_name, bw_val) {
    sample_data <- dt[dt$sample == sample_name, ]
    dens <- density(sample_data$percent_modified, bw = bw_val)
    data.table(x = dens$x, y = dens$y, sample = sample_name)
}

compute_density <- function(data2densify, sample_name, bw_val = FALSE) {
    # computes density per sample
    # data2densify is either the percent_mods list or the betavals list
    if (bw_val) {
        dens <- density(data2densify, bw = bw_val)
    } else {
        dens <- density(data2densify)
    }
    data.table(x = dens$x, y = dens$y, sample = sample_name)
}

for (sample_name in sample_names_list) {
    print(sample_name)
    sample_WG_modBase_percent_modified <- list() # whole genome data
    # sample_WG_modBase_beta_modified <- list() # whole genome data

    for (chrom in chr_names) {
        save_chrom_info_fp <- paste0(
            path2_intermediate_data, "/",
            sample_name, ".", chrom, ".", name_of_test, ".meth_data.Rdata"
        )
        load(file = save_chrom_info_fp)
        sample_WG_modBase_percent_modified[[chrom]] <- modBase_percent_modified
        # sample_WG_modBase_beta_modified[[chrom]] <- beta_vals_list
    }
    remove(
        meth, unmeth, modBase_percent_modified, feature_locs
    )
    suppressMessages(gc())

    # flatten out list - for RAM management
    sample_WG_modBase_percent_modified <- unlist(sample_WG_modBase_percent_modified)
    # sample_WG_modBase_beta_modified <- unlist(sample_WG_modBase_beta_modified)

    if (motif == "A") {
        density_perc_dt <- compute_density(sample_WG_modBase_percent_modified, sample_name, bw_val = 1.0)
        # density_beta_dt <- compute_density(sample_WG_modBase_beta_modified, sample_name, bw_val = 0.003)
    } else if (motif == "C") {
        density_perc_dt <- compute_density(sample_WG_modBase_percent_modified, sample_name, bw_val = 1.0)
        # density_beta_dt <- compute_density(sample_WG_modBase_beta_modified, sample_name, bw_val = 0.003)
    } else {
        if (!is.null(bw_val)) {
            density_perc_dt <- compute_density(sample_WG_modBase_percent_modified, sample_name, bw_val = bw_val)
            # density_beta_dt <- compute_density(sample_WG_modBase_beta_modified, sample_name)
        } else {
            density_perc_dt <- compute_density(sample_WG_modBase_percent_modified, sample_name)
        }
    }

    # save density data
    perc_sample_density_fp <- paste0(
        path2_intermediate_data, "/",
        sample_name, ".", name_of_test, ".perc_density.Rdata"
    )
    save(density_perc_dt, file = perc_sample_density_fp)
}

all_percent_density_data <- rbindlist(lapply(
    sample_names_list, function(sample_name) {
        load(paste0(
            path2_intermediate_data, "/",
            sample_name, ".", name_of_test, ".perc_density.Rdata"
        ))
        return(density_perc_dt)
        # loads density_perc_dt
    }
))
remove(density_perc_dt)

# ================================ #
# Plot all samples percentage data
if (motif == "C") {
    basename <- "Cytosine"
} else if (motif == "A") {
    basename <- "Adenine"
} else {
    basename <- paste0(motif, " site")
# } else if (motif == "AGG") {
#     basename <- "AGG site"
}

# Replace 0 or negative densities with a small positive value (e.g., 1e-10)
all_percent_density_data$y[all_percent_density_data$y <= 0] <- 1e-10

# Add a new column to specify Control or Test
# all_percent_density_data$sample_type <- ifelse(grepl("Control", all_percent_density_data$sample), "Control", "Test")

# Plot the percent density data with custom color mapping and line types
plot_density <- ggplot(all_percent_density_data, aes(x = x, y = y, color = sample)) +  #, linetype = sample_type)) +
    geom_line() + 
    theme_bw() + # Set white background with grey grid
    # theme(
    #     panel.grid.major = element_line(color = "grey80"), # Customize major grid lines
    #     # panel.grid.minor = element_line(color = "grey90")  # Customize minor grid lines
    # )  + 
    # scale_y_log10(limits = c(1e-04, NA)) + # Cap the lower values at 1e-08
    scale_y_log10() + # No capping of values on A data
    scale_color_manual(values = colour_mapping) + # Apply custom colors
    # scale_linetype_manual(values = c("Control" = "solid", "Test" = "dashed")) + # Set line types
    labs(
        title = paste0("Genome Wide ", basename, " Percentage Modification"),
        x = "Percent modified", y = "Density"
    )
filename <- paste0(
    path2plots, "/",
    bebic, ".", name_of_test, ".perc_density.png"
)
ggsave(filename, plot_density, width = 8, height = 6, dpi = 300)
