suppressMessages(library("GenomicRanges"))
suppressMessages(library("tidyr"))
suppressMessages(library("data.table"))
suppressMessages(library("dplyr"))
suppressMessages(library("minfi"))
suppressMessages(library("DSS"))
# suppressMessages(library("DMRcate"))
# suppressMessages(library("edgeR"))


get_vars <- function(name_of_test) {
    if (name_of_test %in% BEBIC_B_tests) {
    bebic <- 'BEBIC_B'
    exp_list <- c("B3", "B4", "B5", "B6")
    } else if (name_of_test %in% BEBIC_C_tests) {
    bebic <- 'BEBIC_C'
    exp_list <- c("C3", "C5", "C6")
    } else {
    stop('test not associated with a bebic')
    }

    if (name_of_test %in% CpG_tests) {
            modBase <- 'C'
            motif <- 'CpG'
        } else if (name_of_test %in% AGG_tests) {
            modBase <- 'A'
            motif <- 'AGG'
        } else if (name_of_test %in% C_tests) {
            modBase <- 'C'
            motif <- 'C'
        } else if (name_of_test %in% A_tests) {
            modBase <- 'A'
            motif <- 'A'
        } else {
            stop('test not associated with a motif')
    }

    # if (name_of_test %in% modkit_score_8) {
    #         modkit_score <- 8
    #     } else if (name_of_test %in% modkit_score_1) {
    #         modkit_score <- 1
    #     } else if (name_of_test %in% modkit_score_6) {
    #         modkit_score <- 6
    #     } else if (name_of_test %in% modkit_score_2) {
    #         modkit_score <- 2
    #     } else {
    #         stop('test not associated with a score')
    # }

    # if (name_of_test %in% minfi_score_50) {
    #         minfi_score <- 50
    #     } else if (name_of_test %in% minfi_score_100) {
    #         minfi_score <- 100
    #     } else if (name_of_test %in% minfi_score_10) {
    #         minfi_score <- 10
    #     } else if (name_of_test %in% minfi_score_20) {
    #         minfi_score <- 20
    #     } else {
    #         stop('test not associated with a score')
    # }

    path2_output = paste0("/home/dejlu879/20240731-BEBIC_dmr/pipeline_", name_of_test, "/output_data/")

    return(list(bebic = bebic, modBase = modBase, motif = motif, path2_output = path2_output, exp_list = exp_list))
}


chr_names <- c(
    "chr1", "chr2", "chr3", "chr4", "chr5",
    "chr6", "chr7", "chr8", "chr9", "chr10",
    "chr11", "chr12", "chr13", "chr14", "chr15",
    "chr16", "chr17", "chr18", "chr19", "chr20",
    "chr21", "chr22", "chrX", "chrY", "chrM"
)

tests_to_check <- c('test2', 'test4', 'test6', 'test7', 'test10', 'test11', 'test12')

CpG_tests <- c('test2', 'test6')
AGG_tests <- c('test4', 'test7')
C_tests <- c('test9', 'test10')
A_tests <- c('test11', 'test12')
BEBIC_B_tests <- c('test2', 'test4', 'test9', 'test11')
BEBIC_C_tests <- c('test6', 'test7', 'test10', 'test12')

# modkit_score_8 <- c('test2', 'test6')
# modkit_score_1 <- c('test10')
# modkit_score_6 <- c('test4')
# modkit_score_2 <- c('test7', 'test11', 'test12')
# minfi_score_50 <- c('test2')
# minfi_score_100 <- c('test6')
# minfi_score_10 <- c('test4', 'test10')
# minfi_score_20 <- c('test7', 'test11', 'test12')

wd <- getwd()
print(paste0("working directory: ", wd))
ad <- "/external/analyses/lucy/bebic_pipeline"


log_level = ""
for (name_of_test in tests_to_check) {
    print(paste0("------ ", name_of_test, " ------"))
    log_level = "  "
    
    ################# Pick variables #################
    vars = get_vars(name_of_test)
    bebic = vars$bebic
    exp_list = vars$exp_list
    modBase = vars$modBase
    motif = vars$motif
    path2_output = vars$path2_output
    # modkit_score = vars$modkit_score
    # minfi_score = vars$minfi_score

    ###################################################
    
    ################# Pipeline paths #################

    path2_intermediate_data <- paste0(ad, "/", name_of_test, "/intermediate_data_files")
    path2_plots <- paste0(getwd(), "/plots/")
    output_dir <- paste0(path2_intermediate_data, "/")

    print(paste0(log_level, "output directory for intermediate files: ", output_dir))
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
    }

    print(paste0(log_level, "output directory for plots: ", path2_plots))
    if (!dir.exists(path2_plots)) {
        dir.create(path2_plots)
    }

    ###################################################

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
    
    # Extract the sample names
    sample_names <- new_colnames
    # Add the dimnames attribute to the MethylSet object
    dimnames(gMset) <- list(NULL, sample_names)
    # Initialize an empty list to store the data frames
    df_list <- list()

    
    # Loop through each sample and create the data frame
    for (sample in sample_names) {
        # Extract the necessary data from the methylset
        chromosome <- as.character(seqnames(gMset))
        position <- start(gMset)
        coverage <- minfi::getMeth(gMset)[, sample] + minfi::getUnmeth(gMset)[, sample]
        methylation <- minfi::getMeth(gMset)[, sample]

        # Create the data frame with the required columns and column names
        df <- data.frame(chr = chromosome, pos = position, N = coverage, X = methylation)
          # Add the data frame to the list with the sample name as the key
        df_list[[sample]] <- df
    }

    # Display the first few rows of the data frame
    # head(df_list[[1]])

    # Get group lists
    control_list <- c()
    test_list <- c()
    for (exp in exp_list) {
        control_list <- c(control_list, paste0(exp, "_Control"))
        test_list <- c(test_list, paste0(exp, "_Test"))
    }

    BS_obj_for_DMR = makeBSseqData(df_list, sample_names)
    data_for_DSS = DMLtest(BS_obj_for_DMR, group1=control_list, group2=test_list)
    dmrs <- callDMR(data_for_DSS, p.threshold = 0.05)

    # Filter to effect size >= 5%
    dmrs_effectsize5 <- dmrs[(dmrs$diff.Methy <= -0.05 | dmrs$diff.Methy >= 0.05), ]
    
    if (is.null(dmrs_effectsize5)) {
        print(paste0(log_level, "No DMRs."))
        flush.console()
        next
    }
    

    ###################################
    print("saving result...")
    
    # Rename the columns
    df <- dmrs_effectsize5 %>%
        rename(
            chrom = chr,
            N_sites = nCG,
            samplea_percents = meanMethy1,
            sampleb_percents = meanMethy2,
            effect_size = diff.Methy
        )

    # Create the new 'location' column
    df <- df %>%
        mutate(location = paste(chrom, start, sep = ":"))
    df <- df %>%
        mutate(location = paste(location, end, sep = "-"))

    # output to rdata
    # save(dmrs_effectsize5, )
    filename = paste0(path2_output, bebic, ".", motif, ".", name_of_test, ".minfiDSS_dmr.csv")
    write.csv(df, filename, row.names = FALSE)
    

    print("done.")

}

# DMRcate didn't return any "individually significant CpGs", used DSS instead
