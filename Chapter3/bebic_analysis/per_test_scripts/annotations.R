library("GenomicRanges")
library("data.table")
library("AnnotationHub")
library("GenomicFeatures")
library("dplyr")
library("stringr")
library(org.Hs.eg.db)

################# Pick variables #################
name_of_test = "test2"
bebic = "BEBIC_B"
# bebic = "BEBIC_C"
modBase <- "C"
# modBase <- "A"
motif <- "CpG"
# motif <- "AGG"

tests_to_check <- c('test2', 'test4', 'test6', 'test7', 'test9', 'test10', 'test11', 'test12')
CpG_tests <- c('test2', 'test6')
AGG_tests <- c('test4', 'test7')
C_tests <- c('test9', 'test10')
A_tests <- c('test11', 'test12')
BEBIC_B_tests <- c('test2', 'test4', 'test9', 'test11')
BEBIC_C_tests <- c('test6', 'test7', 'test10', 'test12')

# up to test7

ad = "/external/analyses/lucy/bebic_pipeline"

################# Functions #################

get_gene_info <- function(subject_idx, keytype) {
    
    geneids = genes_txdb$gene_id[subject_idx]
    # gene_symbols <- mapIds(org.Hs.eg.db, keys = geneids, column = "SYMBOL", keytype = "SYMBOL")
    
    # Sometimes it can't find any gene descriptions 
    result <- try({
        gene_descriptions <- mapIds(org.Hs.eg.db, keys = geneids, column = "GENENAME", keytype = keytype)
    }, silent = TRUE)

    # Check if an error occurred
    if (inherits(result, "try-error")) {
        print("No entries found, making an empty dataframe")
        # Create a dataframe with gene IDs, symbols, and descriptions
        gene_info <- data.frame(
            Symbol = geneids,
            EntrezID = rep(NA, length(subject_idx)),
            Description = rep(NA, length(subject_idx))
        )
    } else {
        print("Success")
        gene_entrez <- mapIds(org.Hs.eg.db, keys = geneids, column = "ENTREZID", keytype = keytype)
        # Create a dataframe with gene IDs, symbols, and descriptions
        gene_info <- data.frame(
            Symbol = geneids,
            EntrezID = gene_entrez,
            Description = gene_descriptions
        )
    }

    gene_info <- gene_info %>%
        mutate(Description = ifelse(startsWith(Symbol, "LOFF"), "rRNA", Description))

    return(gene_info)
}


fill_na_from_df2 <- function(df1, df2) {
  df1 %>% mutate(across(everything(), ~ ifelse(is.na(.), df2[[cur_column()]], .)))
}


add_missing_data <- function(gene_info) {

    # Gene names not quite perfectly annotated in this database from the names I started with, hunted for missing data
    # Found missing data! redo the description column:
    missing_data <- c(
        "RNA5S1_21" = "100169751", "TAF11L5_24" = "646066",
        "SNAR-C3_6" = "100170226", "LOFF_G0010787" = "110255164",
        "LOFF_G0010788" = "110255163", "LOFF_G0010304" = "106631779", 
        "LOFF_G0010292" = "106631780", "LOFF_G0010368" = "106631778", 
        "LOC100996442" = "WASH9P", "LOC283788_1" = "283788", "MGC70870_2" = "403340", 
        "LOC105379016_2" = "105379016", "LOC124901878_11" = "124901878"
    )  # manual additions
    missing_descriptions <- mapIds(org.Hs.eg.db, keys = missing_data, column = "GENENAME", keytype = "ENTREZID")

    # Add missing data 
    gene_info$EntrezID[gene_info$Symbol %in% names(missing_data)] <-
        missing_data[gene_info$Symbol[gene_info$Symbol %in% names(missing_data)]]

    gene_info$Description[gene_info$EntrezID %in% names(missing_descriptions)] <-
        missing_descriptions[gene_info$EntrezID[gene_info$EntrezID %in% names(missing_descriptions)]]

    return(gene_info)
}


################# Import stuff #################
# Annotation file
annotations_T2T_fp = "/home/dejlu879/refs/T2T/chm13v2.0_RefSeq_Liftoff_v5.2.gff3.gz"
# Extract genes from the TxDb
annotations_T2T = txdbmaker::makeTxDbFromGFF(annotations_T2T_fp)
genes_txdb <- genes(annotations_T2T) # 'txdb' is TxDb object

# ENCODE folders
# encode_coverage_dir = "/home/dejlu879/refs/T2T/chm13_ENCODE/coverage"
encode_macs2peak_dir = "/home/dejlu879/refs/T2T/chm13_ENCODE/macs2_peak_bed"
# List all BED files in the directory
bed_files <- list.files(path = encode_macs2peak_dir, pattern = "\\.bed$", full.names = TRUE)

for (name_of_test in tests_to_check) {
    if (name_of_test == 'test9') {
        next  # test9 didn't have any dmps
    }
    # Set up variables
    if (name_of_test %in% BEBIC_B_tests) {
    bebic <- 'BEBIC_B'
    # exp_list <- exp_lists[[1]]
    } else if (name_of_test %in% BEBIC_C_tests) {
    bebic <- 'BEBIC_C'
    # exp_list <- exp_lists[[2]]
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

    print(paste0('Processing ', name_of_test, ': bebic ', bebic, ', base ', modBase, ', motif ', motif))

    ###################################################

    ################# Pipeline paths #################

    wd = paste0("/home/dejlu879/20240731-BEBIC_dmr/pipeline_", name_of_test)
    print(paste0("  working directory: ", wd))
    path2_intermediate_data = paste0(ad, "/", name_of_test, "/intermediate_data_files/")
    path2_output = paste0(wd, "/output_data/")
    path2_plots = paste0(wd, "/plots/")
    print(paste0("  output directory for final files: ", path2_output))
    if (!dir.exists(path2_output)) {
        dir.create(path2_output)
    }

    # Get my data, B is modkit, C is minfi
    path2Bdmr = paste0(path2_intermediate_data, "sig_DMPs.", bebic, ".modkit.noted_minfi.csv")
    path2Cdmr = paste0(path2_intermediate_data, "sig_DMPs.", bebic, ".minfi.noted_modkit.csv")

    Bdmr_results = fread(path2Bdmr)
    Cdmr_results = fread(path2Cdmr)

    if ("end" %in% colnames(Bdmr_results)) {
        print("  End column exists for modkit data")
    } else {
        print("  End column does not exist, being made for minfi data.")
        Bdmr_results$end = Bdmr_results$start + 1
    }
    if ("end" %in% colnames(Cdmr_results)) {
        print("  End column exists for minfi data")
    } else {
        print("  End column does not exist, being made for modkit data.")
        Cdmr_results$end = Cdmr_results$start + 1
    }

    setnames(Bdmr_results, "balanced_map_pvalue", "pval")
    Bdmr_results <- subset(Bdmr_results, select = c(chrom, start, end, balanced_effect_size, score, pval, name, also_in_other))
    Cdmr_results <- subset(Cdmr_results, select = c(chrom, start, end, balanced_effect_size, score, pval, name, also_in_other))

    # Make GRanges from my data
    B_gr = makeGRangesFromDataFrame(Bdmr_results, keep.extra.columns = TRUE)
    C_gr = makeGRangesFromDataFrame(Cdmr_results, keep.extra.columns = TRUE)

    ################# ANNOTATE #################
    print("  Annotating...")

    # Find overlaps between GRanges object and the genes
    B_overlaps <- findOverlaps(B_gr, genes_txdb)
    C_overlaps <- findOverlaps(C_gr, genes_txdb)

    # Extract the indices of overlaps
    B_query_idx <- queryHits(B_overlaps)        # Indices from B_gr
    B_subject_idx <- subjectHits(B_overlaps)    # Indices from genes_txdb

    C_query_idx <- queryHits(C_overlaps)        # Indices from C_gr
    C_subject_idx <- subjectHits(C_overlaps)    # Indices from genes_txdb


    # Get gene symbols, ids, and descriptions

    # geneids = genes_txdb$gene_id[B_subject_idx]
    # gene_descriptions <- mapIds(org.Hs.eg.db, keys = geneids, column = "GENENAME", keytype = "SYMBOL")
    # gene_entrez <- mapIds(org.Hs.eg.db, keys = geneids, column = "ENTREZID", keytype = "SYMBOL")

    # # Create a dataframe with gene IDs, symbols, and descriptions
    # gene_info <- data.frame(
    #     Symbol = geneids,
    #     EntrezID = gene_entrez,
    #     Description = gene_descriptions
    # )

    gene_info_CS = get_gene_info(C_subject_idx, "SYMBOL")
    gene_info_CA = get_gene_info(C_subject_idx, "ALIAS")
    gene_info_C <- fill_na_from_df2(gene_info_CS, gene_info_CA)

    gene_info_BS = get_gene_info(B_subject_idx, "SYMBOL")
    gene_info_BA = get_gene_info(B_subject_idx, "ALIAS")
    gene_info_B <- fill_na_from_df2(gene_info_BS, gene_info_BA)

    # missing_descriptions <- mapIds(org.Hs.eg.db, keys = missing_data, column = "GENENAME", keytype = "ENTREZID")

    # # Add missing data 
    # gene_info$EntrezID[gene_info$Symbol %in% names(missing_data)] <-
    #     missing_data[gene_info$Symbol[gene_info$Symbol %in% names(missing_data)]]

    # gene_info$Description[gene_info$EntrezID %in% names(missing_descriptions)] <-
    #     missing_descriptions[gene_info$EntrezID[gene_info$EntrezID %in% names(missing_descriptions)]]

    gene_info_B <- add_missing_data(gene_info_B)
    gene_info_C <- add_missing_data(gene_info_C)

    # Step 2: Add gene annotations from genes_txdb to B_gr and C_gr

    # BEBIC B
    mcols(B_gr)$gene_id <- NA  # Initialize the gene_id column with NA in B_gr
    mcols(B_gr)$gene_id[B_query_idx] <- genes_txdb$gene_id[B_subject_idx]
    mcols(B_gr)$gene_id_strand <- NA  # Initialize the gene_id column with NA in B_gr
    mcols(B_gr)$gene_id_strand[B_query_idx] <- strand(genes_txdb)[B_subject_idx]

    match_idx <- match(mcols(B_gr)$gene_id, gene_info_B$Symbol) # Match gene_id column to Symbol column
    mcols(B_gr)$gene_description <- NA # Initialize the gene_id column with NA in B_gr
    mcols(B_gr)$gene_description <- gene_info_B$Description[match_idx]
    mcols(B_gr)$EntrezID <- NA # Initialize the EntrezID column with NA in B_gr
    mcols(B_gr)$EntrezID <- gene_info_B$EntrezID[match_idx]

    # BEBIC C
    mcols(C_gr)$gene_id <- NA  # Initialize the gene_id column with NA in C_gr
    mcols(C_gr)$gene_id[C_query_idx] <- genes_txdb$gene_id[C_subject_idx]
    mcols(C_gr)$gene_id_strand <- NA  # Initialize the gene_id column with NA in C_gr
    mcols(C_gr)$gene_id_strand[C_query_idx] <- strand(genes_txdb)[C_subject_idx]

    match_idx <- match(mcols(C_gr)$gene_id, gene_info_C$Symbol) # Match gene_id column to Symbol column
    mcols(C_gr)$gene_description <- NA # Initialize the gene_id column with NA in B_gr
    mcols(C_gr)$gene_description <- gene_info_C$Description[match_idx]
    mcols(C_gr)$EntrezID <- NA # Initialize the EntrezID column with NA in B_gr
    mcols(C_gr)$EntrezID <- gene_info_C$EntrezID[match_idx]

    print('  ...done.')

    ################# ANNOTATE ENCODE #################

    print('  Annotating ENCODE...')
    # Check for overlaps in ENCODE macs2_peak data
    for (file in bed_files) {
        # Extract the filename without the path
        filename <- basename(file)

        # Extract the "name" part from the filename (between "_" and ".")
        name <- str_extract(filename, "(?<=_)[^_]+(?=\\.)")

        # Print the extracted name (for debugging purposes)
        print(paste("  Processing file:", filename, "with name:", name))

        # Load the BED file using bedtorch
        # bed_data <- read_bed(file)

        bed_data <- fread(file)
        bed_data = makeGRangesFromDataFrame(
            bed_data, keep.extra.columns = TRUE, seqnames.field = c('chrom'), start.field = c('chromStart'), end.field = c('chromEnd'), strand.field = c('strand')
            )

        overlaps <- findOverlaps(B_gr, bed_data)
        # Extract the indices of overlaps
        query_idx <- queryHits(overlaps)        # Indices from B_gr
        subject_idx <- subjectHits(overlaps) # Indices from genes_txdb
        mcols(B_gr)[[name]] <- NA # Initialize the encode_name column with NA in B_gr
        mcols(B_gr)[[name]][query_idx] <- str_extract(bed_data$name[subject_idx], "peak_\\d+$")

        overlaps <- findOverlaps(C_gr, bed_data)
        # Extract the indices of overlaps
        query_idx <- queryHits(overlaps)        # Indices from B_gr
        subject_idx <- subjectHits(overlaps) # Indices from genes_txdb
        mcols(C_gr)[[name]] <- NA # Initialize the encode_name column with NA in B_gr
        mcols(C_gr)[[name]][query_idx] <- str_extract(bed_data$name[subject_idx], "peak_\\d+$")

    }

    print('  ...done.')
    ################# SAVE DATA #################

    # Now B_gr and C_gr have the gene annotations from genes_txdb
    # also include non-annotations

    b_fn = paste0(path2_output, "sig_DMPs.", bebic, ".", name_of_test, ".annotated_all.modkit.csv")
    write.csv(B_gr, file = b_fn)
    c_fn = paste0(path2_output, "sig_DMPs.", bebic, ".", name_of_test, ".annotated_all.minfi.csv")
    write.csv(C_gr, file = c_fn)

    print(paste0('  ... Modkit data is saved here: ', b_fn))
    print(paste0('  ... Minfi data is saved here: ', c_fn))

}
