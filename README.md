# Thesis Code
Repository for code referenced in my PhD thesis.  


## BEBIC Differential Modification Analysis (Chapter 3)
Some files here are example files for the SW780 test of Adenine outside of the AGG context. All other tests were performed similarly but with the relevant names and base references changed. The code is presented here in a manner easy to follow, with little repetition or messy development, however raw code for all Chapter 3 analyses is available [here](https://github.com/lucy924/BEBIC_DifferentialMethylation) if you wish to view it.  
Please note that many jupyter notebooks contain both python and R code.  
DMP = Differential Modification Position  
DMR = Differential Modification Region  

**Naming**  
In this repository you will find code referencing samples slightly differently than it appears in the thesis.  
* Cell line B and Cell line C refer to SW780 and RT4 respectively.  
* B3 and C3 refer to Donor 1 from SW780 and RT4 experiments.  
* B4 to Donor 2 (there is no data from Donor 2 in RT4 experiments)
* B5 and C5 refer to Donor 3  
* B6 and C6 refer to Donor 4

### Karyotyping
* The ACE package documentation can be found here: https://github.com/tgac-vumc/ACE  
* `get_karyotype_text.py` is a script to help obtain the correct karyotyping terminology from ACE output.  
* The predicted chromosomal modal number used as input to the ACE program was 4N. 
    * 2N, 3N and 5N were tried on B3B5_Control, however they return poor quality results.  
* There is one R script for each sample karyotype as they were run concurrently.  
* `Karyotyping.ipynb` contains development code and the karyotyping comparisons.  
* The `data` directory contains output from the ACE program for copy number calling. For raw ACE outputs please contact me or Dr. Stevens.

### Preprocessing
Pod5 files output from the P2 Solo were initially aligned during basecalling in one step to the hg38 genome.  
Later these bamfiles were stripped of their alignment and realigned to the CHM13-T2T genome in order to improve alignments across telomeres and centromeres. It was these realigned files that were used as input for [pileups.sh](Chapter3/bebic_analysis/pileups.sh) to generate pileup files for each test. This is an example pileup file for Adenine outside the AGG context.

### Genome-wide modification and MDS analysis
Sequence files were split into per chromosome per sample.  
For these genome wide analyses, data had to first be combined:  [basicstats.getchromdata.R](Chapter3/bebic_analysis/per_test_scripts/basicstats/basicstats.getchromdata.R)  
then run for genome wide methylation density plots: [basicstats.plot_all_sample_data.R](Chapter3/bebic_analysis/per_test_scripts/basicstats/basicstats.plot_all_sample_data.R)  
and finally for MDS plots: [basic_stats.plot_MDS.R](Chapter3/bebic_analysis/per_test_scripts/basicstats/basic_stats.plot_MDS.R).  
These were called sequentially by a bash script [run_basicstats.sh](Chapter3/per_test_scripts/bebic_analysis/basicstats/run_basicstats.sh).   

### Differential modification analysis - DMP
Both `modkit` and `minfi` were used for DMP analysis.  
* Modkit code is here: [run_modkit_dmr.sh](Chapter3/bebic_analysis/per_test_scripts/run_modkit_dmr.sh) and performs single site analysis (DMP) and segmentation analysis (DMR) where it dynamically searches for "segments" along the genome of differential methylation.
* Minfi code is here: [run_minfi_DMP.R](Chapter3/bebic_analysis/per_test_scripts/run_minfi_DMP.R)
* Results from modkit and minfi were restructured in preparation for plotting effect sizes and other plots.
    * modkit restructure: [filter_raw_dmps.py](Chapter3/bebic_analysis/per_test_scripts/filter_raw_dmps.py)
    * minfi restructure: [restructure_minfi_results.R](Chapter3/bebic_analysis/per_test_scripts/restructure_minfi_results.R)
* manhattan plots were first made from each set of data to detemine scoring parameter for plotting the effect sizes  
    * [plot_manhattan.modkit.py](Chapter3/bebic_analysis/per_test_scripts/plot_manhattan.modkit.py)  
    * [plot_manhattan.minfi.py](Chapter3/bebic_analysis/per_test_scripts/plot_manhattan.minfi.py)  
    * [plot_effect_sizes.modkit.py](Chapter3/bebic_analysis/per_test_scripts/plot_effect_sizes.modkit.py)  
    * [plot_effect_sizes.minfi.py](Chapter3/bebic_analysis/per_test_scripts/plot_effect_sizes.minfi.py)  
* each resulting csv file from above scripts was annotated with any DMPs that were found in the other method: [common_DMPs.py](Chapter3/bebic_analysis/per_test_scripts/common_DMPs.py)  


### Differential modification analysis - DMR
`modkit` and `DSS` were used for DMR analysis. As mentioned above, DMR analysis in modkit was performed concurrently with DMP analysis.  
`DSS` analysis, along with trialling other methods, are here: [run_other_DMR.alltests.R](Chapter3/bebic_analysis/run_other_DMR.alltests.R)  
Because there were so few DMR results, all the code from the related files above are condensed into one jupyter notebook: [segments_all_test.ipynb](Chapter3/bebic_analysis/segments_all_test.ipynb)

### Gene analyses
Code for DMP/DMR mapping to genes and ENCODE data, as well as Gene Ontology and KEGGG pathway analyses, can be found in [gene_analyses.ipynb](Chapter3/bebic_analysis/gene_analyses.ipynb) for 5mC and 5hmC, and in some cases for 6mA. Identified genes were input into GO and KEGG enrichment analyses using `goana` and `kegga` from the limma package in R.

#### Genes and ENCODE annotation
References used were from the CHM13-T2T gene and ENCODE annotation data, here: https://github.com/marbl/CHM13 (gene annotation version JHU RefSeqv110 + Liftoff v5.1).
Lolliplots were made to show were DMPs and DMRs are mapped to gene regions, and are located here: [Chapter3/bebic_analysis/lolliplots](Chapter3/bebic_analysis/lolliplots). Only relevant genes or locations with multiple mappings across methods (modkit/minfi/DSS) or cell lines (SW780/RT4) are presented in the main thesis body. 
