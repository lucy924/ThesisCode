# ThesisCode
Repository for all code referenced in my PhD thesis.  

**Naming**  
In this repository you will find code referencing slightly different than it appears in the thesis.  
* Cell line B and Cell line C refer to SW780 and RT4 respectively.  
* B3 and C3 refer to Donor 1 from SW780 and RT4 experiments.  
* B4 to Donor 2 (there is no data from Donor 2 in RT4 experiments)
* B5 and C5 to Donor 3  
* B6 and C6 to Donor 4

## Karyotyping
* The ACE package documentation can be found here: https://github.com/tgac-vumc/ACE  
* `get_karyotype_text.py` is a script to help obtain the correct karyotyping terminology from ACE output.  
* The predicted chromosomal modal number used as input to the ACE program was 4N. 
    * 2N, 3N and 5N were tried on B3B5_Control, however they return poor quality results.  
* There is one R script for each sample karyotype as they were run concurrently.  
* `Karyotyping.ipynb` contains development code and the karyotyping comparisons.  
* The `data` directory contains output from the ACE program for copy number calling. For raw ACE outputs please contact me or Dr. Stevens.
