library(ACE)
library(QDNAseq.hg38)
library(repr)

options(repr.plot.width = 12, repr.plot.height = 7)


# specify the directory containing your bam-files
path2alignedbam <- "/mnt/uow-pathology-dun/lucy/analyses/BEBIC_analysis/CellLineC/demuxed/by_sample/C6_Test/"
output_dir <- "/external/analyses/20241104-BEBIC_Karyotyping/ACE/C6_Test"
dir.create(output_dir)

# if you do not want the output in the same directory, use argument outputdir
runACE(path2alignedbam,
    filetype = "bam",
    ploidies = c(4), imagetype = "png", genome = "hg38", outputdir = output_dir
)
