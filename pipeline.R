#library(Rsamtools)
#library(Rhisat2)
#library(Rsubread)

#unnamed <- function(fastq_path, genome_indexes, output){
  #hisat2(sequences = fastq_path, 
   #      index = genome_indexes,
    #     outfile = output)
  #Rsamtools::sortBam("output.sam", "output.bam")
  
#  Rsubread::featureCounts(output,
#                          annot.ext = "Exosome/Homo_sapiens.GRCh38.106.gtf",
#                          isGTFAnnotationFile = TRUE)
#}

#data <- unnamed("Exosome/data/SRR11412227.fastq", "Exosome/HISAT2/grch38/genome", "Exosome/data/SRR11412227.sam")

file <- read.table("Exosome/GSE147507_RawReadCounts_Human.tsv", sep="\t", header = TRUE)
