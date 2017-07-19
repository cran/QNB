### R code from vignette source 'QNB-Overview.Rnw'

###################################################
### code chunk number 1: options
###################################################
options(width=60)


###################################################
### code chunk number 2: Input file
###################################################
library(QNB)
f1 = system.file("extdata", "control_ip.txt", package="QNB")
f2 = system.file("extdata", "treated_ip.txt", package="QNB")
f3 = system.file("extdata", "control_input.txt", package="QNB")
f4 = system.file("extdata", "treated_input.txt", package="QNB")

meth1 = read.table(f1, header=TRUE)
meth2 = read.table(f2, header=TRUE)
unmeth1 = read.table(f3, header=TRUE)
unmeth2 = read.table(f4, header=TRUE)
head(meth1)
head(unmeth1)


###################################################
### code chunk number 3: per-conditon
###################################################
result = qnbtest(meth1, meth2, unmeth1, unmeth2, mode="per-condition")
head(result)


###################################################
### code chunk number 4: provide size factor
###################################################
total_number_reads_control_ip <- c(3015921,2563976,198530)
total_number_reads_treated_ip <- c(1565101,152389,323569)
total_number_reads_control_input <- c(108561,302534,108123)
total_number_reads_treated_input <- c(301270,208549,308654)

standard_library_size <- exp(mean(log( c(total_number_reads_control_ip,
                                  total_number_reads_treated_ip,
                                  total_number_reads_control_input,
                                  total_number_reads_treated_input))))

size.factor <- list(control_ip = total_number_reads_control_ip/standard_library_size,
                   treated_ip = total_number_reads_treated_ip/standard_library_size,
                   control_input = 
                     total_number_reads_control_input/standard_library_size,
                   treated_input = 
                     total_number_reads_treated_input/standard_library_size)
                   
result <- qnbtest(meth1, meth2, unmeth1, unmeth2,
                 size.factor = size.factor)



###################################################
### code chunk number 5: without replicates
###################################################
f1 = system.file("extdata", "no_rep_controlip.txt", package="QNB")
f2 = system.file("extdata", "no_rep_treatedip.txt", package="QNB")
f3 = system.file("extdata", "no_rep_controlinput.txt", package="QNB")
f4 = system.file("extdata", "no_rep_treatedinput.txt", package="QNB")

no_rep_meth1 = read.table(f1, header=TRUE)
no_rep_meth2 = read.table(f2, header=TRUE)
no_rep_unmeth1 = read.table(f3, header=TRUE)
no_rep_unmeth2 = read.table(f4, header=TRUE)
head(no_rep_meth1)
head(no_rep_unmeth1)
result = qnbtest(no_rep_meth1, 
                 no_rep_meth2,
                 no_rep_unmeth1,
                 no_rep_unmeth2,
                 mode="blind")


###################################################
### code chunk number 6: automatically
###################################################
result = qnbtest(meth1, meth2, unmeth1, unmeth2)


###################################################
### code chunk number 7: The complete processing flow
###################################################
library(exomePeak)
library(QNB)
library(GenomicFeatures)
library(Rsamtools)

#peak calling using exomePeak
GENE_ANNO_GTF = system.file("extdata", "example.gtf", package="exomePeak")
f1 = system.file("extdata", "IP1.bam", package="exomePeak")
f2 = system.file("extdata", "IP2.bam", package="exomePeak")
f3 = system.file("extdata", "treated_IP1.bam", package="exomePeak")
IP_BAM = c(f1,f2,f3)
f4 = system.file("extdata", "Input1.bam", package="exomePeak")
f5 = system.file("extdata", "Input2.bam", package="exomePeak")
f6 = system.file("extdata", "treated_Input1.bam", package="exomePeak")
INPUT_BAM = c(f4,f5,f6)

res = exomepeak(GENE_ANNO_GTF=GENE_ANNO_GTF, IP_BAM=IP_BAM, INPUT_BAM=INPUT_BAM)



#get reads count
peak=res$all_peaks
untreated_ip=matrix(0,nrow=length(peak),ncol=2)
untreated_input=matrix(0,nrow=length(peak),ncol=2)
treated_ip=matrix(0,nrow=length(peak),ncol=1)
treated_input=matrix(0,nrow=length(peak),ncol=1)
txdb <- makeTxDbFromUCSC(genome="hg19")
exonRanges <- exonsBy(txdb, "tx")

#get ip reads count
#f1
aligns <- readGAlignments(f1)
para <- ScanBamParam(what="mapq")
mapq <- scanBam(f1, param=para)[[1]][[1]]
# filter reads with mapq smaller than 30. 
mapq[is.na(mapq)] <- 255  # Note: mapq "NA" means mapq = 255
ID_keep <- (mapq >30)
filtered <- aligns[ID_keep]
id <- countOverlaps(filtered,exonRanges)
transcriptome_filtered_aligns <- filtered[id>0]
counts <- countOverlaps(peak, transcriptome_filtered_aligns)
#counts <- countOverlaps(peak, filtered)
untreated_ip[,1] <- counts

#f2
aligns <- readGAlignments(f2)
para <- ScanBamParam(what="mapq")
mapq <- scanBam(f2, param=para)[[1]][[1]]
# filter reads with mapq smaller than 30. 
mapq[is.na(mapq)] <- 255  # Note: mapq "NA" means mapq = 255
ID_keep <- (mapq >30)
filtered <- aligns[ID_keep]
id <- countOverlaps(filtered,exonRanges)
transcriptome_filtered_aligns <- filtered[id>0]
counts <- countOverlaps(peak, transcriptome_filtered_aligns)
#counts <- countOverlaps(peak, filtered)
untreated_ip[,2] <- counts

#f3
aligns <- readGAlignments(f3)
para <- ScanBamParam(what="mapq")
mapq <- scanBam(f3, param=para)[[1]][[1]]
# filter reads with mapq smaller than 30. 
mapq[is.na(mapq)] <- 255  # Note: mapq "NA" means mapq = 255
ID_keep <- (mapq >30)
filtered <- aligns[ID_keep]
id <- countOverlaps(filtered,exonRanges)
transcriptome_filtered_aligns <- filtered[id>0]
counts <- countOverlaps(peak, transcriptome_filtered_aligns)
#counts <- countOverlaps(peak, filtered)
treated_ip[,1] <- counts

#f4
aligns <- readGAlignments(f4)
para <- ScanBamParam(what="mapq")
mapq <- scanBam(f4, param=para)[[1]][[1]]
# filter reads with mapq smaller than 30. 
mapq[is.na(mapq)] <- 255  # Note: mapq "NA" means mapq = 255
ID_keep <- (mapq >30)
filtered <- aligns[ID_keep]
id <- countOverlaps(filtered,exonRanges)
transcriptome_filtered_aligns <- filtered[id>0]
counts <- countOverlaps(peak, transcriptome_filtered_aligns)
#counts <- countOverlaps(peak, filtered)
untreated_input[,1] <- counts

#get input reads count
#f5
aligns <- readGAlignments(f5)
para <- ScanBamParam(what="mapq")
mapq <- scanBam(f5, param=para)[[1]][[1]]
# filter reads with mapq smaller than 30. 
mapq[is.na(mapq)] <- 255  # Note: mapq "NA" means mapq = 255
ID_keep <- (mapq >30)
filtered <- aligns[ID_keep]
id <- countOverlaps(filtered,exonRanges)
transcriptome_filtered_aligns <- filtered[id>0]
counts <- countOverlaps(peak, transcriptome_filtered_aligns)
#counts <- countOverlaps(peak, filtered)
untreated_input[,2] <- counts


#f6
aligns <- readGAlignments(f6)
para <- ScanBamParam(what="mapq")
mapq <- scanBam(f6, param=para)[[1]][[1]]
# filter reads with mapq smaller than 30. 
mapq[is.na(mapq)] <- 255  # Note: mapq "NA" means mapq = 255
ID_keep <- (mapq >30)
filtered <- aligns[ID_keep]
id <- countOverlaps(filtered,exonRanges)
transcriptome_filtered_aligns <- filtered[id>0]
counts <- countOverlaps(peak, transcriptome_filtered_aligns)
#counts <- countOverlaps(peak, filtered)
treated_input[,1] <- counts


#differential RNA methylation analysis
result = qnbtest(untreated_ip,treated_ip,untreated_input,treated_input)


###################################################
### code chunk number 8: session info
###################################################
sessionInfo()


