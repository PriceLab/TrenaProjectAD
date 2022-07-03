ad.files <- c("06_736_CER_ATAC-6.FCHFGJHBBXY_L5_R1_IACCACTGT.fastq.gz.PE_macs2_peaks.bed",
              "15_871_CER-ATAC-15.FCHFGJHBBXY_L6_R1_ITTGACCCT.fastq.gz.PE_macs2_peaks.bed",
              "21__135_CER_ATAC_21.FCHGG3MBBXY_L7_IGGACTCCT_R1.fastq.gz.PE_macs2_peaks.bed",
              "26__751_CER_ATAC_26.FCHGFLKBBXY_L2_R1_ICGAGGCTG.PE_macs2_peaks.bed",
              "36__1015_CER_ATAC_36.FCHGFLKBBXY_L3_R1_IGTGTGGTG.PE_macs2_peaks.bed",
              "45__1019_CER_ATAC_45.FCHGFLKBBXY_L4_R1_IGGACTCCT.PE_macs2_peaks.bed",
              "49__1166_CER_ATAC-49.FCHGFLKBBXY_L5_R1_IGCTACGCT.PE_macs2_peaks.bed",
              "51__966_CER_ATAC-51.FCHGFLKBBXY_L5_R1_IAAGAGGCA.PE_macs2_peaks.bed")

psp.files <- c("04_11460_CER_ATAC-4.FCHFGJHBBXY_L5_R1_IAAGAGGCA.fastq.gz.PE_macs2_peaks.bed",
               "11_11506_CER-ATAC-11.FCHFGJHBBXY_L6_R1_IAGGTTGGG.fastq.gz.PE_macs2_peaks.bed",
               "14_11259_CER-ATAC-14.FCHFGJHBBXY_L6_R1_ITGGTCACA.fastq.gz.PE_macs2_peaks.bed",
               "18__11374_CER_ATAC_18.FCHGG3MBBXY_L7_ICGTACTAG_R1.fastq.gz.PE_macs2_peaks.bed",
               "27__11258_CER_ATAC_27.FCHGFLKBBXY_L2_R1_IAAGAGGCA.PE_macs2_peaks.bed",
               "34__REPEAT3_ATAC_34.FCHGFLKBBXY_L3_R1_IGAGGGGTT.PE_macs2_peaks.bed",
               "39__11279_CER_ATAC_39.FCHGFLKBBXY_L3_R1_ITTGACCCT.PE_macs2_peaks.bed",
               "43__11491_CER_ATAC_43.FCHGFLKBBXY_L4_R1_IAGGCAGAA.PE_macs2_peaks.bed",
               "50__11497_CER_ATAC-50.FCHGFLKBBXY_L5_R1_ICGAGGCTG.PE_macs2_peaks.bed")

ctl.files <- c("20__11397_CER_ATAC_20.FCHGG3MBBXY_L7_ITCCTGAGC_R1.fastq.gz.PE_macs2_peaks.bed",
               "46__11400_CER_ATAC_46.FCHGFLKBBXY_L4_R1_ITAGGCATG.PE_macs2_peaks.bed",
               "35__11323_CER_ATAC_35.FCHGFLKBBXY_L3_R1_IAGGTTGGG.PE_macs2_peaks.bed",
               "31__11294_CER_ATAC_31.FCHGFLKBBXY_L2_R1_ITGGATCTG.PE_macs2_peaks.bed",
               "02_1942_CER_ATAC-2.FCHFGJHBBXY_L5_R1_IGGACTCCT.fastq.gz.PE_macs2_peaks.bed")

peak.files <- c(ad.files, psp.files, ctl.files)
length(peak.files) # 22
tbls <- list()

for(file in peak.files){
   printf("--- %s", file)
   tbl <- read.table(file.path("../incoming", file), skip=1, sep="\t", as.is=TRUE)
   if(file %in% ad.files) dx <- "AD"
   if(file %in% psp.files) dx <- "PSP"
   if(file %in% ctl.files) dx <- "CTL"
   tbl$dx <- dx
   colnames(tbl) <- c("chrom", "start", "end", "name", "score", "dx")
   tbls[[file]] <- tbl
   }

names(tbls) <- NULL
lapply(tbls, dim)

tbl <- do.call(rbind, tbls)
library(GenomicRanges)
gr <- GRanges(tbl)
gr <- sort(gr)

tbl.sorted <- as.data.frame(gr)
dim(tbl.sorted)
head(tbl.sorted)
colnames(tbl.sorted)[1] <- "chrom"
tbl.sorted$chrom <- as.character(tbl.sorted$chrom)

tbl.atac <- tbl.sorted[, c("chrom", "start", "end", "score", "dx")]
save(tbl.atac, file="~/github/TrenaProjectAD/explore/mayo-epigenetics/atac/mayoAllPeaks.1052789x4.RData")
lapply(tbl.atac, class)
dim(tbl.atac)
fivenum(tbl.atac$score)
gr.atac <- GRanges(tbl.atac)
gr.atac.reduced <- reduce(gr.atac)
tbl.atac.merged <- as.data.frame(gr.atac.reduced)
colnames(tbl.atac.merged)[1] <- "chrom"
tbl.atac.merged$chrom <- as.character(tbl.atac.merged$chrom)
dim(tbl.atac.merged)
dim(tbl.atac)
f <- "~/github/TrenaProjectAD/explore/mayo-epigenetics/atac/mayoAllPeaks.merged.96064x4.RData"
save(tbl.atac.merged, file=f)


#----------------------------------------------------------------------------------------------------
# create piled-up atac regions, first just a small example in around the NDUFS2 tss
#----------------------------------------------------------------------------------------------------

if(!exists("igv")){
    igv <- start.igv("NDUFS2", "hg38")
}

roi <- getGenomicRegion(igv)

for(DX in c("AD", "PSP", "CTL")){
   tbl <- subset(tbl.atac, chrom==roi$chrom & start >= roi$start & end <= roi$end & dx==DX)
   tbl.summedScores <- as.data.frame(GRanges(coverage(GRanges(tbl), weight="score")))
   tbl.summedScores <- subset(tbl.summedScores, score >= 1)
   colnames(tbl.summedScores)[1] <- "chrom"
   tbl.summedScores$chrom <- as.character(tbl.summedScores$chrom)
   track <- DataFrameQuantitativeTrack(DX,
                                       tbl.summedScores[-1, c("chrom", "start", "end", "score")],
                                       autoscale=FALSE, color="random", min=0, max=500)
   displayTrack(igv, track)
   }


#----------------------------------------------------------------------------------------------------
# create a single piled-up coverage, on tbl for each diagnosis
#----------------------------------------------------------------------------------------------------
for(DX in c("AD", "PSP", "CTL")){
   tbl <- subset(tbl.atac, dx==DX)
   tbl.summedScores <- as.data.frame(GRanges(coverage(GRanges(tbl), weight="score")))
   tbl.summedScores <- subset(tbl.summedScores, score >= 1)
   colnames(tbl.summedScores)[1] <- "chrom"
   tbl.summedScores$chrom <- as.character(tbl.summedScores$chrom)
   printf("rows for %s: %d", DX, nrow(tbl.summedScores))
   save(tbl.summedScores, file=sprintf("tbl.summedScoresATAC-seq-%s.RData", DX))
   }
# [1] rows for AD: 661056
# [1] rows for PSP: 786604
# [1] rows for CTL: 372016
