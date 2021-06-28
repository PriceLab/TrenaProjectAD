library(ghdb)
ghdb <- GeneHancerDB()
goi <- scan("Genes_toCheck_100_inFig3.txt", what=character(0), quiet=TRUE)
igv <- start.igv("RXRA", "hg38")

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

#------------------------------------------------------------------------------------------------------------------------
promotersWithOC <- function()
{
   for(gene in goi){
     tbl.gh <- retrieveEnhancersFromDatabase(ghdb, gene, tissues="all")
     if(nrow(tbl.gh) == 0) next;
     tbl.gh <- subset(tbl.gh, elite)
     loc.chrom <- sprintf("chr%s", tbl.gh$chrom[1])
     loc.start <- min(tbl.gh$start)
     loc.end <- max(tbl.gh$end)
     #tbl.gh <- subset(tbl.gh, combinedscore > 500)
     shoulder <- 1000
     tbl.atac <- subset(tbl.dba, seqnames==loc.chrom & start >= (loc.start - shoulder) & end <= (loc.end + shoulder))
     printf("%10s promoter regions: %d  oc: %d", gene, nrow(tbl.gh), nrow(tbl.atac))
     }



} # promotersWithOC
#------------------------------------------------------------------------------------------------------------------------
displayPeak <- function(loc.chrom, loc.start, loc.end)
{
   showGenomicRegion(igv, sprintf("%s:%d-%d", loc.chrom, loc.start, loc.end))
   zoomOut(igv)
   zoomOut(igv)
   removeTracksByName(igv, getTrackNames(igv)[-1])

   loc <- getGenomicRegion(igv)

   tbl.gh <- queryByRegion(ghdb, loc$chrom, loc$start, loc$end)
   dim(tbl.gh)

   tbl.gh$score <- asinh(tbl.gh$combinedscore)
   track <- DataFrameQuantitativeTrack("GH", tbl.gh[, c("chrom", "start", "end", "score")],
                                       autoscale=TRUE, color="brown")
   displayTrack(igv, track)

   for(file in peak.files){
       tbl <- read.table(file.path("../incoming", file), skip=1, sep="\t", as.is=TRUE)
       colnames(tbl) <- c("chrom", "start", "end", "name", "score")
       tbl.sub <- subset(tbl, chrom==loc$chrom & start >= loc$start & end <= loc$end)
       dim(tbl.sub)
       trackName <- substr(file, 1, 14)
       if(file %in% ad.files) color <- "blue"
       if(file %in% psp.files) color <- "red"
       if(file %in% ctl.files) color <- "green"
       track <- DataFrameQuantitativeTrack(trackName, tbl.sub[, c(1:3,5)], autoscale=TRUE, color=color,
                                           trackHeight=25, min=0, max=400)
       displayTrack(igv, track)
       Sys.sleep(0.5)
       } # for file
  Sys.sleep(30)
  } # for gene



} # displayPeak
#------------------------------------------------------------------------------------------------------------------------
