library(DiffBind)

# following the examples in ?
mcf7 <- dba.peakset(NULL,
                    peaks=system.file("extra/peaks/MCF7_ER_1.bed.gz", package="DiffBind"),
                    peak.caller="bed", sampID="MCF7.1",tissue="MCF7",
                    factor="ER",condition="Responsive",replicate=1)


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
ad.dba <- NULL

for(i in seq_len(length(ad.files))){
   ad.dba <- dba.peakset(ad.dba, peaks=file.path("../incoming", ad.files[i]),
                         peak.caller="bed", sampID=sprintf("ad.%d", i), tissue="cerebellum",
                         factor="AD", condition="AD", replicate=i)
   }


for(i in seq_len(length(psp.files))){
   ad.dba <- dba.peakset(ad.dba, peaks=file.path("../incoming", psp.files[i]),
                         peak.caller="bed", sampID=sprintf("psp.%d", i), tissue="cerebellum",
                         factor="PSP", condition="psp", replicate=i)
   }

for(i in seq_len(length(ctl.files))){
   ad.dba <- dba.peakset(ad.dba, peaks=file.path("../incoming", ctl.files[i]),
                         peak.caller="bed", sampID=sprintf("ctl.%d", i), tissue="cerebellum",
                         factor="Control", condition="ctl", replicate=i)
   }


# When adding several peaksets via successive calls to ‘dba.peakset’, it may be
# more efficient to set this parameter to ‘FALSE’ and call ‘dba(DBA)’ after all
# of the peaksets have been added.

ad.dba <- dba(ad.dba)
ad.dba <- dba.contrast(ad.dba, design=TRUE)

ad.dba <- dba.contrast(ad.dba, design=TRUE, contrast=c("psp", "ctl"))


   atacDBA <- dba.peakset(atacDBA, sampID=sampIDs[i], peaks=consensus,
                       factor="ATAC",condition=condition[i],
                       treatment=treatment[i],replicate=replicate[i],
                       counts=counts[,i])

dba.plotMA(ad.dba)
dba.plotVolcano(ad.dba)

dba.plotPCA(ad.dba, DBA_CONDITION, label=DBA_REPLICATE)


#   file <- "~/github/TrenaProjectErythropoiesis/inst/extdata/genomicRegions/tbl.corces.hg38.day0thru4.RData"
#   file.exists(file)
#   load(file)
#   consensus <- GRanges(tbl.corces.hg38[,1:3])                          # all regions
#   counts    <- as.matrix(tbl.corces.hg38[,4:ncol(tbl.corces.hg38)])    # 592665      7
#   sampIDs   <- colnames(counts)                                        # "day0.1" "day0.2" "day0.3" "day0.4" "day2.1" "day2.2" "day4.1"
#   ids       <- strsplit(sampIDs,".",fixed=TRUE)
#   condition <- unlist(lapply(ids,function(x)x[1]))                     # "day0" "day0" "day0" "day0" "day2" "day2" "day4"
#   treatment <- condition
#   treatment[condition != "day0"] <- "dayPlus"                          # "day0"    "day0"    "day0"    "day0"    "dayPlus" "dayPlus" "dayPlus"
#   replicate <- as.numeric(unlist(lapply(ids,function(x)x[2])))         # 1 2 3 4 1 2 1
#
#   atacDBA <- NULL
#   for(i in 1:ncol(counts)) {
#     cat(i,sampIDs[i],"\n")
#     atacDBA <- dba.peakset(atacDBA, sampID=sampIDs[i], peaks=consensus,
#                            factor="ATAC",condition=condition[i],
#                            treatment=treatment[i],replicate=replicate[i],
#                            counts=counts[,i])
#     } # for
#
#   atacDBA$config$doBlacklist <- atacDBA$config$doGreylist <- FALSE
#   atacDBA
#   atacDBA <- dba.count(atacDBA,peaks=NULL,minCount=1,filter=5)
#   atacDBA
#   dba.plotPCA(atacDBA, DBA_CONDITION, label=DBA_REPLICATE)
#   atacDBA <- dba.contrast(atacDBA)
#   atacDBA <- dba.analyze(atacDBA)
#   atacDBA$config$th <- 0.01
#   dba.show(atacDBA, bContrasts = TRUE)
#   dba.save(atacDBA, file="shannon_atac")
#   dba.plotMA(atacDBA)
#   dba.plotVolcano(atacDBA)
#   dba.plotPCA(atacDBA, contrast=1, attributes=DBA_TREATMENT, label=DBA_CONDITION)
#   chrom.diff <- dba.report(atacDBA)  # 13483
#   tbl.rory <- as.data.frame(chrom.diff)
#   save(tbl.rory, file="tbl.diffBind.rory.hg38.day0-4.RData")

   # Note the filtering step, which removed any sites where no sample had at least five (noramlized) reads; this eliminated about 100k sistes from further analysis.
   # Regarding masks, here are some notes:
   #
   #  I'm not sure what pk.retrieved looks like? it is not the same object that you made your mask for. A mask is
   #  specific for a DBA object, you can't re-use a mask unless the object contains exactly the same samples.
   #
   #   Actually you almost never have to call dba.mask(). There are premade masks for each DBA object,
   #  for example peakset$masks$t0 and pk.retrieved$masks$t0. Look at names(pk.retrieved$masks). ​Or for
   #  my example code, try dba.show(atacDBA, mask=atacDBA$masks$dayPlus).
   #
   #  When calling dba.overlap() like this, the mask must identify exactly 2, 3 or 4 samples. Note you
   #  can call dba.plotVenn() and this will return the same peaksets as dba.overlap in addition to
   #  plotting their overlaps.
