library(DiffBind)

dba <- dba(sampleSheet="mayoSamples-all.csv")

# dba <- dba.normalize(dba)

overlap <- dba.overlap(dba, mode=DBA_OLAP_RATE)
# 96127 74273 66250 61089 57192 53997 51112 48609 46288 44139 42034 40053 38066 36211 34250 32296 30178 27872 25387 22455 18317 12274
# peaks appearing in at least 1, 2, 3 peaksets: 68760 42452 30772
plot(overlap, type="b", ylab="# peaks", xlab="Overlap at least this many peaksets")

names(dba$masks)

#  [1] "CER"         "XYZ"         "PSP"         ""            "bed"
#  [6] "Replicate.1" "Replicate.2" "Replicate.3" "All"         "None"

dba.overlap(dba, dba$masks$PSP, mode=DBA_OLAP_RATE)
# 86897 66790 56942 49674 43774 38399 33291 27389 13954

dba.overlap(dba, dba$masks$AD, mode=DBA_OLAP_RATE)
# 77508 59629 51578 45344 39758 33895 26718 18640

dba.overlap(dba, dba$masks$CTL, mode=DBA_OLAP_RATE)
# 67292 47480 38969 32003 24527

dba.plotVenn(dba, dba$masks$CER)

dba.consensus <- dba.peakset(dba, consensus=c(DBA_CONDITION))
dba.consensus
# 25 Samples, 74273 sites in matrix (96127 total):
#        ID Tissue Factor Condition Treatment         Replicate Intervals
# 23     AD    CER    XYZ        AD      None   1-2-3-4-5-6-7-8     59629
# 24    PSP    CER    XYZ       PSP      None 1-2-3-4-5-6-7-8-9     66790
# 25    CTL    CER    XYZ       CTL      None         1-2-3-4-5     47480

dba <- dba.peakset(dba.consensus, bRetrieve=TRUE)
length(dba) # 74273
tbl.dba <- as.data.frame(dba)
dim(tbl.dba)
max.score <- as.numeric(apply(as.matrix(tbl.dba[, c("AD", "PSP", "CTL")]), 1, max))
tbl.dba$max.score <- max.score
save(tbl.dba, file="dbaConsensusRegionsScored.74273x30.RData")
mtx.dba <- as.matrix(tbl.dba[, 6:27])
ad.cols <- 1:8
psp.cols <- 9:17
ctl.cols <- 18:22

peak.names <- sprintf("%s:%d-%d", tbl.dba$seqnames, tbl.dba$start, tbl.dba$end)
tbl.dba$peakName <- peak.names
length(peak.names)  # 74273

pvals.ad.psp <- apply(mtx.dba, 1, function(vec) t.test(vec[ad.cols], vec[psp.cols])$p.value)
names(pvals.ad.psp) <- peak.names
peaks.ad.psp.disagree <- names(pvals.ad.psp[pvals.ad.psp<0.05])
length(peaks.ad.psp.disagree) # 579

pvals.ad.ctl <- apply(mtx.dba, 1, function(vec) t.test(vec[ad.cols], vec[ctl.cols])$p.value)
names(pvals.ad.ctl) <- peak.names
peaks.ad.ctl.disagree <- pvals.ad.ctl[pvals.ad.ctl<0.05] # 2969
length(peaks.ad.ctl.disagree) # 4758

pvals.psp.ctl <- apply(mtx.dba, 1, function(vec) t.test(vec[psp.cols], vec[ctl.cols])$p.value)
names(pvals.psp.ctl) <- peak.names
peaks.psp.ctl.disagree <- pvals.psp.ctl[pvals.psp.ctl<0.05]
length(peaks.psp.ctl.disagree)  # 4764

peaks.both.disagree <- intersect(names(peaks.ad.ctl.disagree), names(peaks.psp.ctl.disagree))
length(peaks.both.disagree)  # 779

save(tbl.dba, mtx.dba, pvals.ad.psp, pvals.ad.ctl, pvals.psp.ctl, file="earlyResultMayoWholeGenome.RData")

# xue requests:
#   overlap ADvCTL and PSBvCTL peaks

peaks.psp.ctl <- names(which(pvals.psp.ctl < 0.01))
peaks.ad.ctl <- names(which(pvals.ad.ctl < 0.01))
peaks.psp.ad.ctl <- intersect(peaks.psp.ctl, peaks.ad.ctl)
length(peaks.psp.ad.ctl)

tbl.dba[unique(which(pvals < 0.05)),]
plot(tbl.dba$PSP, tbl.dba$CTL)

# GRanges object with 51257 ranges and 6 metadata columns:
#         seqnames            ranges strand |     psp.04     psp.11     ctl.20     ctl.46        PSP        CTL
#            <Rle>         <IRanges>  <Rle> |  <numeric>  <numeric>  <numeric>  <numeric>  <numeric>  <numeric>
#       1     chr1       10055-10392      * | 0.00371839 0.00328262 0.00652786 0.01014233 0.00350051 0.00833509
#       2     chr1     180748-180978      * | 0.00371839 0.00000000 0.00187569 0.00318746 0.00000000 0.00253158
#       3     chr1     181245-181629      * | 0.00405665 0.00000000 0.00163788 0.00230434 0.00000000 0.00197111
#       4     chr1     191153-191896      * | 0.03455901 0.01660308 0.01823152 0.01507159 0.02558105 0.01665156
#       5     chr1     629657-630280      * | 0.43255130 0.53244037 0.60523522 0.42531368 0.48249584 0.51527445
#     ...      ...               ...    ... .        ...        ...        ...        ...        ...        ...
#   51253     chrY 56838878-56839164      * | 0.00410688 0.01010500 0.00640315 0.00478931 0.00710594 0.00559623
#   51254     chrY 56847057-56847313      * | 0.00765874 0.00414814 0.00494741 0.00429933 0.00590344 0.00462337
#   51255     chrY 56850312-56850558      * | 0.00390596 0.00367987 0.00000000 0.00000000 0.00379292 0.00000000
#   51256     chrY 56850765-56851018      * | 0.00348970 0.00257160 0.00401717 0.00207868 0.00303065 0.00304792
#   51257     chrY 56870850-56871079      * | 0.00278750 0.00166182 0.00000000 0.00000000 0.00222466 0.00000000
#   -------
#   seqinfo: 24 sequences from an unspecified genome; no seqlengths

#x <- dba.count(dba, peaks=dba$masks$Consensus)
#dba <- dba.count(dba, peaks=NULL,minCount=1,filter=5)

tbl.psp.11 <- get.table(psp.files[2])
rownames(tbl.psp.11) <- NULL
dba <- dba.peakset(dba, peaks=GRanges(tbl.psp.11[, 1:3]),
                   sampID="psp.11", consensus=FALSE,
                   factor="XYZ", tissue="CER",
                   condition="PSP", replicate=2, counts=tbl.psp.11[, 5])
dba$masks
dba$config$doBlacklist<- FALSE
dba$doGreylist <- FALSE



tbl.psp <-

file <- "~/github/TrenaProjectErythropoiesis/inst/extdata/genomicRegions/tbl.corces.hg38.day0thru4.RData"
file.exists(file)
print(load(file))
dim(tbl.corces)  # 592665 10
tbl <- subset(tbl.corces, chrom=="chr19")

consensus <- GRanges(tbl[,1:3])                          # all regions
counts    <- as.matrix(tbl[,4:ncol(tbl)])    # 592665      7
sampIDs   <- colnames(counts)                                        # "day0.1" "day0.2" "day0.3" "day0.4" "day2.1" "day2.2" "day4.1"
ids       <- strsplit(sampIDs,".",fixed=TRUE)
condition <- unlist(lapply(ids,function(x)x[1]))                     # "day0" "day0" "day0" "day0" "day2" "day2" "day4"
treatment <- condition
treatment[condition != "day0"] <- "dayPlus"                          # "day0"    "day0"    "day0"    "day0"    "dayPlus" "dayPlus" "dayPlus"
replicate <- as.numeric(unlist(lapply(ids,function(x)x[2])))         # 1 2 3 4 1 2 1

atacDBA <- NULL
for(i in 1:ncol(counts)) {
  cat(i,sampIDs[i],"\n")
  atacDBA <- dba.peakset(atacDBA, sampID=sampIDs[i], peaks=consensus,
                         factor="ATAC",condition=condition[i],
                         treatment=treatment[i],replicate=replicate[i],
                         counts=counts[,i])
  } # for

atacDBA$config$doBlacklist <- atacDBA$config$doGreylist <- FALSE
atacDBA
atacDBA <- dba.count(atacDBA,peaks=NULL,minCount=1,filter=5)
atacDBA
dba.plotPCA(atacDBA, DBA_CONDITION, label=DBA_REPLICATE)
atacDBA <- dba.contrast(atacDBA)
atacDBA <- dba.analyze(atacDBA)
atacDBA$config$th <- 0.01
dba.show(atacDBA, bContrasts = TRUE)
dba.save(atacDBA, file="shannon_atac")
dba.plotMA(atacDBA)
dba.plotVolcano(atacDBA)
dba.plotPCA(atacDBA, contrast=1, attributes=DBA_TREATMENT, label=DBA_CONDITION)
chrom.diff <- dba.report(atacDBA)  # 13483
tbl.rory <- as.data.frame(chrom.diff)
