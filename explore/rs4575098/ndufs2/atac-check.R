library(EndophenotypeExplorer)
library(ghdb)

targetGene <- "NDUFS2"
tag.snp <- "rs4575098"
ld.snp <-   "rs11585858" #  0.98 LD with 5098


etx <- EndophenotypeExplorer$new(targetGene, "hg19")


print(load("~/github/TrenaProjectAD/explore/mayo-epigenetics/atac/earlyResultMayoWholeGenome.RData"))
colnames(tbl.dba)[1] <- "chrom"
tbl.dba$chrom <- as.character(tbl.dba$chrom)
igv <- start.igv("ADAMTS4", "hg38")
zoomOut(igv)
roi <- getGenomicRegion(igv)
roi
tbl.track <- subset(tbl.dba, chrom==roi$chrom & start >= roi$start & end <= roi$end)
track <- DataFrameAnnotationTrack("atac", tbl.track, color="brown")
displayTrack(igv, track)


rsid.list <- c(tag.snp, ld.snp)
tbl.locs <- etx$rsidToLoc(rsid.list)

tbl.track <- tbl.locs[, c("chrom", "hg38", "hg38")]
colnames(tbl.track) <- c("chrom", "start", "end")
tbl.track$start <- tbl.track$start - 1
track <- DataFrameAnnotationTrack("snps", tbl.track, "red", trackHeight=24)
displayTrack(igv, track)


ghdb <- GeneHancerDB()
tbl.gh <- retrieveEnhancersFromDatabase(ghdb, "NDUFS2", tissues="all")
tbl.gh$score <- asinh(tbl.gh$combinedscore)
track <- DataFrameQuantitativeTrack("GH", tbl.gh[, c("chrom", "start", "end", "score")],
                                    autoscale=TRUE, color="brown")
displayTrack(igv, track)

tbl.gh <- retrieveEnhancersFromDatabase(ghdb, "ADAMTS4", tissues="all")
tbl.gh$score <- asinh(tbl.gh$combinedscore)
track <- DataFrameQuantitativeTrack("adamts4.GH", tbl.gh[, c("chrom", "start", "end", "score")],
                                    autoscale=TRUE, color="black")
displayTrack(igv, track)

