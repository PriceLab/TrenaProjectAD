library(MotifDb)
library(motifbreakR)
library(BSgenome.Hsapiens.UCSC.hg38)

loc <- 88884592
snp.name <- sprintf("%s:%d:%s:%s", "chr5", loc, "A", "C")
tbl.snp <- data.frame(chrom="chr5", start=loc-1, end=loc, name=snp.name, score=0, strand="+",
                      stringsAsFactors=FALSE)

bedFile <- "snp.bed"
write.table(tbl.snp, file=bedFile, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
getSeq(BSgenome.Hsapiens.UCSC.hg38, "chr5", loc, loc)
snps.gr <- snps.from.file(bedFile,
                          search.genome = BSgenome.Hsapiens.UCSC.hg38,
                          format = "bed")

motifs <- query(MotifDb, c("sapiens", "TBR1"))
motifs <- query(MotifDb, c("sapiens", "jaspar2018"), c("EGR3", "TBR1"))
results <- motifbreakR(snpList = snps.gr,
                          filterp = TRUE,
                          pwmList = motifs,
                          show.neutral=FALSE,
                          method = c("ic", "log", "notrans")[2],
                          bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                          BPPARAM = BiocParallel::bpparam(),
                          verbose=TRUE)

plotMB(results, rsid=snp.name, effect=c("weak", "strong"))

