library(EndophenotypeExplorer)
targetGene <- "LDLR"
etx <- EndophenotypeExplorer$new(targetGene, "hg38")
tbl.eqtls <- etx$getEQTLsForGene()
dim(tbl.eqtls)
tbl.eqtls.sub <- subset(tbl.eqtls, study=="ampad-mayo" & tissue == "tcx" & pvalue <= 0.05)
dim(tbl.eqtls.sub)  # 94 10

tbl.eqtls.sub <- tbl.eqtls.sub[, c("chrom", "hg38", "hg38", "rsid", "pvalue")]
colnames(tbl.eqtls.sub) <- c("chrom", "start", "end", "rsid", "pvalue")
shoulder <- 10
tbl.eqtls.sub$start <- tbl.eqtls.sub$start - shoulder
tbl.eqtls.sub$end <- tbl.eqtls.sub$end + shoulder

loc.chrom <- tbl.eqtls.sub$chrom[1]
loc.start <- min(tbl.eqtls.sub$start)
loc.end   <- max(tbl.eqtls.sub$end)
loc.end - loc.start   # about 2M

f <- "~/github/TrenaProjectAD/explore/mayo-epigenetics/atac/mayoAllPeaks.merged.96064x4.RData"
tbl.atac <- get(load(f))
dim(tbl.atac)  # 96064  5
tbl.atac <- subset(tbl.atac, chrom==loc.chrom & start >= loc.start & end <= loc.end)
dim(tbl.atac)  # 113 5

igv <- start.igv(targetGene, "hg38")
showGenomicRegion(igv, targetGene)
zoomOut(igv); zoomOut(igv)
track <- DataFrameAnnotationTrack("atac", tbl.atac, color="random")
displayTrack(igv, track)

track <- DataFrameAnnotationTrack("eqtl", tbl.eqtls.sub, color="random")
displayTrack(igv, track)

gr.eqtls <- GRanges(tbl.eqtls.sub)
gr.atac  <- GRanges(tbl.atac)

gr.merged <- reduce(c(gr.eqtls, gr.atac))
gr.merged   # 230

tbl.merged <- as.data.frame(gr.merged)
colnames(tbl.merged)[1] <- "chrom"
tbl.merged$chrom <- as.character(tbl.merged$chrom)

track <- DataFrameAnnotationTrack("merged", tbl.merged, color="random")
displayTrack(igv, track)

sum(width(gr.eqtls))  # 1974
sum(width(gr.atac))   # 79106
sum(width(union(gr.eqtls, gr.atac)))  # 81957

tbl.oc <- tbl.merged[, c("chrom", "start", "end")]
head(tbl.oc)
save(tbl.oc, file="tbl.oc.ldlr.atac.eqtl.merged.RData")

