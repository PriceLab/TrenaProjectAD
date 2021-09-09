library(EndophenotypeExplorer)
library(ghdb)
targetGene <- "PTK2B"
ghdb <- GeneHancerDB()
tbl.gh <- retrieveEnhancersFromDatabase(ghdb, targetGene, tissues="all")

etx <- EndophenotypeExplorer$new("PTK2B", "hg38")
igv <- start.igv("PTK2B", "hg38")
zoomOut(igv);zoomOut(igv);
track <- DataFrameQuantitativeTrack("GH", tbl.gh[, c("chrom", "start", "end", "combinedscore")],
                                       autoscale=TRUE, color="brown")
displayTrack(igv, track)


dir <- "~/github/TrenaProjectAD/inst/extdata/gwasLoci"
tbl.wms <- get(load(file.path(dir, "williams-natureNeuroscience2020.RData")))
tbl.pos <- get(load(file.path(dir, "tbl.posthuma-38-loci-curated.RData")))
tbl.posAssoc <- get(load(file.path(dir, "tbl.posthuma-38-loci-curated.RData")))
subset(tbl.wms, locusOrGene=="PTK2B")

#    locusOrGene       rsid
# 24       PTK2B rs28834970

tbl.eQTL <- etx$getEQTLsForGene()
dim(tbl.eQTL)   # 16693 10
tbl.eQTL.sig <- subset(tbl.eQTL, pvalue < 0.5)
dim(tbl.eQTL.sig) # 17


tag.snp <- "rs28834970"
subset(tbl.eQTL, rsid==tag.snp)

#       chrom     hg19     hg38       rsid    pvalue            ensg genesymbol        study tissue   assay
# 1617   chr8 27195121 27337604 rs28834970 0.7766134 ENSG00000120899      PTK2B   ampad-mayo    tcx unknown
# 7207   chr8 27195121 27337604 rs28834970 0.1501057 ENSG00000120899      PTK2B   ampad-mayo    cer unknown
# 12792  chr8 27195121 27337604 rs28834970 0.5057638 ENSG00000120899      PTK2B ampad-rosmap  dlpfc unknown

track <- DataFrameAnnotationTrack("tag snp, rs28834970",
                                  data.frame(chrom="chr8", start=27337603, end=27337604, stringsAsFactors=FALSE),
                                  color="red", trackHeight=24)
displayTrack(igv, track)

tbl.hap <- read.table("haplotype-rs29934970.tsv", sep="\t", header=TRUE, as.is=TRUE)
tbl.hap$end <- tbl.hap$start
tbl.hap$start <- tbl.hap$start - 1
tbl.track <- tbl.hap[, c("chrom", "start", "end", "rSquared", "rsid")]

track <- DataFrameQuantitativeTrack("haploreg", tbl.track, color="brown", autoscale=TRUE)
displayTrack(igv, track)

tbl.track <- unique(tbl.eQTL.sig[, c("chrom", "hg38", "hg38", "pvalue")])
colnames(tbl.track)[c(2,3)] <- c("start", "end")
tbl.track$score <- -log10(tbl.track$pvalue)
tbl.track <- tbl.track[, c("chrom", "start", "end", "score")]
track <- DataFrameQuantitativeTrack("eqtl", tbl.track, autoscale=TRUE, color="darkBlue")
displayTrack(igv, track)


f <- "~/github/TrenaProjectAD/explore/mayo-epigenetics/atac/dbaConsensusRegionsScored.74273x30.RData"
tbl.atac <- get(load(f))
roi <- getGenomicRegion(igv)

tbl.atac.sub <- subset(tbl.atac, chrom==roi$chrom[1] & start >= roi$start & end <= roi$end)
tbl.track <- tbl.atac.sub[, c("chrom", "start", "end", "max.score")]
track <- DataFrameQuantitativeTrack("atac", tbl.track, autoscale=TRUE, color="darkgreen")
displayTrack(igv, track)



f <- "~/github/TrenaProjectAD/explore/mayo-epigenetics/atac/mayoAllPeaks.1052789x4.RData"
tbl.atac <- get(load(f))
roi <- getGenomicRegion(igv)

tbl.atac.sub <- subset(tbl.atac, chrom==roi$chrom[1] & start >= roi$start & end <= roi$end)
tbl.track <- tbl.atac.sub[, c("chrom", "start", "end", "score")]
track <- DataFrameQuantitativeTrack("atac-2", tbl.track, autoscale=TRUE, color="darkgreen")
displayTrack(igv, track)

library(RPostgreSQL)

db <- dbConnect(PostgreSQL(), user= "trena", password="trena", dbname="brain_hint_16", host="khaleesi")
          # browser()
query <- sprintf("select * from regions where chrom='%s' and start >= %d and endpos <= %d",
                 roi$chrom, roi$start, roi$end)

tbl.fp <- dbGetQuery(db, query)
dim(tbl.fp)
dbDisconnect(db)
tbl.fp <- tbl.fp[, c("chrom", "start", "endpos")]
colnames(tbl.fp) <- c("chrom", "start", "end")
gr.fp <- GRanges(tbl.fp)

gr.fp.reduced <- reduce(gr.fp)
tbl.fp.reduced <- as.data.frame(gr.fp.reduced)
colnames(tbl.fp.reduced)[1] <- "chrom"
tbl.fp.reduced$chrom <- as.character(tbl.fp.reduced$chrom)

track <- DataFrameAnnotationTrack("fp", tbl.fp.reduced, color="black")
displayTrack(igv, track)
