library(trena)
library(RPostgreSQL)
library(ghdb)
#------------------------------------------------------------------------------------------------------------------------
tbl.38 <- get(load("~/github/TrenaProjectAD/inst/extdata/gwasLoci/tbl.posthuma-38-loci-curated.RData"))
tbl.assoc <- get(load("~/github/TrenaProjectAD/inst/extdata/gwasLoci/tbl.posthuma-38-geneAssociations-curated-3828x12.RData"))
tbl.atac <- get(load("~/github/TrenaProjectAD/explore/mayo-epigenetics/atac/dbaConsensusRegionsScored.74273x30.RData"))
tbl.eqtls.all <- get(load("tbl.eqtls.all38genes.RData"))

#----------------------------------------------------------------------------------------------------
get.eqtls <- function(gene, pval.threshold)
{
   subset(tbl.eqtls.all, genesymbol==gene & pvalue <= pval.threshold)

} # get.eqtls
#----------------------------------------------------------------------------------------------------
atac.filter.chrom.start.end <- function(tbl, tbl.name)
{
    stopifnot(all(c("chrom", "start", "end") %in% colnames(tbl)))
    gr <- GRanges(seqnames=tbl$chrom, IRanges(start=tbl$start, end=tbl$end))

    suppressWarnings(tbl.ov <- as.data.frame(findOverlaps(gr, gr.atac)))

    indices <- unique(tbl.ov[,1])
    printf("-- found %d entries in atac for %s", length(indices), tbl.name)

    return(tbl[indices,])

} # atac.filter.chrom.start.end
#----------------------------------------------------------------------------------------------------
atac.filter <- function(tbl, tbl.name)
{
    stopifnot(all(c("chrom", "hg38") %in% colnames(tbl)))

    gr <- GRanges(seqnames=tbl$chrom, IRanges(start=tbl$hg38, end=tbl$hg38))

    suppressWarnings(tbl.ov <- as.data.frame(findOverlaps(gr, gr.atac)))

    indices <- unique(tbl.ov[,1])
    printf("-- found %d entries in atac for %s", length(indices), tbl.name)

    return(tbl[indices,])

} # atac.filter
#----------------------------------------------------------------------------------------------------
gr.atac <- GRanges(seqnames=tbl.atac$chrom, IRanges(start=tbl.atac$start, end=tbl.atac$end))

    #------------------------------------------------------------
    # optional: fatten up the atac-peaks, realistically I hope
    #------------------------------------------------------------

atac.widths <- width(gr.atac)
head(atac.widths)
atac.padding <- as.integer(round(atac.widths/3))
start(gr.atac) <- (start(gr.atac) - atac.padding)
end(gr.atac) <- (end(gr.atac) + atac.padding)
atac.widths <- width(gr.atac)
head(atac.widths)

ghdb <- GeneHancerDB()

if(!exists("eqtl.db"))
    eqtl.db <- dbConnect(PostgreSQL(), user="trena", password="trena", dbname="genereg2021", host="khaleesi")

targetGene <- "BIN1"

igv <- start.igv(targetGene, "hg38")

tbl.gh <- retrieveEnhancersFromDatabase(ghdb, targetGene, tissues="all")
tbl.gh$chrom <- paste0("chr", tbl.gh$chrom)
dim(tbl.gh)
tbl.track <- tbl.gh[, c("chrom", "start", "end", "combinedscore")]
tbl.track$combinedscore <- asinh(tbl.track$combinedscore)
track <- DataFrameQuantitativeTrack("gh", tbl.track, color="grey", autoscale=FALSE, min=0, max=7.5)
displayTrack(igv, track)
shoulder <- 10000
showGenomicRegion(igv, sprintf("%s:%d-%d", tbl.gh$chrom[1], min(tbl.gh$start)-shoulder, max(tbl.gh$end)+shoulder))

roi <- getGenomicRegion(igv)
tbl.atac.sub <- subset(tbl.atac, chrom==roi$chrom & start >= roi$start & end <= roi$end)
dim(tbl.atac.sub)

track <- DataFrameAnnotationTrack("ATAC", tbl.atac.sub[, c("chrom", "start", "end")], color="red", trackHeight=25)
displayTrack(igv, track)

tbl.lead <- subset(tbl.38, geneSymbol==targetGene)[, c("chrom", "hg38", "hg38", "leadVariant")]
track <- DataFrameAnnotationTrack("leadSNP", tbl.lead, color="blue", trackHeight=25)
displayTrack(igv, track)

tbl.assoc.sub <- subset(tbl.assoc, gene==targetGene)[, c("chrom", "hg38", "hg38", "rsid")]
dim(tbl.assoc.sub)
track <- DataFrameAnnotationTrack("assoc", tbl.assoc.sub, color="blue", trackHeight=25)
displayTrack(igv, track)

tbl.eqtls <- get.eqtls(targetGene, 0.01)
dim(tbl.eqtls)
printf("------ eqtls found: %d", nrow(tbl.eqtls))

tbl.eqtls$score <- -log10(tbl.eqtls$pvalue)
tbl.sub <- tbl.eqtls
 #  tbl.sub <- subset(tbl.eqtls, tissue=="tcx", score >=0)
dim(tbl.sub)
track <- DataFrameQuantitativeTrack("eqtl scored", tbl.sub[, c("chrom", "hg38", "hg38", "score")], color="black", autoscale=TRUE, min=0, max=40, trackHeight=100)
displayTrack(igv, track)
track <- DataFrameAnnotationTrack("eqtl rsids", tbl.sub[, c("chrom", "hg38", "hg38", "rsid")], color="black", trackHeight=25)
displayTrack(igv, track)

tbl.sub.filtered <- atac.filter(tbl.sub, "eqtl")
track <- DataFrameQuantitativeTrack("eqtl scored filtered", tbl.sub.filtered[, c("chrom", "hg38", "hg38", "score")], color="black", autoscale=TRUE, min=0, max=40, trackHeight=100)
displayTrack(igv, track)

print(load("tbl.bin1.haplotype.snps.RData"))
track <- DataFrameAnnotationTrack("haploreg", tbl.hr[, c("chrom", "hg38", "hg38", "rsid")], color="random")
displayTrack(igv, track)

tbl.eqtl.oneAtacRegion <- subset(tbl.eqtls.all, chrom=="chr2" & hg38 > 127161565 & hg38 < 127162421)
tbl.eqtl.oneAtacRegion <- unique(tbl.eqtl.oneAtacRegion[, c("chrom", "hg38", "hg38", "rsid")])
dim(tbl.eqtl.oneAtacRegion)
track <- DataFrameAnnotationTrack("eQTL", tbl.eqtl.oneAtacRegion, color="random")
displayTrack(igv, track)

print(load("bin1.trena.fimo.230kroi.RData"))
#----------------------------------------------------------------------------------------------------
trena.and.fimo <- function()
{
    tbl.atacPlusGH <- reduce(GRanges(tbl.atac.sub), GRanges(tbl.gh), drop.empty.ranges=TRUE, min.gapwidth=10)

    tbl.roi <- data.frame(chrom="chr2", start=127015068, end=127245690, stringsAsFactors=FALSE)

    if(!file.exists("bin1.trena.fimo.230kroi.RData")){
        x <- runTrena(targetGene, tbl.roi, 1e-3)
        tbl.trena <- x$trena
        tbl.fimo <- x$fimo
        tbl.trena <- tbl.trena[order(abs(tbl.trena$spearmanCoeff), decreasing=TRUE),]
        tfbs.03 <- unlist(lapply(tbl.trena$gene, function(gene) nrow(subset(tbl.fimo, tf==gene & p.value <= 0.001))))
        tfbs.04 <- unlist(lapply(tbl.trena$gene, function(gene) nrow(subset(tbl.fimo, tf==gene & p.value <= 0.0001))))
        tfbs.05 <- unlist(lapply(tbl.trena$gene, function(gene) nrow(subset(tbl.fimo, tf==gene & p.value <= 0.00001))))
        tfbs.06 <- unlist(lapply(tbl.trena$gene, function(gene) nrow(subset(tbl.fimo, tf==gene & p.value <= 0.000001))))
        tbl.trena$tfbs.06 <- tfbs.06
        tbl.trena$tfbs.05 <- tfbs.05
        tbl.trena$tfbs.04 <- tfbs.04
        tbl.trena$tfbs.03 <- tfbs.03
        }

    tbl.fimo.sub <- subset(tbl.fimo, tf %in% tbl.trena$gene[1:10])
    dim(tbl.fimo.sub)

} # trena.and.fimo
#----------------------------------------------------------------------------------------------------
tbl.fimo.sub <- subset(tbl.fimo, p.value <= 0.001 & tf %in% tbl.trena$gene[1:10])
dim(tbl.fimo.sub)
as.data.frame(sort(table(tbl.fimo.sub$tf), decreasing=TRUE))

for(TF in head(tbl.trena$gene, n=10)){
    tbl.fimo.sub.tf <- subset(tbl.fimo.sub, tf==TF)[, c("chrom", "start", "end")]
    track <- DataFrameAnnotationTrack(TF, tbl.fimo.sub.tf, color="random", trackHeight=25)
    displayTrack(igv, track)
    } # for

list(igv=igv, trena=tbl.trena, fimo=tbl.fimo)

tbl.rs4439906 <- read.table("rs4439906.tsv", sep="\t", as.is=TRUE, header=FALSE)
colnames(tbl.rs4439906) <- c("chrom", "hg38", "rsq", "dprime", "rsid")
tbl.rs4439906$querySnp <- "rs4439906"
track <- DataFrameAnnotationTrack("rs4439906.kids", tbl.rs4439906[, c("chrom", "hg38", "hg38", "rsid")],
                                  color="random")
displayTrack(igv, track)

library(motifbreakR)
library(SNPlocs.Hsapiens.dbSNP151.GRCh38)
library(BSgenome.Hsapiens.UCSC.hg38)
#----------------------------------------------------------------------------------------------------
motifbreakrGOF <- function()
{
   rsids <- subset(tbl.assoc, gene=="BIN1")$rsid
   length(rsids) # 30
   snps.gr <- snps.from.rsid(rsids,
                             dbSNP=SNPlocs.Hsapiens.dbSNP151.GRCh38,
                             search.genome=BSgenome.Hsapiens.UCSC.hg38)

   motifs.selected <- query(MotifDb, "sapiens", c("jaspar2018", "hocomoco-core-A"))

   results <- motifbreakR(snpList = snps.gr,
                          filterp = TRUE,
                          pwmList = motifs.selected,
                          show.neutral=FALSE,
                          method = c("ic", "log", "notrans")[1],
                          bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                          BPPARAM = BiocParallel::bpparam(),
                          verbose=TRUE)

   save(results, file="27jun2021-motifbreakR-30eqtl-associated.RData")
   tbl.results <- as.data.frame(results, row.names=NULL)

   tbl.strong <- subset(tbl.results, effect=="strong" & pctAlt > 0.95)
   as.data.frame(sort(table(tbl.strong$SNP_id)))

} # motifbrearGOF
#----------------------------------------------------------------------------------------------------
diffBindSearch <- function()
{
   print(load("~/github/TrenaProjectAD/explore/mayo-epigenetics/atac/earlyResultMayoWholeGenome.RData"))
   tbl.dba$pvals.ad.psp <- pvals.ad.psp
   tbl.dba$pvals.ad.ctl <- pvals.ad.ctl
   tbl.dba$pvals.psp.ctl <- pvals.psp.ctl
   colnames(tbl.dba)[1] <- "chrom"
   tbl.dba$chrom <- as.character(tbl.dba$chrom)

   roi <- getGenomicRegion(igv)
   tbl.sub <- subset(tbl.dba, chrom==roi$chrom & start > roi$start & end < roi$end)
   dim(tbl.sub)   # 23 33

} # diffBindSearch
#----------------------------------------------------------------------------------------------------
displayPeak <- function()
{
   roi <- getGenomicRegion(igv)
   tbl.sub <- subset(tbl.dba, chrom==roi$chrom & start > roi$start & end < roi$end)
   dim(tbl.sub)

   for(sample in sort(colnames(tbl.sub)[6:27])){
       trackName <- sample
       if(grepl("^ad", trackName)) color <- "lightblue"
       if(grepl("^psp", trackName)) color <- "red"
       if(grepl("^ctl", trackName)) color <- "green"
       track <- DataFrameQuantitativeTrack(trackName,
                                           tbl.sub[, c("chrom", "start", "end", sample)],
                                           autoscale=FALSE, color=color,
                                           trackHeight=25, min=0, max=0.01)
       displayTrack(igv, track)
       } # for file

} # displayPeak
#------------------------------------------------------------------------------------------------------------------------
displayTFBS <- function()
{
  print(load("bin1-gh+atac-roi-trenaAndFimo.RData"))
  tbl.trena <- x$trena
  tbl.fimo  <- x$fimo

  tfbs.03 <- unlist(lapply(tbl.trena$gene, function(gene) nrow(subset(tbl.fimo, tf==gene & p.value <= 0.001))))
  tfbs.04 <- unlist(lapply(tbl.trena$gene, function(gene) nrow(subset(tbl.fimo, tf==gene & p.value <= 0.0001))))
  tfbs.05 <- unlist(lapply(tbl.trena$gene, function(gene) nrow(subset(tbl.fimo, tf==gene & p.value <= 0.00001))))
  tfbs.06 <- unlist(lapply(tbl.trena$gene, function(gene) nrow(subset(tbl.fimo, tf==gene & p.value <= 0.000001))))
  tbl.trena$tfbs.06 <- tfbs.06
  tbl.trena$tfbs.05 <- tfbs.05
  tbl.trena$tfbs.04 <- tfbs.04
  tbl.trena$tfbs.03 <- tfbs.03

  head(tbl.trena, n=10)
  tfs <- head(tbl.trena$gene, n=10)
  threshold <- 1e-4
  for(TF in tfs){
      tbl.track <- subset(tbl.fimo, tf==TF & p.value <= threshold)[, c("chrom", "start", "end", "tf")]
      tbl.track <- atac.filter.chrom.start.end(tbl.track, "fimo")
      track <- DataFrameAnnotationTrack(TF, tbl.track, color="black", trackHeight=20)
      displayTrack(igv, track)
      }


} # displayTFBS
#------------------------------------------------------------------------------------------------------------------------
tfeb <- function()
{
  print(load("bin1-gh+atac-roi-trenaAndFimo.RData"))
  tbl.trena <- x$trena
  tbl.trena <- tbl.trena[order(tbl.trena$pearsonCoeff, decreasing=TRUE),]
  rownames(tbl.trena) <- NULL
  head(tbl.trena)
  tbl.fimo  <- x$fimo
  tbl.fimo.tfeb <- subset(tbl.fimo, tf=="TFEB" & p.value < 1e-4)
  track <- DataFrameAnnotationTrack("tfeb.5", tbl.fimo.tfeb[, c("chrom", "start", "end")], color="red", trackHeight=25)
  displayTrack(igv, track)

} # tfeb
#------------------------------------------------------------------------------------------------------------------------
