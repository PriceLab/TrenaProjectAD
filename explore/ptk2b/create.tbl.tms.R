library(TrenaMultiScore)
library(TrenaProjectAD)
library(trena)
library(EndophenotypeExplorer)
targetGene <- "PTK2B"
#----------------------------------------------------------------------------------------------------
get.precalculated.fp.fimo <- function()
{
   invisible(get(load("tbl.fimo.PTK2B-footprints.RData")))
}
#----------------------------------------------------------------------------------------------------
get.rna.matrices <- function()
{
    dir <- "~/github/TrenaProjectAD/prep/rna-seq-counts-from-synapse/eqtl"
    files <- c(mayo="mtx.mayo.rnaseq-residual-eqtl-geneSymbols-patients-17009x523.RData",     # 68826146 Aug 7
               sinai="mtx.mssm.rnaseq-residual-eqtl-geneSymbols-patients-16346x753.RData",    # 93294924 Aug 7
               rosmap="mtx.rosmap.rnaseq-residual-eqtl-geneSymbols-patients-15582x632.RData") # 75390699 Aug 7

    mtx.rosmap <- get(load(file.path(dir, files[["rosmap"]])))
    mtx.mayo <- get(load(file.path(dir, files[["mayo"]])))
    mtx.sinai <- get(load(file.path(dir, files[["sinai"]])))

    mtx.rosmap[is.na(mtx.rosmap)] <- 0
    mtx.mayo[is.na(mtx.mayo)] <- 0
    mtx.sinai[is.na(mtx.sinai)] <- 0

    goi <- intersect(rownames(mtx.rosmap), intersect(rownames(mtx.mayo), rownames(mtx.sinai)))

    mtx.rosmap <- mtx.rosmap[goi,]
    mtx.mayo <- mtx.mayo[goi,]
    mtx.sinai <- mtx.sinai[goi,]

    mtx.rna <- cbind(mtx.rosmap, mtx.mayo, mtx.sinai)
    list(all=mtx.rna, rosmap=mtx.rosmap, mayo=mtx.mayo, sinai=mtx.sinai)

} # get.rna.matrices
#----------------------------------------------------------------------------------------------------
get.atac.oc <- function()
{
   get(load("~/github/TrenaProjectAD/explore/mayo-epigenetics/atac/mayoAllPeaks.1052789x4.RData"))

} # get.atac.oc
#----------------------------------------------------------------------------------------------------
create.tbl.tms <- function()
{
   tpad <- TrenaProjectAD()
   tms <- TrenaMultiScore(tpad, targetGene, tbl.fimo=get.precalculated.fp.fimo(), tbl.oc=get.atac.oc(), quiet=FALSE)

   suppressWarnings(
     db.access.test <- try(system("/sbin/ping -c 1 khaleesi", intern=TRUE, ignore.stderr=TRUE)))
   if(length(db.access.test) == 0)
     stop("khaleesi database server unavailable")

   x <- get.rna.matrices()
   addGeneExpressionCorrelations(tms, x$all, "cor.all")
   addGeneExpressionCorrelations(tms, x$rosmap, "cor.rosmap")
   addGeneExpressionCorrelations(tms, x$sinai, "cor.sinai")
   addGeneExpressionCorrelations(tms, x$mayo, "cor.mayo")

   scoreMotifHitsForGeneHancer(tms)
   scoreMotifHitsForConservation(tms)
   scoreMotifHitsForOpenChromatin(tms)
   addDistanceToTSS(tms)
   addGenicAnnotations(tms)
   addChIP(tms)
   tbl.tms <- getMultiScoreTable(tms)
   dim(tbl.tms)
   save(tbl.tms, x, file="tbl.tms.RData")

} # create.tbl.tms
#----------------------------------------------------------------------------------------------------


