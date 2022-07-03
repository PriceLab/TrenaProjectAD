#------------------------------------------------------------------------------------------------------------------------
library(EndophenotypeExplorer)
library(TrenaProjectAD)
library(RPostgreSQL)
library(plyr)
library(ghdb)
#source("~/github/endophenotypeExplorer/R/getExpressionMatrices.R")
source("~/github/TrenaMultiScore/tools/runner/v2/tmsCore.R")
library(RUnit)
state <- new.env(parent=emptyenv())
load("~/github/TrenaProjectAD/explore/enhancers/nott/tbl.nott-microglial-chromatin-interactions.RData")
#------------------------------------------------------------------------------------------------------------------------
targetGene <- "NDUFS2"
if(!exists("etx")){
   etx <- EndophenotypeExplorer$new(targetGene, "hg38", initialize.snpLocs=TRUE)
   ghdb <- GeneHancerDB()
   }

tag.snp <- "rs4575098"
tag.snp.loc <- 161185602
ld.snp  <- "rs11585858"

if(!exists("igv")){
   igv <- start.igv(targetGene, "hg38")
   zoomOut(igv)
   zoomOut(igv)
   }
#------------------------------------------------------------------------------------------------------------------------
# redraw igv from a fresh start
reviz <- function()
{
  display.snps.and.LD.partners()
  display.eqtls(c("NDUFS2"))
  # display.nott.interactomes

} # reviz
#------------------------------------------------------------------------------------------------------------------------
combine.trena.results.from.three.mrna.matrices <- function()
{
   mtx.old  <- get.rna.matrix("old-rosmap")
   mtx.max <- get.rna.matrix("max-rosmap")
   mtx.eqtl <- get.rna.matrix("sage-eqtl-rosmap")
   mtx.eqtl[is.na(mtx.eqtl)] <- 0

   tbl.old <- tms$build.trena.model(tfs, mtx.old)
   tbl.max <- tms$build.trena.model(tfs, mtx.max)
   tbl.eqtl <- tms$build.trena.model(tfs, mtx.eqtl)

   tbl.old$targetGene <- targetGene
   tbl.max$targetGene <- targetGene
   tbl.eqtl$targetGene <- targetGene

   new.order.old <- order(abs(tbl.old$spearmanCoeff), decreasing=TRUE)
   new.order.max <- order(abs(tbl.max$spearmanCoeff), decreasing=TRUE)
   new.order.eqtl <- order(abs(tbl.eqtl$spearmanCoeff), decreasing=TRUE)

   tbl.old <- tbl.old[new.order.old,]
   tbl.max <- tbl.max[new.order.max,]
   tbl.eqtl <- tbl.eqtl[new.order.eqtl,]

   rownames(tbl.old) <- NULL
   rownames(tbl.max) <- NULL
   rownames(tbl.eqtl) <- NULL

   tbl.old.pretty <- head(tbl.old, n=20)[, -(c(8,9))]
   for(column.name in colnames(tbl.old.pretty)[-c(1,8)])
       tbl.old.pretty[[column.name]] <- round(tbl.old.pretty[[column.name]], digits=3)

   tbl.max.pretty <- head(tbl.max, n=20)[, -(c(8,9))]
   for(column.name in colnames(tbl.max.pretty)[-c(1,8)])
       tbl.max.pretty[[column.name]] <- round(tbl.max.pretty[[column.name]], digits=3)

   tbl.eqtl.pretty <- head(tbl.eqtl, n=20)[, -(c(8,9))]
   for(column.name in colnames(tbl.eqtl.pretty)[-c(1,8)])
       tbl.eqtl.pretty[[column.name]] <- round(tbl.eqtl.pretty[[column.name]], digits=3)


   head(tbl.old, n=20)
   head(tbl.max, n=20)
   head(tbl.eqtl, n=20)

   rownames(tbl.trena) <- NULL


} # combine.trena.results.from.three.mrna.matrices
#------------------------------------------------------------------------------------------------------------------------
display.snps.and.LD.partners <- function()
{
      #------------------------------------------------------------
      # first: a track with rsids, visible when clicked upon
      #------------------------------------------------------------

   tbl.snp <- data.frame(chrom="chr1", start=tag.snp.loc-1, end=tag.snp.loc,
                         stringsAsFactors=FALSE)
   track <- DataFrameAnnotationTrack("rs4575098", tbl.snp, color="black")
   displayTrack(igv, track)

   tbl.hap <- read.table("../haploreg-rs4575098-0.2.tsv", sep="\t", as.is=TRUE, header=TRUE, nrow=-1)
   tbl.track <- tbl.hap[, c("chrom", "hg38", "hg38", "rsid")]
   colnames(tbl.track)[c(2,3)] <- c("start", "end")
   tbl.track$start <- tbl.track$start - 1
   track <- DataFrameAnnotationTrack("hapName", tbl.track, color="brown")
   #displayTrack(igv, track)


   shoulder <- 5000
   #roi <- sprintf("%s:%d-%d", tbl.track$chrom[1], min(tbl.track$start)-shoulder, max(tbl.track$end)+shoulder)
   #showGenomicRegion(igv, roi)

      #------------------------------------------------------------
      # now a track with rSquared magnitudes
      #------------------------------------------------------------

   tbl.track <- tbl.hap[, c("chrom", "hg38", "hg38", "rSquared")]
   colnames(tbl.track)[c(2,3)] <- c("start", "end")
   tbl.track$start <- tbl.track$start - 1
   track <- DataFrameQuantitativeTrack("eur LD", tbl.track, color="darkRed", autoscale=FALSE, min=0, max=1)
   displayTrack(igv, track)

} # display.snps.and.LD.partners
#------------------------------------------------------------------------------------------------------------------------
display.mayo.cer.ndufs2.eqtls <- function()
{
    tbl.eqtl.mayo <- subset(state$tbl.eqtl.all, study=="ampad-mayo" & tissue=="cer" &
                            genesymbol=="NDUFS2" & pvalue < 0.05)
    dim(tbl.eqtl.mayo)
    tbl.track <- tbl.eqtl.mayo[, c("chrom", "hg38", "hg38", "pvalue")]
    colnames(tbl.track) <- c("chrom", "start", "end", "score")
    tbl.track$start <- tbl.track$start - 1
    tbl.track$score <- -log10(tbl.track$score)
    track <- DataFrameQuantitativeTrack(sprintf("%s-mayo-cer--eQTL", targetGene),
                                        tbl.track, autoscale=FALSE,
                                        min=0, max=20, color="random")
    displayTrack(igv, track)



} # display.mayo.cer.ndufs2.eqtls
#------------------------------------------------------------------------------------------------------------------------
display.eqtls <- function(targets=NA)
{
    db <- dbConnect(PostgreSQL(), user="trena", password="trena", dbname="genereg2021", host="khaleesi")
    dbGetQuery(db, "select * from eqtls limit 3")
    getGenomicRegion(igv)

    roi <- getGenomicRegion(igv)

    shoulder <- 1000000
    shoulder <- 500000
    roi.string <- sprintf("chr1:%d-%d", tag.snp.loc - shoulder, tag.snp.loc + shoulder)
    showGenomicRegion(igv, roi.string)
    roi <- getGenomicRegion(igv)
    query <- with(roi, sprintf("select * from eqtls where chrom='%s' and hg38 > %d and hg38 < %d",
                               chrom, start, end))
    if(!("tbl.eqtl.all" %in% names(state))){
       tbl.eqtl.all <- dbGetQuery(db, query)
       state$tbl.eqtl.all <- tbl.eqtl.all
       dim(tbl.eqtl.all)
       }
    tbl.eqtl.sub <- subset(state$tbl.eqtl.all, study=="ampad-rosmap")
    tbl.eqtl.sub <- subset(tbl.eqtl.sub, pvalue < 0.05)
    dim(tbl.eqtl.sub)

    goi <- sort(unique(subset(tbl.eqtl.sub, abs(hg38-tag.snp.loc) < 1000 & pvalue < 0.05)$genesymbol))
    length(goi)  # 13
    removeTracksByName(igv, grep("-eQTL", getTrackNames(igv), v=TRUE))
    deleters <- which(nchar(goi) == 0)
    if(length(deleters) > 0)
        goi <- goi[-deleters]
    hand.chosen <- c("DUSP12", "FCER1G", "NDUFS2", "NOS1AP", "PEA15", "PPOX", "USF1")
    goi <- hand.chosen
    if(!is.na(targets))
        goi <- targets
    set.seed(17)
    for(gene in goi){ # "NDUFS2"){
       tbl.sub <- subset(tbl.eqtl.sub, genesymbol==gene)
       tbl.track <- tbl.sub[, c("chrom", "hg38", "hg38", "pvalue")]
       colnames(tbl.track) <- c("chrom", "start", "end", "score")
       tbl.track$start <- tbl.track$start - 1
       tbl.track$score <- -log10(tbl.track$score)
       track <- DataFrameQuantitativeTrack(sprintf("%s-eQTL", gene),
                                           tbl.track, autoscale=FALSE,
                                           min=0, max=20, color="random")
       displayTrack(igv, track)
       }

} # display.eqtls
#------------------------------------------------------------------------------------------------------------------------
human.brain.mayo.oc <- function()
{
   dir <- "~/github/TrenaProjectAD/inst/extdata/genomicRegions"
   file <- "mayoAllPeaks.merged.96064x4.RData"
   full.path <- file.path(dir, file)
   tbl.atac <- get(load(full.path))
   tbl.oc <- tbl.atac[, c("chrom", "start", "end")]

   invisible(tbl.oc)

} # human.brain.mayo.oc
#---------------------------------------------------------------------------------------------------
human.brain.consensus.boca.oc <- function()
{
    dir <- "~/github/TrenaProjectAD/inst/extdata/genomicRegions"
    f <- "boca-hg38-consensus-ATAC.RData"
    full.path <- file.path(dir, f)
    stopifnot(file.exists(full.path))
    tbl.oc <- get(load(full.path))
    dim(tbl.oc)   # 125430

    invisible(tbl.oc)

} # human.brain.consensus.boca.oc
#---------------------------------------------------------------------------------------------------
display.gh <- function()
{
    tbl.gh <- retrieveEnhancersFromDatabase(ghdb, targetGene, tissues="all")
    tbl.gh$score <- tbl.gh$combinedscore
    track <- DataFrameQuantitativeTrack("GH", tbl.gh[, c("chrom", "start", "end", "score")],
                                    autoscale=TRUE, color="brown")
    displayTrack(igv, track)
    state$tbl.gh <- tbl.gh

} # display.gh
#---------------------------------------------------------------------------------------------------
display.openChromatin <- function()
{
   tbl.atac <- human.brain.mayo.oc()
   track <- DataFrameAnnotationTrack("MayoATAC", tbl.atac, color="darkgreen")
   displayTrack(igv, track)

   tbl.atac <- human.brain.consensus.boca.oc()
   track <- DataFrameAnnotationTrack("bocaATAC", tbl.atac, color="darkblue")
   displayTrack(igv, track)

   data.dir <- "~/github/TrenaProjectAD/explore/mayo-epigenetics/atac"
   filename <- "mayoAllPeaks.1052789x4.RData"
   full.path <- file.path(data.dir, filename)
   file.exists(full.path)
   tbl.atac <- get(load(full.path))
   head(tbl.atac)
   dim(tbl.atac) # 1052789 4
   roi <- getGenomicRegion(igv)
   tbl.atac.sub <- subset(tbl.atac, chrom==roi$chrom & start >= roi$start & end <= roi$end)
   dim(tbl.atac.sub)  # 159 5

   tbl.atac.ad <- subset(tbl.atac, dx=="AD")
   tbl.atac.psp <- subset(tbl.atac, dx=="PSP")
   tbl.atac.ctl <- subset(tbl.atac, dx=="CTL")

   track <- DataFrameAnnotationTrack("Mayo-AD", tbl.atac.ad, color="red")
   displayTrack(igv, track)

   track <- DataFrameAnnotationTrack("Mayo-PSP", tbl.atac.psp, color="purple")
   displayTrack(igv, track)

   track <- DataFrameAnnotationTrack("Mayo-CTL", tbl.atac.ctl, color="green")
   displayTrack(igv, track)


} # display.openChromatin
#---------------------------------------------------------------------------------------------------
ndufs2 <- function()
{
    tbl.gh <- retrieveEnhancersFromDatabase(ghdb, targetGene, tissues="all")
    tbl.gh$score <- tbl.gh$combinedscore
    track <- DataFrameQuantitativeTrack("GH", tbl.gh[, c("chrom", "start", "end", "score")],
                                    autoscale=TRUE, color="brown")
    displayTrack(igv, track)

    tbl.atac <- human.brain.mayo.oc()
    track <- DataFrameAnnotationTrack("MayoATAC", tbl.atac, color="darkgreen")
    displayTrack(igv, track)

    tbl.atac <- human.brain.consensus.boca.oc()
    track <- DataFrameAnnotationTrack("bocaATAC", tbl.atac, color="darkblue")
    displayTrack(igv, track)


    tbl.sub <- subset(tbl.fimo, p.value < 1e-4)
    tbl.sub$score <- -log10(tbl.sub$p.value)
    head(as.data.frame(sort(table(tbl.sub$tf), decreasing=TRUE)))
    track <- DataFrameQuantitativeTrack("fimo", tbl.sub[, c("chrom", "start", "end", "score")],
                                    autoscale=TRUE, color="blue")
    displayTrack(igv, track)

} # ndufs2
#----------------------------------------------------------------------------------------------------
logicFS.tracks <- function()
{
   rsids.top.10 <- c("rs2070901", "rs2070902", "rs4575098", "rs17356051", "rs145282062")

   rsids.top.100 <- c("rs79181095", "rs4348741", "rs4489574", "rs148185634",
                       "rs57213460", "rs2501873", "rs112300525",
                      "rs11579514", "rs6671288", "rs2070902", "rs3003596")

   tbl.track <- unique(subset(tbl.eqtl, rsid %in% rsids.top.100)[, c("chrom", "hg38", "hg38", "rsid")])
   dim(unique(tbl.track))
   colnames(tbl.track) <- c("chrom", "start", "end", "rsid")
   tbl.track$start <- tbl.track$start - 1
   track <- DataFrameAnnotationTrack("logicFS-100", tbl.track, color="red")
   displayTrack(igv, track)

} # logicFS.tracks
#----------------------------------------------------------------------------------------------------
ndufs2.fimo <- function()
{
  tbl.fimo <- get(load("tbl.fimo.NDUFS2.RData"))
  subset(tbl.fimo, start < 161185602 & end > 161185602)

} # ndufs2.fimo
#----------------------------------------------------------------------------------------------------
tms.and.trena <- function(mtx.rna)
{
   mtx.rna <- get.rna.matrix("old-rosmap")
   mtx.rna <- get.rna.matrix("max-rosmap")
   mtx.rna <- get.rna.matrix("sage-eqtl-rosmap")
   state$mtx.rna <- mtx.rna

   #mtx.rna <- get.rna.matrix("sage-eqtl-rosmap")
   fimo.file <- "tbl.fimo.1e-3.1Mspan.RData"
   tbl.fimo <- get(load(fimo.file))
   dim(tbl.fimo)
   state$tbl.fimo <- tbl.fimo
   #tbl.atac <-  human.brain.consensus.boca.oc()
   tbl.atac <- human.brain.mayo.oc()
   state$tbl.mayo.atac <- tbl.atac

   trenaProject <- TrenaProjectAD()
   tms <- TMS$new(trenaProject, targetGene, tbl.fimo, tbl.atac, quiet=FALSE)
   tms$scoreFimoTFBS()   # chip, conservation, genehancer, genic annotations, distance to tss
   tms$add.tf.mrna.correlations(mtx.rna, featureName="cor.all")

   tbl.tms <- unique(tms$getTfTable())
   dim(tbl.tms)
   checkTrue(nrow(tbl.tms)/nrow(tbl.fimo) > 0.99)
   save(tbl.tms, file="tbl.tms.ndufs2.mayo-atac.RData")

   fivenum(tbl.tms$cor.all)    # old rosm: -0.43 -0.17 -0.05  0.06  0.43
                               # max rosmap: -0.37 -0.13 -0.04  0.04  0.32
   dim(subset(tbl.tms, abs(cor.all) > 0.42))

   dim(subset(tbl.tms, fimo_pvalue <= 1e-6 & chip))
   dim(subset(tbl.tms, fimo_pvalue < 1e-5 & abs(cor.all) > 0.3 & chip & gh > 200))
   roi <- getGenomicRegion(igv)
   dim(subset(tbl.tms, start >= roi$start & end <= roi$end & abs(cor.all) > 0.2)) # & chip))

   fimo.threshold <- 1e-3
   tbl.tms.final <- subset(tbl.tms, fimo_pvalue < fimo.threshold & abs(cor.all) > 0.2 & gh > 0)
   dim(tbl.tms.final)
   tfs <- unique(tbl.tms.final$tf)
   printf("candidate regulators: %d", length(tfs))
   tbl.trena <- tms$build.trena.model(tfs, mtx.rna)
   tbl.trena$targetGene <- targetGene

   new.order <- order(abs(tbl.trena$spearmanCoeff), decreasing=TRUE)
   new.order <- order(abs(tbl.trena$rfScore), decreasing=TRUE)

   tbl.trena <- tbl.trena[new.order,]
   rownames(tbl.trena) <- NULL
   state$tbl.trena <- tbl.trena
   head(tbl.trena, n=20)

   tbl.pretty <- head(tbl.trena, n=20)[, -(c(8,9))]
   for(column.name in colnames(tbl.pretty)[-c(1,8)])
       tbl.pretty[[column.name]] <- round(tbl.pretty[[column.name]], digits=3)

   gr.pec.gh <- reduce(c(GRanges(state$tbl.gh), c(GRanges(state$tbl.pec.enhancers))))
   sum(width(gr.pec.gh))   # 77k
   tbl.pec.gh <- as.data.frame(gr.pec.gh)
   dim(tbl.pec.gh)   # 35 5
   tbl.ov <- as.data.frame(findOverlaps(GRanges(tbl.tms), gr.pec.gh))
   tbl.tms.gh.pec <- tbl.tms[tbl.ov$queryHits,]
   dim(tbl.tms)          #  180440
   dim(tbl.tms.gh.pec)   # 2101878
   viz <- FALSE

   tfbs.counts <- unlist(lapply(tbl.trena$gene, function(TF) nrow(subset(tbl.tms.gh.pec, tf==TF))))
   tbl.trena$tfbs <- tfbs.counts

   if(viz){
      tfs <- head(tbl.trena$gene, n=20)
      set.seed(31)
      fimo.threshold <- 1e-3
      for(TF in tfs){
         tbl.track <- subset(tbl.tms.gh.pec, tf==TF & fimo_pvalue < fimo.threshold)[, c("chrom", "start", "end")]
         track <- DataFrameAnnotationTrack(TF, tbl.track, color="random")
         displayTrack(igv, track)
         } # for tf
      } # if viz

} # trena
#----------------------------------------------------------------------------------------------------
run.many.trenas <- function()
{
   tbl.map <- etx$getIdMap()
   models <- list()


   key <- "old-rosmap"
   mtx.old <- get.rna.matrix(key)
   mtx.old.ptColnames <- mtx.old
   indices <- match(colnames(mtx.old), tbl.map$sample)
   patient.names <- tbl.map$patient[indices]
   colnames(mtx.old.ptColnames) <- patient.names

   models[[key]] <- tms.and.trena(mtx.old.ptColnames)

   key <- "sage-eqtl-rosmap"
   mtx.eqtl <- get.rna.matrix(key)
   # mtx.eqtl[which(is.na(mtx.eqtl))] <- 0

   models[[key]] <- tms.and.trena(mtx.eqtl)

  x <- etx$splitExpressionMatrixByMutationStatusAtRSID(mtx.old.ptColnames, tag.snp, study.name="rosmap")
  names(x)
  x$genotypes.rna
  key <- "old-rosmap-wtOnly"
  models[[key]] <- tms.and.trena(x$wt)

  key <- "old-rosmap-mutOnly"
  models[[key]] <- tms.and.trena(x$mut)

     # these old-rosmap-(mut|wt)Only have very different models
     #         gene   betaLasso    betaRidge spearmanCoeff pearsonCoeff   rfScore      xgboost class tfbs
     #  wt  2  HSF2  0.14907282   0.066276095   0.4096720  0.4333165    21.262856  0.101314518   tf    1
     # mut 19  HSF2  0.126890751  0.055976412   0.3481374  0.4078219     6.510664  0.037723865   tf    1
     # t.test(x$mut["NDUFS2",], x$wt["NDUFS2",])$p.value  # [1] 0.002199789
     # t.test(x$mut["HSF2",], x$wt["HSF2",])$p.value      # [1] 0.06444814

   save(models, file="ndufs2.rosmp.trena.models.RData")
   # next up: look at broken motifs for these TFs

} # run.many.trenas
#----------------------------------------------------------------------------------------------------
browse.breaks <- function()
{
   print(load("motif.breaks.2119.variantsno.pvals.Sun-Sep-26-10:10:16-2021.RData"))
   tbl.sub <- subset(tbl.breaks, SNP_id=="rs4575098")
   as.data.frame(sort(table(tbl.sub$geneSymbol), decreasing=TRUE))

} # browse.breaks
#----------------------------------------------------------------------------------------------------
# cory suggests: split the samples based on the variant and looked at differences in
# correlation with HSF2?
#  chr1  161185602 rs4575098 4.837459e-06 ENSG00000158864     NDUFS2   ampad-mayo    cer unknown
#  chr1  161185602 rs4575098 8.217797e-08 ENSG00000158864     NDUFS2 ampad-rosmap  dlpfc unknown

test_eQTL_split <- function()
{
   mtx <- get.rna.matrix("sage-eqtl-cer")
   #mtx <- get.rna.matrix("old-mayo-tcx")
   dim(mtx)
   cor(mtx["NDUFS2",], mtx["HSF2",], method="spearman", use="pairwise.complete") # 0.8777
   col.names <- colnames(mtx)
   col.names <- sub("^X", "", col.names)
   col.names <- sub("_TCX", "", col.names)
   colnames(mtx) <- col.names

   x <- etx$splitExpressionMatrixByMutationStatusAtRSID(mtx, rsid, study.name="mayo")
   x$genotypes.rna

   cor(x$mut["HSF2",], x$mut["NDUFS2",], method="spearman", use="pairwise.complete")
   cor(x$wt["HSF2",], x$wt["NDUFS2",], method="spearman", use="pairwise.complete")

   dir <- "~/github/TrenaProjectAD/prep/rna-seq-counts-from-synapse/eqtl"
   file.rosmap <- "mtx.rosmap.rnaseq-residual-eqtl-geneSymbols-patients-15582x632.RData"
   mtx <- get(load(file.path(dir, file.rosmap)))

   targetGene <- "NDUFS2"
   rsid <- "rs4575098"

   etx <- EndophenotypeExplorer$new(targetGene, "hg38")
   x <- etx$splitExpressionMatrixByMutationStatusAtRSID(mtx, rsid, study.name="rosmap")
   x$genotypes.rna

   dim(x$mut)
   dim(x$wt)
   dim(mtx)
   cor(x$mut["HSF2",], x$mut["NDUFS2",], method="spearman", use="pairwise.complete")
   cor(x$wt["HSF2",], x$wt["NDUFS2",], method="spearman", use="pairwise.complete")

} # trena_eQTL_split
#----------------------------------------------------------------------------------------------------
# intersect potent eQTLs with atac bearing motifs from high ranked tfs in ndufs2 trena model
eqtl.atac.tf.motif <- function()
{
   if(!exists("models"))
     load("ndufs2.rosmp.trena.models.RData")
   tbl.fimo <- get(load("tbl.fimo.NDUFS2.RData"))
   tbl.atac <- human.brain.mayo.oc()
   tbl.tms <- get(load("tbl.tms.ndufs2.mayo-atac.RData"))
   #tbl.atac <-  human.brain.consensus.boca.oc()

   roi <- getGenomicRegion(igv)

   tbl.atac.roi <- subset(tbl.atac, chrom==roi$chrom & start >= roi$start & end <= roi$end)
   track <- DataFrameAnnotationTrack("mayo atac", tbl.atac.roi, color="darkgreen")
   displayTrack(igv, track)

   tbl.trena <- models[["old-rosmap-wtOnly"]]
   top.tfs <- tbl.trena$gene[1:5]
   top.tfs   # "MEIS2" "HSF2"  "NFIA"  "TBR1"  "RFX5"

   tbl.tms.hsf2 <- subset(tbl.tms, tf=="HSF2")   # 29 18

   tbl.eqtl.ndufs2 <- subset(tbl.eqtl, genesymbol=="NDUFS2")[, c("chrom", "hg38", "hg38", "pvalue")]
   colnames(tbl.eqtl.ndufs2) <- c("chrom", "start", "end", "score")
   tbl.eqtl.ndufs2$start <-    tbl.eqtl.ndufs2$start - 1
   tbl.eqtl.ndufs2$score <- -log10(tbl.eqtl.ndufs2$score)
   gr.eqtl.ndufs2 <- GRanges(tbl.eqtl.ndufs2)

   gr.atac <- GRanges(tbl.atac.roi)
   tbl.ov <- as.data.frame(findOverlaps(gr.eqtl.ndufs2, gr.atac))
   dim(tbl.ov)
   tbl.atac.with.ndufs2.eqtls <- unique(tbl.atac.roi[tbl.ov$subjectHits,])

   track <- DataFrameAnnotationTrack("atac.ndufs.eqtls", tbl.atac.with.ndufs2.eqtls, col="red")
   displayTrack(igv, track)


   tbl.ov <-
       as.data.frame(findOverlaps(GRanges(tbl.tms.hsf2), GRanges(tbl.atac.with.ndufs2.eqtls)))

   tbl.ov <-
       as.data.frame(findOverlaps(GRanges(tbl.tms.hsf2), GRanges(tbl.atac.with.ndufs2.eqtls[8,])))

   tbl.tms.hsf2.atac.eqtl <- tbl.tms.hsf2[tbl.ov$queryHits,]
   tbl.track <- subset(tbl.tms.hsf2.atac.eqtl,  phast7 >= 0)
   track <- DataFrameAnnotationTrack("hsf2 tfbs", tbl.track, color="purple")
   displayTrack(igv, track)

   track <- DataFrameAnnotationTrack("hsf2 all", tbl.tms.hsf2, color="orange")
   displayTrack(igv, track)


} # eqtl.atac.tf.motif
#----------------------------------------------------------------------------------------------------
identify.tfbs.in.atac.regions.with.eQTLS <- function(tbl.tms, tbl.atac, tbl.eqtls, TF, target,
                                                     display=FALSE)
{
    tbl.tms.tf <- subset(tbl.tms, tf==TF)   # 29 18
    dim(tbl.tms.tf)

   tbl.eqtl.target <- subset(tbl.eqtl, genesymbol==target)[, c("chrom", "hg38", "hg38", "pvalue")]
   dim(tbl.eqtl.target)
   colnames(tbl.eqtl.target) <- c("chrom", "start", "end", "score")
   tbl.eqtl.target$start <-    tbl.eqtl.target$start - 1
   tbl.eqtl.target$score <- -log10(tbl.eqtl.target$score)
   gr.eqtl.target <- GRanges(tbl.eqtl.target)

   gr.atac <- GRanges(tbl.atac)
   tbl.ov <- as.data.frame(findOverlaps(gr.eqtl.target, gr.atac))
   dim(tbl.ov)
   tbl.atac.with.target.eqtls <- unique(tbl.atac[tbl.ov$subjectHits,])

   #if(display){
   #    track <- DataFrameAnnotationTrack(sprintf("atac.%s.eqtls", target),
   #                                      tbl.atac.with.target.eqtls, col="red")
   #    displayTrack(igv, track)
   #    }

   tbl.ov <-
       as.data.frame(findOverlaps(GRanges(tbl.tms.tf), GRanges(tbl.atac.with.target.eqtls)))

   tbl.tms.tfbs.in.atac.containing.eqtl <- tbl.tms.tf[tbl.ov$queryHits,]

   if(display){
      tbl.track <- subset(tbl.tms.tfbs.in.atac.containing.eqtl)
      trackName <- sprintf("%s tfbs", TF)
      track <- DataFrameAnnotationTrack(trackName, tbl.track, color="purple", trackHeight=24)
      displayTrack(igv, track)
      }

    tbl.tms.tfbs.in.atac.containing.eqtl

} # identify.tfbs.in.atac.regions.with.eQTLS
#----------------------------------------------------------------------------------------------------
test_identify.tfbs.in.atac.regions.with.eQTLS <- function()
{
   message(sprintf("--- test_identify.tfbs.in.atac.regions.with.eQTLS"))

   if(!exists("models"))
     load("ndufs2.rosmp.trena.models.RData")

   tbl.fimo <- get(load("tbl.fimo.NDUFS2.RData"))
   tbl.atac <- human.brain.mayo.oc()
   tbl.tms <- get(load("tbl.tms.ndufs2.mayo-atac.RData"))
   targetGene <- "NDUFS2"
   tbl.eqtl <- get(load("tbl.eqtl+-600kb-rs4575098.RData"))

   tbl.tms.tfbs.hsf2 <- identify.tfbs.in.atac.regions.with.eQTLS(tbl.tms, tbl.atac, tbl.eqtls,
                                                            "HSF2", "NDUFS2")
   checkEquals(nrow(tbl.tms.tfbs.hsf2), 14)

   tbl.tms.tfbs.meis2 <- identify.tfbs.in.atac.regions.with.eQTLS(tbl.tms, tbl.atac, tbl.eqtls,
                                                                  "MEIS2", "NDUFS2")
   checkEquals(nrow(tbl.tms.tfbs.meis2), 8)


   tbl.tms.tfbs.nfix <- identify.tfbs.in.atac.regions.with.eQTLS(tbl.tms, tbl.atac, tbl.eqtls,
                                                                  "NFIX", "NDUFS2", display=TRUE)
   checkEquals(nrow(tbl.tms.tfbs.nfix), 8)

   tfs <- c("MEIS2", "HSF2", "NFIA", "TBR1", "RFX5", "NFIX")
   for(tf in tfs)
      tbl.tms.tfbs <- identify.tfbs.in.atac.regions.with.eQTLS(tbl.tms, tbl.atac, tbl.eqtls,
                                                                  tf, "NDUFS2", display=TRUE)
  tfs.all <- unique(mcols(query(MotifDb, "sapiens", c("jaspar2018")))$geneSymbol)


} # test_identify.tfbs.in.atac.regions.with.eQTLS
#----------------------------------------------------------------------------------------------------
display.all.precalucated.motifs.in.region <- function()
{
   fimo.file <- "tbl.fimo.NDUFS2.RData"
   tbl.fimo <- get(load(fimo.file))
   dim(tbl.fimo)
   state$tbl.fimo <- tbl.fimo

   roi <- getGenomicRegion(igv)
   tbl.fimo.roi <- subset(tbl.fimo, start >= roi$start & end <= roi$end)
   dim(tbl.fimo.roi)

} # display.all.precalucated.motifs.in.region
#----------------------------------------------------------------------------------------------------
query.all.motifs.in.region <- function()
{
  source("~/github/fimoService/batchMode/fimoBatchTools.R")
  meme.file <- "jaspar2018.meme"
  motifs <- query(MotifDb, c("sapiens"), c("jaspar2018", "HOCOMOCOv11-core-A"))

  length(motifs)
  export(motifs, con=meme.file, format="meme")
  roi <- getGenomicRegion(igv)
  tbl.roi <- with(roi, data.frame(chrom=chrom, start=start, end=end, stringsAsFactors=FALSE))
  printf("roi width: %d", 1 + roi$end - roi$start)
  tbl.fimo.roi <- fimoBatch(tbl.roi, matchThreshold=1e-3, genomeName="hg38", pwmFile=meme.file)
  dim(tbl.fimo.roi)
  state$tbl.fimo.roi <- tbl.fimo.roi
  roi.small <- getGenomicRegion(igv)
  tfs <- sort(unique(subset(tbl.fimo.roi, start >= roi.small$start & end <= roi.small$end)$tf))
  tbl.fimo.small <- subset(tbl.fimo.roi, start >= roi.small$start & end <= roi.small$end)
  length(tfs)
  set.seed(117)
  tfs.also.trena <- intersect(tfs, state$tbl.trena$gene)
  length(tfs.also.trena)
  for(TF in tfs.also.trena){
      tbl.track <- subset(tbl.fimo.small, tf==TF)[, c("chrom", "start", "end", "p.value")]
      tbl.track$p.value <- -log10(tbl.track$p.value)
      track <- DataFrameQuantitativeTrack(TF, tbl.track, autoscale=FALSE, min=0, max=10, color="random")
      displayTrack(igv, track)
      }

  save(tbl.fimo, file="tbl.fimo.1e-3.jaspar2018.RData")
  threshold <- 1e-3
  tbl.fimo.sub <- subset(tbl.fimo, p.value<=threshold)
  tbl.freq <- as.data.frame(sort(table(tbl.fimo.sub$tf), decreasing=TRUE), drop=FALSE)
  tbl.freq
  colnames(tbl.freq) <- c("tf", "count")
  top.tfs <- head(models[["old-rosmap"]]$gene, n=10)
  dim(tbl.freq)
  subset(tbl.freq, tf %in% top.tfs)

} # query.all.motifs.in.region
#----------------------------------------------------------------------------------------------------
display.nott.interactomes <- function()
{
  data.dir <- "~/github/TrenaProjectAD/explore/enhancers/nott"
  file <- "interactomes-oligo-neuronal-microglial.RData"
  full.path <- file.path(data.dir, file)
  stopifnot(file.exists(full.path))
  print(load(full.path))
     # "tbl.microglialInteractome" "tbl.neuronalInteractome"   "tbl.oligoInteractome"

  roi <- getGenomicRegion(igv)
  interactomes <- list(mg=tbl.microglialInteractome,
                       n=tbl.neuronalInteractome,
                       o=tbl.oligoInteractome)
  for(name in names(interactomes)){
      tbl.data <- interactomes[[name]]
      title <- name
      tbl.i <- subset(tbl.data,
                     chr1==roi$chrom & start1 >= roi$start & end1 <= roi$end &
                     chr2==roi$chrom & start2 >= roi$start & end2 <= roi$end)
      dim(tbl.i)
      track <- BedpeInteractionsTrack(name, tbl.i, trackHeight=200)
      displayTrack(igv, track)
      tbl.p <- unique(tbl.i[, c("chr2", "start2", "end2")])
      title <- sprintf("%s.promoter", name)
      track <- DataFrameAnnotationTrack(title, tbl.p, color="darkGreen")
      displayTrack(igv, track)
      }

} # display.nott.interactomes
#----------------------------------------------------------------------------------------------------
display.pec.enhancers <- function()
{

  f <- "~/github/TrenaProjectAD/explore/psychEncode/psychEncode.enhancers.all.RData"
  tbl.pec.enhancers <- get(load(f))
  roi <- getGenomicRegion(igv)
  tbl.pe.sub <- subset(tbl.pec.enhancers, chrom==roi$chrom & start > roi$start & end < roi$end)
  state$tbl.pec.enhancers <- tbl.pe.sub
  track <- DataFrameAnnotationTrack("pec enhancers", tbl.pe.sub, color="blue")
  displayTrack(igv, track)

  tbl.e <- unique(tbl.i.mg[, c("chr1", "start1", "end1")])
  track <- DataFrameAnnotationTrack("mg enhancer", tbl.e, color="red")
  displayTrack(igv, track)

} # display.microglial.enhancers
#----------------------------------------------------------------------------------------------------
# mayo atac seq, 22 samples, foound here:
# ~/github/TrenaProjectAD/explore/mayo-epigenetics/atac/readAllPeaks.R
# diffbind-earlyExploration.R
tbl.adoc <- get(load("~/github/TrenaProjectAD/explore/mayo-epigenetics/atac/tbl.summedScoresATAC-seq-AD.RData"))
tbl.pspoc <- get(load("~/github/TrenaProjectAD/explore/mayo-epigenetics/atac/tbl.summedScoresATAC-seq-PSP.RData"))
tbl.ctloc <- get(load("~/github/TrenaProjectAD/explore/mayo-epigenetics/atac/tbl.summedScoresATAC-seq-CTL.RData"))

tbls <- list(tbl.adoc, tbl.pspoc, tbl.ctloc)
names(tbls) <- c("AD", "PSP", "CTL")

roi <- getGenomicRegion(igv)

for(name in names(tbls)){
    tbl.sub <- subset(tbls[[name]], chrom==roi$chrom & start >= roi$start & end <= roi$end)
    track <- DataFrameQuantitativeTrack(name,
                                        tbl.sub[, c("chrom", "start", "end", "score")],
                                        autoscale=FALSE, color="random", min=0, max=500)
    #browser()
    displayTrack(igv, track)
    }


