#------------------------------------------------------------------------------------------------------------------------
library(EndophenotypeExplorer)
library(TrenaProjectAD)
library(RPostgreSQL)
library(plyr)
library(ghdb)
source("~/github/endophenotypeExplorer/R/getExpressionMatrices.R")
source("~/github/TrenaMultiScore/tools/runner/v2/tmsCore.R")
library(RUnit)
#------------------------------------------------------------------------------------------------------------------------
targetGene <- "NDUFS2"
if(!exists("etx")){
   etx <- EndophenotypeExplorer$new(targetGene, "hg38", initialize.snpLocs=TRUE)
   ghdb <- GeneHancerDB()
   }

tag.snp <- "rs4575098"
ld.snp  <- "rs11585858"

if(!exists("igv")){
   igv <- start.igv(targetGene, "hg38")
   zoomOut(igv)
   zoomOut(igv)
   }
#------------------------------------------------------------------------------------------------------------------------
display.snps <- function()
{
      #------------------------------------------------------------
      # first: a track with rsids, visible when clicked upon
      #------------------------------------------------------------

   tbl.hap <- read.table("../haploreg-rs4575098-0.2.tsv", sep="\t", as.is=TRUE, header=TRUE, nrow=-1)
   tbl.track <- tbl.hap[, c("chrom", "hg38", "hg38", "rsid")]
   colnames(tbl.track)[c(2,3)] <- c("start", "end")
   tbl.track$start <- tbl.track$start - 1
   track <- DataFrameAnnotationTrack("hapName", tbl.track, color="brown")
   displayTrack(igv, track)
   shoulder <- 5000
   roi <- sprintf("%s:%d-%d", tbl.track$chrom[1], min(tbl.track$start)-shoulder, max(tbl.track$end)+shoulder)
   showGenomicRegion(igv, roi)

      #------------------------------------------------------------
      # now a track with rSquared magnitudes
      #------------------------------------------------------------

   tbl.track <- tbl.hap[, c("chrom", "hg38", "hg38", "rSquared")]
   colnames(tbl.track)[c(2,3)] <- c("start", "end")
   tbl.track$start <- tbl.track$start - 1
   track <- DataFrameQuantitativeTrack("hapScore", tbl.track, color="darkRed", autoscale=FALSE, min=0, max=1)
   displayTrack(igv, track)

} # display.snps
#------------------------------------------------------------------------------------------------------------------------
get.eqtls <- function()
{
    db <- dbConnect(PostgreSQL(), user="trena", password="trena", dbname="genereg2021", host="khaleesi")
    dbGetQuery(db, "select * from eqtls limit 3")
    getGenomicRegion(igv)

    roi <- getGenomicRegion(igv)
    printf("%d", with(roi, end-start))
    query <- with(roi, sprintf("select * from eqtls where chrom='%s' and hg38 > %d and hg38 < %d",
                               chrom, start, end))
    tbl.eqtl <- dbGetQuery(db, query)
    dim(tbl.eqtl)
    tbl.eqtl <- subset(tbl.eqtl, study=="ampad-rosmap")
    tbl.eqtl <- subset(tbl.eqtl, pvalue < 0.05)
    tbl.freq <- as.data.frame(sort(table(tbl.eqtl$genesymbol), decreasing=TRUE))
    all.genes <- as.character(tbl.freq$Var1)
    for(gene in "NDUFS2"){
       tbl.sub <- subset(tbl.eqtl, genesymbol==gene)
       tbl.track <- tbl.sub[, c("chrom", "hg38", "hg38", "pvalue")]
       colnames(tbl.track) <- c("chrom", "start", "end", "score")
       tbl.track$start <- tbl.track$start - 1
       tbl.track$score <- -log10(tbl.track$score)
       track <- DataFrameQuantitativeTrack(sprintf("%s-eQTL", gene),
                                           tbl.track, autoscale=FALSE,
                                           min=0, max=10, color="random")
       displayTrack(igv, track)
       }

   #    dim(tbl.eqtl)
   #    head(tbl.eqtl)
   #    tbl.eqtl.hg38 <- unique(tbl.eqtl[, c("chrom", "hg38", "rsid", "genesymbol", "pvalue", "study", "tissue")])

    # target.genes <- sort(unique(tbl.eqtl.hg38$genesymbol))
    # deleters <- which(nchar(target.genes) == 0)
    # deleters
    # if(length(deleters) > 0){
    #     printf("deleting %d empty string gene names", length(deleters))
    #     target.genes <- target.genes[-deleters]
    #     }
    # length(target.genes)

   #    pval.max <- 1e-5
   #
   #    for(gene in targetGene){
   #       tbl.track <- subset(tbl.eqtl, genesymbol==gene & pvalue < pval.max)
   #       if(nrow(tbl.track) == 0) next;
   #       tbl.track$end <- tbl.track$hg38
   #       tbl.track$start <- tbl.track$hg38 -1
   #       tbl.track <- tbl.track[, c("chrom", "start", "end", "rsid")]
   #       track <- DataFrameAnnotationTrack(gene, tbl.track, color="random", trackHeight=25)
   #       displayTrack(igv, track)
   #       } # for target.gene
   #
   #    genes.assoc <- c("PPOX","APOA2","ADAMTS4","B4GALT3","USP21","NDUFS2","TOMM40L")
   #    intersect(target.genes, genes.assoc)
   #
   #    tbl.eqtls.oi <- subset(tbl.eqtl.hg38, pvalue <= 1e-6)  # NDUFS2 PCP4L1 TSTD1
   #    new.order <- order(tbl.eqtls.oi$rsid)
   #    tbl.eqtls.oi <- tbl.eqtls.oi[new.order,]
   #    dim(tbl.eqtls.oi)  # 32 7
   #    save(tbl.eqtls.oi, file="AD-eqtls-near-rs4575098-threshold-1e-6.RData")
   #
   #    tbl.eqtls.oi <- subset(tbl.eqtl.hg38, pvalue <= 1e-5)  # F11R NDUFS2 PCP4L1 PPOX  TSTD1
   #    new.order <- order(tbl.eqtls.oi$rsid)
   #    tbl.eqtls.oi <- tbl.eqtls.oi[new.order,]
   #    dim(tbl.eqtls.oi)  # 90 7
   #    save(tbl.eqtls.oi, file="AD-eqtls-near-rs4575098-threshold-1e-5.RData")
   #
   #    tbl.eqtls.oi <- subset(tbl.eqtl.hg38, pvalue <= 1e-4)  # F11R NDUFS2 PCP4L1 PPOX  TSTD1
   #    new.order <- order(tbl.eqtls.oi$rsid)
   #    tbl.eqtls.oi <- tbl.eqtls.oi[new.order,]
   #    dim(tbl.eqtls.oi)  # 117 7
   #    save(tbl.eqtls.oi, file="AD-eqtls-near-rs4575098-threshold-1e-4.RData")

      # next up on khaleesi: build tissue/study specific trena models for each set of 3 then 5 genes
      # break motifs for all rsids

} # get.eqtls
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
    track <- DataFrameAnnotationTrack("bocaATAC", tbl.atac, color="darkgreen")
    displayTrack(igv, track)

    tbl.fimo <- get(load("tbl.fimo.NDUFS2.RData"))
    dim(tbl.fimo)  # 29115

    tbl.sub <- subset(tbl.fimo, p.value < 1e-4)
    tbl.sub$score <- -log10(tbl.sub$p.value)
    head(as.data.frame(sort(table(tbl.sub$tf), decreasing=TRUE)))
    track <- DataFrameQuantitativeTrack("fimo", tbl.sub[, c("chrom", "start", "end", "score")],
                                    autoscale=TRUE, color="blue")
    displayTrack(igv, track)

} # ndufs2
#----------------------------------------------------------------------------------------------------
ndufs2.fimo <- function()
{
  tbl.fimo <- get(load("tbl.fimo.NDUFS2.RData"))
  subset(tbl.fimo, start < 161185602 & end > 161185602)

} # ndufs2.fimo
#----------------------------------------------------------------------------------------------------
tms.and.trena <- function(mtx.rna)
{
   #mtx.rna <- get.rna.matrix("sage-eqtl-rosmap")
   tbl.fimo <- get(load("tbl.fimo.NDUFS2.RData"))
   tbl.atac <-  human.brain.consensus.boca.oc()
   #tbl.atac <- human.brain.mayo.oc()

   trenaProject <- TrenaProjectAD()
   tms <- TMS$new(trenaProject, targetGene, tbl.fimo, tbl.atac, quiet=FALSE)
   tms$scoreFimoTFBS()   # chip, conservation, genehancer, genic annotations, distance to tss
   tms$add.tf.mrna.correlations(mtx.rna, featureName="cor.all")

   tbl.tms <- tms$getTfTable()
   checkEquals(nrow(tbl.tms), nrow(tbl.fimo))

   subset(tbl.tms, abs(cor.all) > 0.42)

   dim(subset(tbl.tms, fimo_pvalue <= 1e-6 & chip))
   dim(subset(tbl.tms, fimo_pvalue < 1e-5 & abs(cor.all) > 0.4 & chip & gh > 200))

   tbl.tms.final <- subset(tbl.tms, fimo_pvalue < 1e-4 & abs(cor.all) > 0.1 & gh > 1)
   tfs <- unique(tbl.tms.final$tf)
   printf("candidate regulators: %d", length(tfs))
   tbl.trena <- tms$build.trena.model(tfs, mtx.rna)
   new.order <- order(abs(tbl.trena$spearmanCoeff), decreasing=TRUE)
   tfbs.counts <- unlist(lapply(tbl.trena$gene, function(TF) nrow(subset(tbl.tms.final, tf==TF))))
   tbl.trena$tfbs <- tfbs.counts
   tbl.trena$targetGene <- targetGene
   tbl.trena <- tbl.trena[new.order,]
   tbl.trena

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

   # next up: look at broken motifs for these TFs

} # run.many.trenas
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


} # trena
#----------------------------------------------------------------------------------------------------
