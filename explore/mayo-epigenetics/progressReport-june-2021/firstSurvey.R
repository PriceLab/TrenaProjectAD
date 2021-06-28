library(RPostgreSQL)
if(!exists("eqtl.db"))
    eqtl.db <- dbConnect(PostgreSQL(), user="trena", password="trena", dbname="genereg2021", host="khaleesi")

tbl.38 <- get(load("~/github/TrenaProjectAD/inst/extdata/gwasLoci/tbl.posthuma-38-loci-curated.RData"))
tbl.assoc <- get(load("~/github/TrenaProjectAD/inst/extdata/gwasLoci/tbl.posthuma-38-geneAssociations-curated-3828x12.RData"))
tbl.atac <- get(load("~/github/TrenaProjectAD/explore/mayo-epigenetics/atac/dbaConsensusRegionsScored.74273x30.RData"))

    #------------------------------------------------------------
    # fatten up the atac-peaks, realistically I hope
    #------------------------------------------------------------

gr.atac <- GRanges(seqnames=tbl.atac$chrom, IRanges(start=tbl.atac$start, end=tbl.atac$end))
atac.widths <- width(gr.atac)
head(atac.widths)
atac.padding <- as.integer(round(atac.widths/10))
start(gr.atac) <- (start(gr.atac) - atac.padding)
end(gr.atac) <- (end(gr.atac) + atac.padding)
atac.widths <- width(gr.atac)
head(atac.widths)



library(ghdb)
ghdb <- GeneHancerDB()
                                        #targetGene <- "BIN1"
stopifnot(exists("targetGene"))

if(!exists("igv"))
   igv <- start.igv(targetGene, "hg38")

tbl.eqtls.all <- get(load("tbl.eqtls.all38genes.RData"))
#----------------------------------------------------------------------------------------------------
get.all.eqtls <- function()
{
   selector <- sprintf("('%s')", paste(tbl.38$gene, collapse="','"))
   query <- sprintf("select * from eqtls where genesymbol in %s", selector)
   tbl.eqtls <- dbGetQuery(eqtl.db, query)
   tbl.eqtls$score <- -log10(tbl.eqtls$pvalue)
   as.data.frame(sort(table(tbl.eqtls$genesymbol)))
   save(tbl.eqtls, file="tbl.eqtls.all38genes.RData")

} # get.all.eqtls
#----------------------------------------------------------------------------------------------------
get.eqtls <- function(gene, pval.threshold)
{
   subset(tbl.eqtls.all, genesymbol==gene & pvalue <= pval.threshold)

} # get.eqtls
#----------------------------------------------------------------------------------------------------
get.eqtls.db <- function(gene, pval.threshold)
{
  query <- sprintf("select * from eqtls where genesymbol='%s' and pvalue <= %f", gene, pval.threshold)
  tbl <- dbGetQuery(eqtl.db, query)
  invisible(tbl)

} # get.eqtls.db
#----------------------------------------------------------------------------------------------------
find.variants.in.atac <- function(targetGene, viz=FALSE)
{
  if(viz) showGenomicRegion(igv, targetGene)

  tbl.gh <- retrieveEnhancersFromDatabase(ghdb, targetGene, tissues="all")

  if(viz){
      tbl.track <- tbl.gh[, c("chrom", "start", "end", "combinedscore")]
      colnames(tbl.track)[4] <- "score"
      tbl.track$score <- asinh(tbl.track$score)
      track <- DataFrameQuantitativeTrack("gh", tbl.track, autoscale=TRUE, color="red")
      displayTrack(igv, track)
      }

  shoulder <- 10000
  roi <- sprintf("%s:%d-%d",
                 tbl.gh$chrom[1],
                 min(tbl.gh$start) - shoulder,
                 max(tbl.gh$end) + shoulder)

  if(viz) showGenomicRegion(igv, roi)

  tbl.leadVariant <- subset(tbl.38, geneSymbol==targetGene)

  if(viz){
      track <- with(tbl.leadVariant, DataFrameAnnotationTrack(leadVariant,
                                                              data.frame(chrom=chrom, start=hg38-1, end=hg38,
                                                                         name=leadVariant,
                                                                         stringsAsFactors=FALSE),
                                                              color="black",
                                                              trackHeight=25))
      displayTrack(igv, track)
      }

  #locusNumber <- grep(targetGene, tbl.38$geneSymbol)
  tbl.assoc.sub <- subset(tbl.assoc, gene==targetGene)[, c("chrom", "hg38", "hg38", "rsid")]
  if(nrow(tbl.assoc.sub) == 0) return()
  colnames(tbl.assoc.sub) <- c("chrom", "start", "end", "rsid")
  deleters <- which(is.na(tbl.assoc.sub$start))
  if(length(deleters) > 0){
      printf("deleting %d/%d bad rows in tbl.assoc for locusNumber %d, geneSymbol %s",
             length(deleters), nrow(tbl.assoc.sub), locusNumber, targetGene)
      tbl.assoc.sub <- tbl.assoc.sub[-(deleters),]
      }
  tbl.assoc.sub$start <- tbl.assoc.sub$start - 1;


  if(viz){
      track <- DataFrameAnnotationTrack(sprintf("%s.assoc", targetGene),
                                        tbl.assoc.sub, color="black", trackHeight=25)

      displayTrack(igv, track)
      }

  if(targetGene == "HLA-DRB1") browser()

  gr.atac <- GRanges(tbl.atac)
  gr.assoc <- GRanges(tbl.assoc.sub)

  roi <- sprintf("chr%s:%d-%d", tbl.gh$chrom[1], min(tbl.gh$start), max(tbl.gh$end))
  gr.region <- GRanges(roi)
  tbl.ov <- as.data.frame(findOverlaps(gr.atac, gr.assoc))
  printf("%s associated variants in atac: %d/%d", targetGene, nrow(tbl.ov), nrow(tbl.assoc.sub))

  if(viz){
      track <- DataFrameAnnotationTrack("atac", tbl.atac[tbl.ov$queryHits,], color="brown", trackHeight=25)
      displayTrack(igv, track)
      }

} # find.variants.in.atac
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
runTrena <- function(targetGene, tbl.atac.gene, fimoThreshold)
{
  source("~/github/fimoService/batchMode/fimoBatchTools.R")
  mtx <- get(load("~/github/TrenaProjectAD/inst/extdata/expression/mayo.tcx.16969x262.covariateCorrection.log+scale.RData"))
  meme.file <- "jaspar2018.meme"
  motifs <- query(MotifDb, c("sapiens", "jaspar2018"))
  length(motifs)
  export(motifs, con=meme.file, format="meme")

  tbl.fimo <- fimoBatch(tbl.atac.gene[, c("chrom", "start", "end")],
                        matchThreshold=fimoThreshold, genomeName="hg38", pwmFile=meme.file)
  dim(tbl.fimo)

  tf.candidates <- intersect(tbl.fimo$tf, rownames(mtx))
  length(tf.candidates)

  library(trena)
  solver <- EnsembleSolver(mtx,
                           targetGene=targetGene,
                           candidateRegulators=tf.candidates,
                           solverNames=c("lasso", "Ridge", "Spearman", "Pearson", "RandomForest", "xgboost"))
  tbl.out <- run(solver)

  return(list(trena=tbl.out, fimo=tbl.fimo))

} # runTrena
#----------------------------------------------------------------------------------------------------
gene.demo <- function(targetGene)
{
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

  browser()
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
  if(nrow(tbl.eqtls) > 0){
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
     }

  tbl.atacPlusGH <- reduce(GRanges(tbl.atac.sub), GRanges(tbl.gh), drop.empty.ranges=TRUE, min.gapwidth=10)

  tbl.roi <- data.frame(chrom="chr2", start=127015068, end=127245690, stringsAsFactors=FALSE)

  x <- runTrena(targetGene, tbl.roi, 1e-3)
  tbl.trena <- x$tren
  #x <- runTrena(targetGene, tbl.atac.sub, 1e-2)
  #x <- runTrena(targetGene, tbl.atac.sub, 1e-2)
  tbl.trena <- x$trena
  tbl.fimo <- x$fimo
  tbl.trena <- tbl.trena[order(abs(tbl.trena$spearmanCoeff), decreasing=TRUE),]

  tbl.fimo.sub <- subset(tbl.fimo, tf %in% tbl.trena$gene[1:10])
  dim(tbl.fimo.sub)

  tbl.tfs <- as.data.frame(sort(table(tbl.fimo.sub$tf), decreasing=TRUE))
  for(TF in head(tbl.trena$gene, n=10)){
      tbl.fimo.sub <- subset(tbl.fimo, tf==TF)[, c("chrom", "start", "end")]
      track <- DataFrameAnnotationTrack(TF, tbl.fimo.sub, color="random", trackHeight=25)
      displayTrack(igv, track)
      } # for

  list(igv=igv, trena=tbl.trena, fimo=tbl.fimo)

} # gene.demo
#----------------------------------------------------------------------------------------------------
for(gene in tbl.38$geneSymbol){
    find.variants.in.atac(gene)
    }

#----------------------------------------------------------------------------------------------------
haploreg.bin1.rsids <- function()
{
   targetGene <- "BIN1"
   lead <- subset(tbl.38, geneSymbol==targetGene)$leadVariant
   assoc <- subset(tbl.assoc, gene==targetGene)$rsid
   tbl.eqtls <- get.eqtls(targetGene, 0.01)
   eqtls <- subset(tbl.eqtls, tissue=="tcx")$rsid
   length(lead)    # 1
   length(assoc)   # 30
   length(eqtls)   # 85

   rsids <- sort(unique(c(lead, assoc, eqtls)))
   length(rsids)   # 115
   write(rsids, sep="\n", file="~/stage/bin1.rsids")

   tbl.hr <- read.table("bin1-haploreg-cooked-01.txt", sep="\t", as.is=TRUE, header=TRUE, nrow=-1, fill=TRUE)
   tbl.hr$chrom <- paste0("chr", tbl.hr$chrom)
   save(tbl.hr, file="tbl.bin1.haplotype.snps.RData")
   dim(tbl.hr)
   track <- DataFrameAnnotationTrack("hr", tbl.hr[, c("chrom", "hg38", "hg38", "rsid")],
                                     color="blue", trackHeight=25)
   displayTrack(igv, track)


   gr.hr <- GRanges(seqnames=tbl.hr$chrom, ranges=IRanges(start=tbl.hr$hg38, end=tbl.hr$hg38))
   length(gr.hr)
   tbl.ov <- as.data.frame(findOverlaps(gr.hr, gr.atac))
   dim(tbl.ov)
   tbl.ov
   tbl.hr.ov <- tbl.hr[tbl.ov$queryHits,]
   dim(tbl.hr.ov)

   igv <- start.igv("BIN1", "hg38")
   tbl.gh <- retrieveEnhancersFromDatabase(ghdb, targetGene, tissues="all")
   showGenomicRegion(igv, sprintf("%s:%d-%d", tbl.gh$chrom[1], min(tbl.gh$start)-shoulder, max(tbl.gh$end)+shoulder))

   roi <- getGenomicRegion(igv)
   tbl.atac.sub <- subset(tbl.atac, chrom==roi$chrom & start >= roi$start & end <= roi$end)
   dim(tbl.atac.sub)  # 23 31
   track <- DataFrameAnnotationTrack("ATAC", tbl.atac.sub[, c("chrom", "start", "end")], color="red", trackHeight=25)
   displayTrack(igv, track)

   track <- DataFrameAnnotationTrack("hr-atac", tbl.hr.ov[, c("chrom", "hg38", "hg38", "rsid")],
                                     color="black", trackHeight=25)
   displayTrack(igv, track)
   tbl.fimo <- calcs[["BIN1"]]$fimo
   tbl.trena <- calcs[["BIN1"]]$trena

   tfs <- tbl.trena$gene[1:10]
   for(TF in head(tbl.trena$gene, n=10)){
      tbl.fimo.sub <- subset(tbl.fimo, tf==TF)[, c("chrom", "start", "end")]
      track <- DataFrameAnnotationTrack(TF, tbl.fimo.sub, color="random", trackHeight=25)
      displayTrack(igv, track)
      } # for


} # haploreg.bin1.rsids
#----------------------------------------------------------------------------------------------------


# AGRN variants in atac: 2/68
# deleting 1/22 bad rows in tbl.assoc for locusNumber 22, geneSymbol CR1
# CR1 variants in atac: 0/21
# FHL2 variants in atac: 7/97
# BIN1 variants in atac: 0/30
# INPP5D variants in atac: 1/28
# CLNK variants in atac: 1/150
# HAVCR2 variants in atac: 0/58
# TNIP1 variants in atac: 1/39
# TREM2 variants in atac: 0/176
# CD2AP variants in atac: 1/93
# EPHA1 variants in atac: 1/88
# TMEM106B variants in atac: 2/80
# PILRA variants in atac: 4/36
# CLU variants in atac: 3/14
# deleting 2/201 bad rows in tbl.assoc for locusNumber 22, geneSymbol SHARPIN
# SHARPIN variants in atac: 12/199
# CCDC6 variants in atac: 1/33
# ECHDC3 variants in atac: 2/147
# SORL1 variants in atac: 6/79
# MS4A6A variants in atac: 0/34
# MADD variants in atac: 2/38
# PICALM variants in atac: 0/31
# SLC24A4 variants in atac: 0/47
# FERMT2 variants in atac: 0/60
# deleting 2/152 bad rows in tbl.assoc for locusNumber 22, geneSymbol APH1B
# APH1B variants in atac: 3/150
# ADAM10 variants in atac: 1/73
# RNF43 variants in atac: 1/33
# ABI3 variants in atac: 2/105
# deleting 2/335 bad rows in tbl.assoc for locusNumber 22, geneSymbol ACE
# ACE variants in atac: 5/333
# GRN variants in atac: 54/838
# CHRNE variants in atac: 2/199
# ABCA7 variants in atac: 6/158
# CD33 variants in atac: 0/13
# deleting 1/75 bad rows in tbl.assoc for locusNumber 22, geneSymbol KIR3DL2
# KIR3DL2 variants in atac: 1/74
# NTN5 variants in atac: 5/79
# CASS4 variants in atac: 0/17
# APP variants in atac: 1/102


# zoomOut(igv)
# zoomOut(igv)
# zoomOut(igv)

# tbl.atac <- get(load("~/github/TrenaProjectAD/explore/mayo-epigenetics/atac/dbaConsensusRegionsScored.74273x30.RData"))
# roi <- getGenomicRegion(igv)
# tbl.atac.sub <- subset(tbl.atac, chrom==roi$chrom & start >= roi$start & end <= roi$end)
# dim(tbl.atac.sub)
#
# for(col in c("AD", "PSP", "CTL")){
#    tbl.track <- tbl.atac.sub[, c("chrom", "start", "end", col)]
#    fivenum(tbl.track$AD)
#    track <- DataFrameQuantitativeTrack(col, tbl.track, color="random", autoscale=FALSE, min=0, max=0.025)
#    displayTrack(igv, track)
#    }
#
#
#
# source("~/github/fimoService/batchMode/fimoBatchTools.R")
# meme.file <- "jaspar2018-hocomocoCore.meme"
# motifs <- query(MotifDb, c("sapiens"), c("jaspar2018", "HOCOMOCOv11-core-A"))
# length(motifs) # 718
#
# export(motifs, con=meme.file, format="meme")
#
# tbl.fimo <- fimoBatch(tbl.atac.sub[, c("chrom", "start", "end")],
#                       matchThreshold=1e-3,
#                       genomeName="hg38",
#                       pwmFile=meme.file)
# dim(tbl.fimo)
# head(tbl.fimo)
# tail(as.data.frame(sort(table(tbl.fimo$tf))))
#
# library(trena)
# candidate.tfs <- unique(tbl.fimo$tf)
# length(candidate.tfs)   # 535
# mtx <- get(load("~/github/TrenaProjectAD/inst/extdata/expression/mayo.tcx.16969x262.covariateCorrection.log+scale.RData"))
# dim(mtx)   # 16969 262
# length(setdiff(c(candidate.tfs, targetGene), rownames(mtx))) # 17/80
# candidate.tfs <- intersect(candidate.tfs, rownames(mtx))
# length(candidate.tfs)  # 370
#
# solver <- EnsembleSolver(mtx,
#                          targetGene=targetGene,
#                          candidateRegulators=candidate.tfs,
#                          solverNames=c("lasso", "Ridge", "Spearman", "Pearson", "RandomForest", "xgboost"),
#                          geneCutoff=1.0)
# tbl.out <- run(solver)
# tbl.out <- tbl.out[order(abs(tbl.out$spearmanCoeff), decreasing=TRUE),]
# dim(tbl.out)
# head(tbl.out, n=10)
#
# unlist(lapply(candidate.tfs, function(tf) length(grep(tf, tbl.fimo$tf))))
# dim(tbl.out)
#
# tbl.out$bs.3 <- unlist(lapply(tbl.out$gene,
#                               function(TF)
#                                   nrow(subset(tbl.fimo, tf==TF & p.value <= 1e-3))))
# tbl.out$bs.4 <- unlist(lapply(tbl.out$gene,
#                               function(TF)
#                                   nrow(subset(tbl.fimo, tf==TF & p.value <= 1e-4))))
#
# tbl.out$bs.5 <- unlist(lapply(tbl.out$gene,
#                               function(TF)
#                                   nrow(subset(tbl.fimo, tf==TF & p.value < 1e-5))))
# tbl.out$bs.6 <- unlist(lapply(tbl.out$gene,
#                               function(TF)
#                                   nrow(subset(tbl.fimo, tf==TF & p.value < 1e-6))))
#
# tbl.out.12 <- head(tbl.out, n=12)
#
# tfs <- tbl.out.12$gene
#
# for(TF in tfs){
#     tbl.fimo.sub <- subset(tbl.fimo, tf==TF)
#     track <- DataFrameAnnotationTrack(TF, tbl.fimo.sub[, c("chrom", "start", "end")], color="random",
#                                       trackHeight=25)
#     displayTrack(igv, track)
#     }
#
#
#
#
