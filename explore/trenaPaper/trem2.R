library(TrenaProjectAD)
library(RUnit)
library(futile.logger)
library(org.Hs.eg.db)
library(trenaSGM)
library(igvR)
library(MotifDb)
library(FimoClient)
library(RPostgreSQL)

source("~/github/fimoService/batchMode/fimoBatchTools.R")
#------------------------------------------------------------------------------------------------------------------------
genome <- "hg38"
targetGene <- "TREM2"
chromosome <- "chr6"
TSS <- 41163186
#------------------------------------------------------------------------------------------------------------------------
if(!exists("tp")){
   tp <- TrenaProjectAD();
   setTargetGene(tp, targetGene)
   }

if(!exists("fpdb")){
   fpdb <- dbConnect(PostgreSQL(), user= "trena", password="trena", dbname="brain_hint_20", host="khaleesi")
   }

if(!exists("fimo")){
   FIMO_HOST <- "localhost"
   FIMO_PORT <- 60000
   fimo <-FimoClient(FIMO_HOST, FIMO_PORT, quiet=FALSE)
   }

#------------------------------------------------------------------------------------------------------------------------
createGrid <- function(models, tfs.to.use)
{
   model.count <- length(models)
   all.tfs <- unique(unlist(lapply(models[1:model.count], function(x) head(x$model$gene, n=tfs.to.use))))
   model.names <- names(models)[1:model.count]
   grid <- matrix(0, nrow=length(all.tfs), ncol=model.count, dimnames=list(all.tfs, model.names))

   for(i in seq_len(model.count)){
      model.name <- names(models)[i]
      topRankedTFs <- head(models[[i]]$model$gene, n=tfs.to.use)
        #
      rank <- seq_len(tfs.to.use)   # 1, 2, 3, ... tfs.to.use
        # make high rank correspond to high numbers:
      rank.inverted <- (tfs.to.use + 1) - rank
      grid[topRankedTFs, model.name] <- rank.inverted
      }

   grid

} # createGrid
#------------------------------------------------------------------------------------------------------------------------
test_createGrid <- function()
{
   print("--- test_createGrid")

   load("models.RData")
   grid <- createGrid(models, 5)
   grid <- createGrid(models, 10)
   grid <- createGrid(models, 15)
   colnames(grid) <- c("Gene Expression Only",
                       "FIMO < 1e-4, TSS +1500bp -500bp",
                       "FIMO < 1e-3, TSS +1500bp -500bp",
                       "FIMO < 1e-4, TSS +/-5kb",
                       "FIMO < 1e-3, TSS +/-5kb",
                       "Footprints (fp), FIMO < 1e-4, +1500bp -500bp",
                       "fp, FIMO < 1e-3, +1500bp -500bp",
                       "fp, FIMO < 1e-4, +/-5kb",
                       "fp, FIMO < 1e-3, +/-5kb",
                       "fp in all promoters & enhancers, FIMO < 1e-4",
                       "fp in all promoters & enhancers, FIMO < 1e-3")

   heatmap.2(t(grid), trace="none", margins=c(20,40), col=c("#dddddd", rev(heat.colors(15))),
             density.info="none", dendrogram="none", Rowv=FALSE,
             key=FALSE,
             cexRow=2,
             cexCol=2)

} # test_createGrid
#------------------------------------------------------------------------------------------------------------------------
orig.getFootprintRegions <- function(chromosome, start, end)
{
   query <- sprintf("select * from regions where chrom='%s' and start >= %d and endpos <= %d", chromosome, start, end)
   tbl.fp.regions <- dbGetQuery(fpdb, query)
   printf("%d fp regions found", nrow(tbl.fp.regions))
   colnames(tbl.fp.regions) <- c("loc", "chrom", "start", "end")

   tbl.fp.regions <- tbl.fp.regions[, c("chrom", "start", "end")]
   new.order <- order(tbl.fp.regions$start, decreasing=FALSE)
   tbl.fp.regions <- tbl.fp.regions[new.order,]
   tbl.fp.regions <- as.data.frame(union(GRanges(tbl.fp.regions), GRanges(tbl.fp.regions)))
   colnames(tbl.fp.regions) <- c("chrom", "start", "end", "width", "strand")
   tbl.fp.regions

} # orig.getFootprintRegions
#------------------------------------------------------------------------------------------------------------------------
getFootprintRegions <- function(tbl.regions)
{
   x <- lapply(seq_len(nrow(tbl.regions)), function(i) {
      # browser()
      chromosome <- tbl.regions$chrom[i]
      start <- tbl.regions$start[i]
      end <- tbl.regions$end[i]
      query <- sprintf("select * from regions where chrom='%s' and start >= %d and endpos <= %d", chromosome, start, end)
      tbl.fp.regions <- dbGetQuery(fpdb, query)
      printf("%d fp regions found", nrow(tbl.fp.regions))
      result <- data.frame()
      if(nrow(tbl.fp.regions) > 0){
        colnames(tbl.fp.regions) <- c("loc", "chrom", "start", "end")
        tbl.fp.regions <- tbl.fp.regions[, c("chrom", "start", "end")]
        result <- tbl.fp.regions
        } # if db query > 0
      })


   nothing.found <- sum(unlist(lapply(x, length))) == 0
   printf("nothing.found? %s", nothing.found)

   if(nothing.found)
      return(data.frame())

   tbl.out <- do.call(rbind, x)
   new.order <- order(tbl.out$start, decreasing=FALSE)
   tbl.out <- tbl.out[new.order,]
   tbl.out <- as.data.frame(union(GRanges(tbl.out), GRanges(tbl.out)))
   colnames(tbl.out) <- c("chrom", "start", "end", "width", "strand")

   tbl.out

} # orig.getFootprintRegions
#------------------------------------------------------------------------------------------------------------------------
test_getFootprintRegions <- function()
{
   printf("--- test_getFootprintRegions")

   tbl.enhancers <- getEnhancers(tp)
   tbl.1 <- getFootprintRegions(tbl.enhancers[1,])
   checkEquals(dim(tbl.1), c(11, 5))
   checkEquals(colnames(tbl.1), c("chrom", "start", "end", "width", "strand"))
   checkEquals(order(tbl.1$start), seq_len(nrow(tbl.1)))

   tbl.3 <- getFootprintRegions(tbl.enhancers[c(3),])
   checkEquals(nrow(tbl.3), 25)

   tbl.13 <- getFootprintRegions(tbl.enhancers[c(1,3),])
   checkEquals(nrow(tbl.13), 36)
   checkEquals(order(tbl.13$start), seq_len(nrow(tbl.13)))

   checkEquals(nrow(rbind(tbl.1, tbl.3)), nrow(tbl.13))

      # no footprints from enhancer 2
   tbl.2 <- getFootprintRegions(tbl.enhancers[c(2),])
   checkEquals(nrow(tbl.2), 0)

   tbl.12 <- getFootprintRegions(tbl.enhancers[c(1,2),])
   checkEquals(dim(tbl.12), dim(tbl.1))

   tbl.all <- getFootprintRegions(tbl.enhancers)
   checkEquals(dim(tbl.all), c(50, 5))
   checkEquals(dim(unique(tbl.all)), c(50, 5))
   checkEquals(order(tbl.all$start), seq_len(nrow(tbl.all)))

} # test_getFootprintRegions
#------------------------------------------------------------------------------------------------------------------------
displayRegulatoryRegions <- function(tbl.regulatoryRegions)
{
   tbl.regions <- tbl.regulatoryRegions[, c("chrom", "start", "end", "combinedscore")]
   na.scores.rows <- which(is.na(tbl.regions$combinedscore))
   if(length(na.scores.rows) > 0)
      tbl.regions[na.scores.rows, "combinedscore"] <- 100
   track <- DataFrameQuantitativeTrack("gh.score", tbl.regions, color="blue", autoscale=FALSE, min=0, max=50)
   displayTrack(igv, track)

   tbl.anno <- tbl.regulatoryRegions[, c("chrom", "start", "end", "tissue")]
   track <- DataFrameAnnotationTrack("gha.regions", tbl.anno, color="blue", displayMode="EXPANDED", trackHeight=100)
   displayTrack(igv, track)


} # displayRegulatoryRegions
#------------------------------------------------------------------------------------------------------------------------
degenerateModel <- function()
{
   mtx.name <- "temporalCortex.15167x264"
   mtx <- getExpressionMatrix(tp, mtx.name)

   genome <- "hg38"
   targetGene <- "TREM2"
   setTargetGene(tp, targetGene)
   candidate.tfs <- intersect(rownames(mtx), allKnownTFs())
   checkTrue(length(candidate.tfs) > 1150)

   recipe.noDNA <- list(title="TREM2.noDNA.allTFs",
                        type="noDNA.tfsSupplied",
                        matrix=mtx,
                        candidateTFs=candidate.tfs,
                        tfPool=allKnownTFs(),
                        tfPrefilterCorrelation=0.1,
                        annotationDbFile=dbfile(org.Hs.eg.db),
                        orderModelByColumn="spearman",
                        solverNames=c("lasso", "lassopv", "ridge", "pearson", "spearman", "randomForest", "xgboost"),
                        quiet=FALSE)

   builder <- NoDnaModelBuilder(genome, targetGene, recipe.noDNA, quiet=FALSE)
   x <- build(builder)
   x

} # degenerateModel
#------------------------------------------------------------------------------------------------------------------------
fimoBatchModel <- function(tbl.regions, fimoThreshold)
{
   mtx.name <- "temporalCortex.15167x264"
   mtx <- getExpressionMatrix(tp, mtx.name)

   tbl.match <- fimoBatch(tbl.regions, fimoThreshold)
   fimo.tfs <- sort(unique(tbl.match$tf))

   genome <- "hg38"
   targetGene <- "TREM2"
   setTargetGene(tp, targetGene)
   candidate.tfs <- intersect(rownames(mtx), allKnownTFs())
   checkTrue(length(candidate.tfs) > 1150)

   recipe.noDNA <- list(title="TREM2.noDNA.allTFs",
                        type="noDNA.tfsSupplied",
                        matrix=mtx,
                        candidateTFs=fimo.tfs,
                        tfPool=allKnownTFs(),
                        tfPrefilterCorrelation=0.1,
                        annotationDbFile=dbfile(org.Hs.eg.db),
                        orderModelByColumn="spearman",
                        solverNames=c("lasso", "lassopv", "ridge", "pearson", "spearman", "randomForest", "xgboost"),
                        quiet=FALSE)

   builder <- NoDnaModelBuilder(genome, targetGene, recipe.noDNA, quiet=FALSE)
   x <- build(builder)
   tbl.match.filtered <- subset(tbl.match, tf %in% x$model$gene)
   tbl.bindingSiteCounts <- as.data.frame(table(tbl.match.filtered$tf), stringsAsFactors=FALSE)
   colnames(tbl.bindingSiteCounts) <- c("tf", "count")
   rownames(tbl.bindingSiteCounts) <- tbl.bindingSiteCounts$tf
   tbl.bindingSiteCounts <- tbl.bindingSiteCounts[,2, drop=FALSE]
   x$model$bindingSites <- tbl.bindingSiteCounts[x$model$gene, "count"]
   x$regulatoryRegions <- tbl.match.filtered
   x

} # fimoBatchModel
#------------------------------------------------------------------------------------------------------------------------
fimoRegionsModel <- function(tbl.regions, fimoThreshold=1e-4)
{
      # strand-aware start and end: trem2 is on the minus strand
      # tbl.regions <- data.frame(chrom=chromosome, start=TSS-200, end=TSS+200, stringsAsFactors=FALSE)

   mtx.name <- "temporalCortex.15167x264"
   mtx <- getExpressionMatrix(tp, mtx.name)

   build.spec <- list(title="trem2.rmm.2000up.200down",
                      type="regions.motifMatching",
                      tss=TSS,
                      regions=tbl.regions,
                      matrix=mtx,
                      pfms=query(MotifDb, "sapiens", "jaspar2018"),
                      matchThreshold=1e-4,
                      motifDiscovery="fimo",
                      tfPool=allKnownTFs(),
                      tfMapping="MotifDB",
                      annotationDbFile=dbfile(org.Hs.eg.db),
                      quiet=TRUE,
                      fimoThreshold=fimoThreshold,
                      tfPrefilterCorrelation=0.2,
                      orderModelByColumn="rfScore",
                      solverNames=c("lasso", "lassopv", "ridge", "spearman", "pearson", "randomForest", "xgboost"))

   builder <- RegionsFimoModelBuilder(genome, targetGene,  build.spec, fimo, quiet=TRUE)

   x <- build(builder)

} # fimoRegionsModel
#------------------------------------------------------------------------------------------------------------------------
displayTF <- function(tf, tbl.regulatoryRegions, tbl.fp)
{
   pfms <- query(MotifDb, c("sapiens", tf))
   pfms.jaspar2018 <- query(pfms, "jaspar2018")
   pfms.hocomoco <- query(pfms, "hocomoco")

   #if(length(pfms.jaspar2018) > 0){
   #   pfms <- pfms.jaspar2018[1]
   #} else if (length(pfms.hocomoco) > 0) {
   #   pfms <- pfms.hocomoco[1]
   #} else {
   #   pfms <- matrix()
   #}

   pfm.count <- length(pfms)
   if(pfm.count > 0){
      browser()
      for(i in seq_len(pfm.count)){
         pfm <- pfms[i]
         motifMatcher <- MotifMatcher(genomeName="hg38", pfms=as.list(pfm), quiet=TRUE)
         tbl.regions.all <- with(tbl.regulatoryRegions, data.frame(chrom=chrom[1], start=min(start), end=max(end), stringsAsFactors=FALSE))
         tbl.hits <- findMatchesByChromosomalRegion(motifMatcher, tbl.regions.all, pwmMatchMinimumAsPercentage=80)
         dim(tbl.hits)
         if(nrow(tbl.hits) == 0)
            printf("    *** no motif matches for %s", tf)
         if(nrow(tbl.hits) > 0){
            title <- sprintf("motif.%d.%s", i, tf)
            track <- DataFrameQuantitativeTrack(title, tbl.hits[, c("chrom", "motifStart", "motifEnd", "motifRelativeScore")],
                                                color="darkgreen", autoscale=TRUE)
            displayTrack(igv, track)
            } # if tbl.hits
         } # for pfm
      } # if pfms

   setTargetGene(tp, tf)
   tbl.chipSeq <- getChipSeq(tp, chrom="chr6", start=min(tbl.regulatoryRegions$start) - shoulder,
                             end=max(tbl.regulatoryRegions$end + shoulder), tfs=tf)
   if(nrow(tbl.chipSeq) > 0){
      title <- sprintf("ChIP.%s", tf)
      track <- DataFrameAnnotationTrack(title, tbl.chipSeq[, c("chrom", "peakStart", "peakEnd", "name")], color="orange")
      displayTrack(igv, track)
      }

   tbl.fp.filtered <- tbl.fp[grep(tf, tbl.fp$motifName),]
   tbl.fp.fixed <- tbl.fp.filtered[, c("fp_start", "fp_end", "score1")]
   tbl.fp.fixed$chrom <- rep(tbl.regulatoryRegions$chrom[1], nrow(tbl.fp.fixed))
   colnames(tbl.fp.fixed) <- c("start", "end", "score", "chrom")
   tbl.fp.fixed <- tbl.fp.fixed[, c("chrom", "start", "end", "score")]
   browser()

   title <- sprintf("fp.%s", tf)
   track <- DataFrameAnnotationTrack(title, tbl.fp.fixed, color="red")
   displayTrack(igv, track)

} # displayTF
#------------------------------------------------------------------------------------------------------------------------
# intended for interactive line-by-line exectuion
explore_TREM2 <- function()
{
   genome <- "hg38"
   targetGene <- "TREM2"
   setTargetGene(tp, targetGene)
   tbl.info <- getTranscriptsTable(tp)
   tss <- tbl.info$tss
   chrom <- tbl.info$chrom

      # get genehancer regulatory regions, all tissues, then demonstrate
      # how further filtering can be done if you wish to do so

   tbl.regulatoryRegions <- getGeneRegulatoryRegions(tp, tissues="all",
                                                     geneHancerSupplemental.promoter.upstream=2000,
                                                     geneHancerSupplemental.promoter.downstream=1000)
   checkTrue(nrow(tbl.regulatoryRegions) > 8)
   supplemental.row <- grep("TrenaProjectHG38", tbl.regulatoryRegions$source)
   stopifnot(length(supplemental.row) == 1)
   tbl.regulatoryRegions$combinedScore[supplemental.row] <- 500
   tbl.regulatoryRegions$tissue[supplemental.row] <- "supplemental"
   getExpressionMatrixNames(tp)
   mtx <- getExpressionMatrix(tp, "temporalCortex.15167x264")

      #----------------------------------------------------------------------------------------------------
      # first, build a model with "placenta2", an early version of the placenta footprint database
      #----------------------------------------------------------------------------------------------------


   browser()
   tbl.all <- with(tbl.regulatoryRegions, data.frame(chrom=chrom[1], start=min(start), end=max(end), stringsAsFactors=FALSE))

   recipe <- list(title="gh all tissue demo on TREM2",
                  type="footprint.database",
                  regions=tbl.all,
                  #regions=tbl.regulatoryRegions,
                  geneSymbol=targetGene,
                  tss=tss,
                  matrix=mtx,
                  db.host="khaleesi.systemsbiology.net",
                  db.port=5432,
                  databases=list("brain_hint_20"),
                  annotationDbFile=dbfile(org.Hs.eg.db),
                  motifDiscovery="builtinFimo",
                  tfPool=allKnownTFs(),
                  tfMapping="MotifDB",
                  tfPrefilterCorrelation=0.1,
                  orderModelByColumn="rfScore",
                  solverNames=c("lasso", "pearson", "randomForest", "ridge", "spearman", "xgboost"))

   fpBuilder <- FootprintDatabaseModelBuilder(genome, targetGene, recipe, quiet=FALSE)
   x.trem2 <- build(fpBuilder)

   tbl.model <- x.trem2$model[1:20,]
   dim(tbl.model)
   tbl.fp.all <- x.trem2$regulatoryRegions
   tbl.fp <- subset(tbl.fp.all, geneSymbol %in% tbl.model$gene)
   dim(tbl.fp)


   igv <- igvR()
   setGenome(igv, "hg38")
   shoulder <- 1000
   loc <- sprintf("%s:%d-%d", tbl.info$chrom, min(tbl.regulatoryRegions$start) - shoulder,
                  max(tbl.regulatoryRegions$end + shoulder))
   showGenomicRegion(igv, loc)
   displayRegulatoryRegions(tbl.regulatoryRegions)

   for(tf in tbl.model$gene[1:3]){
      displayTF(tf, tbl.regulatoryRegions, tbl.fp)
      }

} # explore_TREM2
#------------------------------------------------------------------------------------------------------------------------
buildFimoModels <- function()
{
   x.0 <- degenerateModel()
   lapply(x.0, dim)

   tbl.regions.1 <- data.frame(chrom=chromosome, start=TSS-500, end=TSS+1500, stringsAsFactors=FALSE)

   #x.1 <- fimoRegionsModel(tbl.regions.1, 1e-4)
   #lapply(x.1, dim)
   x.1 <- fimoBatchModel(tbl.regions.1, 1e-4)
   lapply(x.1, dim)

   x.2 <- fimoBatchModel(tbl.regions.1, 1e-3)
   lapply(x.2, dim)

   tbl.regions.2 <- data.frame(chrom=chromosome, start=TSS-5000, end=TSS+5000, stringsAsFactors=FALSE)
   x.3 <- fimoBatchModel(tbl.regions.2, 1e-4)
   lapply(x.3, dim)

   x.4 <- fimoBatchModel(tbl.regions.2, 1e-3)
   lapply(x.4, dim)

   tbl.fp.1 <- getFootprintRegions(tbl.regions.1) # ith(tbl.regions.1, getFootprintRegions(chrom, start, end))
   dim(tbl.fp.1)
   x.5 <- fimoBatchModel(tbl.fp.1, 1e-4)
   lapply(x.5, dim)

   x.6 <- fimoBatchModel(tbl.fp.1, 1e-3)
   lapply(x.6, dim)

   tbl.fp.2 <- getFootprintRegions(tbl.regions.2)
   dim(tbl.fp.2)
   x.7 <- fimoBatchModel(tbl.fp.2, 1e-4)
   lapply(x.7, dim)

   x.8 <- fimoBatchModel(tbl.fp.2, 1e-3)
   lapply(x.8, dim)

      # now bring in genehancer.  first, use the "all tissues" default
   tbl.enhancers <- getEnhancers(tp)
   tbl.fp.gh.all <- getFootprintRegions(tbl.enhancers)
   dim(tbl.fp.gh.all)
   x.9 <- fimoBatchModel(tbl.fp.gh.all, 1e-4)
   lapply(x.9, dim)

   x.10 <- fimoBatchModel(tbl.fp.gh.all, 1e-3)
   lapply(x.10, dim)


   models <- list(noGenome=x.0,
                  m2kb.f4=x.1,
                  m2kb.f3=x.2,
                  m10kb.f4=x.3,
                  m10kb.f3=x.4,
                  m2kb.fp.f4=x.5,
                  m2kb.fp.f3=x.6,
                  m10kb.fp.f4=x.7,
                  m10kb.fp.f3=x.8,
                  mgh.fp.f4=x.9,
                  mgh.fp.f3=x.10)
   save(models, file="models.RData")



   tfs.degenerate <- head(x.0$model$gene, n=15)
   tfs.promoter.short.fimo4 <- head(x.1$model$gene, n=15)
   tfs.promoter.short.fimo3 <- head(x.2$model$gene, n=15)

   tfs.promoter.10k.fimo4 <- head(x.3$model$gene, n=15)
   tfs.promoter.10k.fimo3 <- head(x.4$model$gene, n=15)

   tfs.promoter.fp.short.fimo4 <- head(x.5$model$gene, n=15)
   tfs.promoter.fp.short.fimo3 <- head(x.6$model$gene, n=15)
   tfs.promoter.fp.10k.fimo4 <- head(x.7$model$gene, n=15)
   tfs.promoter.fp.10k.fimo3 <- head(x.8$model$gene, n=15)

   tfs.enhancers.fp.fimo4 <- head(x.9$model$gene, n=15)
   tfs.enhancers.fp.fimo3 <- head(x.10$model$gene, n=15)

   all.tf.sets <- list(tf0=tfs.degenerate,
                       tf1=tfs.promoter.short.fimo4,
                       tf2=tfs.promoter.short.fimo3,
                       tf3=tfs.promoter.10k.fimo4,
                       tf4=tfs.promoter.10k.fimo3,
                       tf5=tfs.promoter.fp.short.fimo4,
                       tf6=tfs.promoter.fp.short.fimo3,
                       tf7=tfs.promoter.fp.10k.fimo4,
                       tf8=tfs.promoter.fp.10k.fimo3,
                       tf9=tfs.enhancers.fp.fimo4,
                       tf10=tfs.enhancers.fp.fimo3)
   model.names <-list(tf0="noGenome",
                       tf1="2kb,f4",
                       tf2="2kb.f3",
                       tf3="10kb.f4",
                       tf4="10kb.f3",
                       tf5="2kb,fp,f4",
                       tf6="2kb,fp,f3",
                       tf7="10kb,fp,f4",
                       tf8="10kb,fp,f3",
                       tf9="gh,fp,f4",
                       tf10="gh,fp,f3")

   tfs.all <- unique(Reduce(c, all.tf.sets))
   tfs.shared <- Reduce(intersect, all.tf.sets)

   save(all.tf.sets, model.names, tfs.all, tfs.shared, file="11model-summary.RData")


   list(tfs.degenerate,
                                        tfs.promoter.short.fimo4,
                                        tfs.promoter.short.fimo3,
                                        tfs.promoter.10k.fimo4,
                                        tfs.promoter.10k.fimo3,
                                        tfs.promoter.fp.short.fimo4,
                                        tfs.promoter.fp.short.fimo3,
                                        tfs.promoter.fp.10k.fimo4,
                                        tfs.promoter.fp.10k.fimo3,
                                        tfs.enhancers.fp.fimo4,
                                        tfs.enhancers.fp.fimo3))

   a <- tfs.degenerate
   b <- tfs.enhancers.fp.fimo3
   c <- tfs.promoter.short.fimo3

   draw.triple.venn(length(a),
                    length(b),
                    length(c),
                   length(intersect(a, b)),
                   length(intersect(b, c)),
                   length(intersect(a, c)),
                   length(intersect(a, intersect(b, c))),
                   category=c("degenerate", "enhancers.3", "promoter.short.4"),
                   col=c("red", "green", "blue"))




     # genehancer regions for trem2 restricted to brain-related tissues - nothing there for TREM2.  odd
   tissues.of.interest <- grep("brain", getEnhancerTissues(tp), ignore.case=TRUE, value=TRUE)
   tbl.enhancers <- getEnhancers(tp, tissues=tissues.of.interest)
   dim(tbl.enhancers)  # 0 0
   tissues.of.interest <- grep("neural", getEnhancerTissues(tp), ignore.case=TRUE, value=TRUE)
   tbl.enhancers <- getEnhancers(tp, tissues=tissues.of.interest)
   dim(tbl.enhancers) # 0 0

   #tbl.fp.gh.brain <- getFootprintRegions(tbl.enhancers)
   #dim(tbl.fp.gh.brain)
   #x.9 <- fimoBatchModel(tbl.fp.gh.brain, 1e-4)
   #lapply(x.10, dim)

   #x.10 <- fimoBatchModel(tbl.fp.gh.brain, 1e-3)
   #lapply(x.10, dim)

} # buildFimoModels
#------------------------------------------------------------------------------------------------------------------------
