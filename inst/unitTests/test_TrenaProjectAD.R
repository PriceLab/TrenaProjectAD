library(TrenaProjectAD)
library(RUnit)
library(futile.logger)
library(org.Hs.eg.db)
library(trenaSGM)
#------------------------------------------------------------------------------------------------------------------------
flog.appender(appender.file("timing-TreneProjectAD.log"))

logTimingInfo <- function(msg, timingInfo){
   string <- sprintf("%50s: %5.2f %5.2f %5.2f %5.2f %5.2f", msg,
                     timingInfo[[1]], timingInfo[[2]], timingInfo[[3]], timingInfo[[4]], timingInfo[[5]])
   flog.info(string, "timingLog")
    }
#------------------------------------------------------------------------------------------------------------------------
if(!exists("tp"))
   tp <- TrenaProjectAD();
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_basic()
   test_buildSingleGeneModel_withNewGeneHancerClass()
   test_buildSingleGeneModel_noDNA()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
test_basic <- function()
{
   message(sprintf("---- test_basic"))

   expected <-  c("ABCA7", "APOE", "BIN1", "CASS4", "CD2AP", "CD33", "CELF1", "CLU", "CR1", "DSG2", "EPHA1",
                  "FERMT2", "HLA-DRB1", "HLA-DRB5", "INPP5D", "MEF2C", "MS4A6A", "NME8", "PICALM", "PTK2B",
                  "RIN3", "SLC24A4", "SORL1", "TOMM40", "ZCWPW1", "TREM2")

   checkTrue(all(expected %in% getSupportedGenes(tp)))

   expected <- c("brain_wellington_16", "brain_wellington_20", "brain_hint_16", "brain_hint_20")
   checkTrue(all(expected %in% getFootprintDatabaseNames(tp)))
   checkEquals(getFootprintDatabaseHost(tp), "khaleesi.systemsbiology.net")

   expected <- c("cerebellum.15167x263", "rosmap.14235x632", "temporalCortex.15167x264")
   checkTrue(all(expected %in% getExpressionMatrixNames(tp)))
   mtx <- getExpressionMatrix(tp, "temporalCortex.15167x264")
   checkEquals(dim(mtx), c(15167, 264))

      # note this next test only works wtih tcx.
      # TODO: get covariate data from Mayo cerebellum and rosmap expression matrices

   setTargetGene(tp, "TREM2")
   checkEquals(getTargetGene(tp), "TREM2")

   tbl.transcripts <- getTranscriptsTable(tp)
   checkTrue(nrow(tbl.transcripts) == 1)

   tbl.enhancers <- getGeneRegulatoryRegions(tp)
   dim(tbl.enhancers)
   checkEquals(head(colnames(tbl.enhancers)), c("chrom", "start", "end", "gene", "eqtl", "hic"))
   checkTrue(nrow(tbl.enhancers) > 8)

   tbl.dhs <- getEncodeDHS(tp)
   checkTrue(nrow(tbl.dhs) > 100)
   checkEquals(colnames(tbl.dhs), c("chrom", "chromStart", "chromEnd", "count", "score"))

   chromosome <- unique(c(tbl.dhs$chrom, tbl.enhancers$chrom))
   checkEquals(length(chromosome), 1)
   start <- min(tbl.dhs$chromStart)
   end   <- max(tbl.dhs$chromEnd)

   tbl.chipSeq <- getChipSeq(tp, chrom=chromosome, start=start, end=end, tfs=NA)
   dim(tbl.chipSeq)
   checkTrue(nrow(tbl.chipSeq) > 4000)
   checkEquals(colnames(tbl.chipSeq), c("chrom", "start", "endpos", "tf", "name", "strand", "peakStart", "peakEnd"))

} # test_basic
#------------------------------------------------------------------------------------------------------------------------
test_buildSingleGeneModel_withNewGeneHancerClass <- function()
{
   printf("--- test_buildSingleGeneModel_withNewGeneHancerClass")

   genome <- "hg38"
   targetGene <- "TREM2"
   setTargetGene(tp, targetGene)
   tbl.info <- getTranscriptsTable(tp)
   tss <- tbl.info$tss
   chrom <- tbl.info$chrom

      # get genehancer regulatory regions, all tissues, then demonstrate
      # how further filtering can be done if you wish to do so

   logTimingInfo("gh tissues all, TREM2",  system.time(tbl.regulatoryRegions <- getGeneRegulatoryRegions(tp, tissues="all")))
   checkTrue(nrow(tbl.regulatoryRegions) > 8)

   getExpressionMatrixNames(tp)
   mtx <- getExpressionMatrix(tp, "temporalCortex.15167x264")

      #----------------------------------------------------------------------------------------------------
      # first, build a model with "placenta2", an early version of the placenta footprint database
      #----------------------------------------------------------------------------------------------------


   recipe <- list(title="gh all tissue demo on TREM2",
                  type="footprint.database",
                  regions=tbl.regulatoryRegions,
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
                  solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))

   fpBuilder <- FootprintDatabaseModelBuilder(genome, targetGene, recipe, quiet=FALSE)

   msg <- sprintf("gh all tissues demo on TREM2, fp model %d regions", nrow(tbl.regulatoryRegions))
   logTimingInfo(msg, system.time(x.trem2 <- build(fpBuilder)))

   checkTrue("SPI1" %in% head(x.trem2$model$gene))

} # test_buildSingleGeneModel_withNewGeneHancerClass
#------------------------------------------------------------------------------------------------------------------------
test_buildSingleGeneModel_noDNA <- function(reportTiming=FALSE)
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
                        tfPrefilterCorrelation=0.5,
                        annotationDbFile=dbfile(org.Hs.eg.db),
                        orderModelByColumn="rfScore",
                        solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"),
                        quiet=FALSE)

   builder <- NoDnaModelBuilder(genome, targetGene, recipe.noDNA, quiet=FALSE)
   logTimingInfo("building noDNA model, cor prefilter 0.3", system.time(x.noDNA <- build(builder)))

   checkEquals(x.noDNA$regulatoryRegions, data.frame())
   tbl.model <- head(x.noDNA$model, n=10)
   checkTrue(all(abs(tbl.model$pearsonCoeff) > 0.5))
   checkTrue("SPI1" %in% tbl.model$gene)

} # test_buildSingleGeneModel_noDNA
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()
