library(trenaSGM)
library(TrenaProjectAD)
library(org.Hs.eg.db)
library(RUnit)

genome <- "hg38"
targetGene <- "TREM2"

tp <- TrenaProjectAD();

setTargetGene(tp, targetGene)
tbl.info <- getTranscriptsTable(tp)

chromosome <- tbl.info$chrom
tss <- tbl.info$tss

tbl.gh <- getEnhancers(tp)[, c("chrom", "start", "end")]

# by default, merge genehancers with tss +/- 2k
# if no genehancers, set to tss +/- 5k

#shoulder <- 2000
#if(nrow(tbl.gh) == 0)
shoulder <- 5000

tbl.proxPromoter <- data.frame(chrom=chromosome, start=tss-shoulder, end=tss+shoulder, stringsAsFactors=FALSE)

if(nrow(tbl.gh) == 0){
   tbl.regions <- tbl.proxPromoter
}else{
   gr.gh <- GRanges(tbl.gh)
   gr.pp <- GRanges(tbl.proxPromoter)
   gr.both <- sort(c(gr.gh, gr.pp))
   tbl.regions <- as.data.frame(reduce(gr.both))
   colnames(tbl.regions)[1] <- "chrom"
   tbl.regions$chrom <- as.character(tbl.regions$chrom)
   }

mtx.1 <- getExpressionMatrix(tp, "corrected_PMI_AgeAtDeath_astrocyte_excitatory_inhibitory_microglia_oligodendrocyte")
mtx.2 <- getExpressionMatrix(tp, "mayo.tcx.16969x262.covariateCorrection.log+scale")
mtx.3 <- getExpressionMatrix(tp, "temporalCortex.15167x264")
fpdbs <- c("brain_hint_16", "brain_hint_20", "brain_wellington_16", "brain_wellington_20")[c(2,4)]

build.spec <- list(title="unit test on TREM2",
                   type="footprint.database",
                   regions=tbl.regions,
                   geneSymbol=targetGene,
                   tss=tss,
                   matrix=mtx.1,
                   db.host="khaleesi.systemsbiology.net",
                   db.port=5432,
                   databases=fpdbs,
                   annotationDbFile=dbfile(org.Hs.eg.db),
                   motifDiscovery="builtinFimo",
                   tfPool=allKnownTFs(),
                   tfMapping="MotifDB",
                   tfPrefilterCorrelation=0.15,
                   solverNames=c("spearman", "pearson"),
                   orderModelByColumn="spearmanCoeff"
                   )

fpBuilder <- FootprintDatabaseModelBuilder(genome, targetGene,  build.spec, quiet=TRUE)
suppressWarnings(x <- build(fpBuilder))
lapply(x, dim)
tbl.model <- x$model
tbl.reg   <- x$regulatoryRegions

checkTrue(all(c("IKZF1","SPI1","IRF8","RUNX1","NFATC2","FLI1","TAL1","RUNX3","CEBPA","ELK3") %in%
              tbl.model$gene[1:20]))
checkTrue(all(tbl.model$spearman[1:10] > 0.5))
