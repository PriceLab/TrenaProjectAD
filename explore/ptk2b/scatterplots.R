#----------------------------------------------------------------------------------------------------
load.matrices <- function()
{
   dir <- "~/github/TrenaProjectAD/prep/rna-seq-counts-from-synapse/eqtl"
   file.tcx <- "mtx.mayo.tcx.eqtl-optimized-geneSymbols-sampleIDs-with-vcf17009x257.RData"
   file.cer <- "mtx.mayo.cer.eqtl-optimized-geneSymbols-sampleIDs-with-vcf17009x255.RData"
   file.exists(file.path(dir, file.cer))

   mtx.tcx <- get(load(file.path(dir, file.tcx)))
   mtx.cer <- get(load(file.path(dir, file.cer)))
   mtx.rosmap <- get(load(file.path(dir,  "mtx.rosmap.rnaseq-counts-geneSymbols-15582x632.RData")))
   mtx.sinai  <- get(load(file.path(dir, "mtx.mssm.rnaseq-residual-eqtl-geneSymbols-patients-16346x753.RData")))
   mtx.sinai[which(is.na(mtx.sinai))] <- 0
   mtx.tcx[which(is.na(mtx.tcx))] <- 0
   mtx.cer[which(is.na(mtx.cer))] <- 0

   dir <- "~/github/TrenaProjectAD/inst/extdata/expression"
   file.covariateCorrected <- "mayo.tcx.16969x262.covariateCorrection.log+scale.RData"
   file.windsorized <- "Scaled_Winsorized_MayoRNAseq_TCX-ENSG.RData"
   mtx.old.tcx.cov <- get(load(file.path(dir, file.covariateCorrected)))
   mtx.old.tcx.win <- get(load(file.path(dir, file.windsorized)))
   mtx.old.rosmap <- get(load("~/github/TrenaProjectAD/inst/extdata/expression/rosmap.14235x632.RData"))



} # load.matrices
#----------------------------------------------------------------------------------------------------
load.matrices()
plot(mtx.tcx["EGR4",], mtx.tcx["PTK2B",], col="blue", ylim=c(-11,11), xlim=c(-11, 11),
     main="EGR4 vs PTK2B expression in 6 datasets", xlab="EGR4", ylab="PTK2B")
points(mtx.cer["EGR4",], mtx.cer["PTK2B",], col="green")
points(mtx.rosmap["EGR4",], mtx.rosmap["PTK2B",], col="red")
points(mtx.sinai["EGR4",], mtx.sinai["PTK2B",], col="orange")
points(mtx.old.rosmap["EGR4",], mtx.old.rosmap["PTK2B",], col="pink")
points(mtx.old.tcx.cov["EGR4",], mtx.old.tcx.cov["PTK2B",], col="black")
legend(5, 9,
       c("mayo eqtl tcx", "mayo eqtl cer", "rosmap eqtl", "sinai eqtl", "old rosmap", "old mayo tcx"),
       c("blue", "green", "red", "orange", "pink", "black"))

# saved plot as ptk2b-egr4-scatter-plot-6-datasets.png
