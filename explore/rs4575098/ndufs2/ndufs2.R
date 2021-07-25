# ndufs2.R
#----------------------------------------------------------------------------------------------------
source("../tmsRunner.R")
source("../determineRegionForModel.R")

targetGene <- "NDUFS2"
mtx.rna <- get(load("~/github/TrenaProjectAD/inst/extdata/expression/mayo.tcx.16969x262.covariateCorrection.log+scale.RData"))

# tbl.region <- determineRegionForModel(targetGene)

library(ghdb)
ghdb <- GeneHancerDB()
tbl.gh <- retrieveEnhancersFromDatabase(ghdb, targetGene, tissues="all")

chromosome <- tbl.gh$chrom[1]
start <- min(tbl.gh$start) - 1000
end <- max(tbl.gh$end) + 1000

chromosome <- "chr1"
start <- 161184428
end <-  161204399

# promoter + open atac
chromosome <- "chr1"
start <- 161201350
end <- 161203624

printf("full region for motif search: %dkb", as.integer((end-start)/1000))

x <- tryCatch(
    runTMS(targetGene, chromosome, start, end, mtx.rna, fimo.threshold=1e-6),
    error=function(e){
        printf("failure with %s", targetGene)
        print(e)
        return(NA)
        }
    )
if(!is.na(x)){
   lm.tables <- x$get.lm.tables()
   lm.rsquareds <- x$get.lm.Rsquareds()
   filename <- sprintf("%s-model-%s:%d-%d.RData", targetGene, chromosome, start, end)
   save(x, file=filename)
   lm.tables <- x$get.lm.tables()
   lapply(names(lm.tables), function(name) {
       printf("-------- %s", name)
       print(subset(lm.tables[[name]])) # , p.value < 0.25))
       })
   print(x$get.lm.Rsquareds())
  } # if proceed

