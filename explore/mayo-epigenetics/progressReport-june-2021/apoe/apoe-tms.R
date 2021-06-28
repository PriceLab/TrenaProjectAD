source("../tmsCore.R")
tp <- TrenaProjectAD()
targetGene <- "APOE"
mtx.rna <- get(load("~/github/TrenaProjectAD/explore/mayo-epigenetics/progressReport-june-2021/mayo.tcx.16969x262.covariateCorrection.log+scale.RData"))

tms <- TMS$new("APOE",
               tp,

