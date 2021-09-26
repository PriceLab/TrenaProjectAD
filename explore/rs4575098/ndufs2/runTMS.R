library(RUnit)
source("~/github/TrenaMultiScore/tools/runner/v2/tmsCore.R")
targetGene <- "NDUFS2"
tbl.fimo <- get(load("tbl.fimo.NDUFS2.RData"))
#----------------------------------------------------------------------------------------------------
test_inclusiveFimoTable <- function()
{
   message(sprintf("--- test_inclusiveFimoTable"))

   tms <- TMS$new(trenaProject, targetGene, tbl.fimo, tbl.atac, quiet=FALSE)
   tms$scoreFimoTFBS()   # chip, conservation, genehancer, genic annotations, distance to tss
   tms$add.tf.mrna.correlations(mtx.rna, featureName="cor.all")

   tbl.tms <- tms$getTfTable()
   checkEquals(nrow(tbl.tms), nrow(tbl.fimo))

   dim(subset(tbl.tms, fimo_pvalue <= 1e-6 & chip))

   dim(subset(tbl.tms, fimo_pvalue < 1e-5 & abs(cor.all) > 0.8 & chip & gh > 600))

   tfs <- subset(tbl.tms, fimo_pvalue < 1e-5 & abs(cor.all) > 0.7 & gh > 1)$tf
   length(tfs)

   tms$addRBP()
   tms$add.rbp.mrna.correlations(mtx.rna, featureName="cor.all")   # added to tbl.rbp
   tbl.rbp <- tms$getRbpTable()
   dim(tbl.rbp)

   rbps <- unique(subset(tbl.rbp, abs(cor.all) > 0.3)$gene)
   printf("candidate rbps: %d", length(rbps))
   tbl.trena.tf <- tms$build.trena.model(tfs, list(), mtx.rna)

   checkEquals(unique(tbl.trena.tf$class), "tf")
   checkTrue(all(c("NR3C1", "ELF3", "MAFG", "BCL6", "KLF13", "IRF1") %in% tbl.trena.tf$gene[1:10]))

   tbl.trena.both <- tms$build.trena.model(tfs, rbps, mtx.rna)

   checkTrue(all(c("tf", "rbp") %in% tbl.trena.both$class))
   checkTrue(all(c("STAU1", "NR3C1", "ELF3", "RBM15", "MAFG", "BCL6") %in% tbl.trena.both$gene[1:10]))

} # test_inclusiveFimoTable
#----------------------------------------------------------------------------------------------------
