library(trena)

if(!exists("tbl.tms")){
    load("tbl.tms.RData")
    mtx.all    <- x$all
    mtx.rosmap <- x$rosmap
    mtx.sinai  <- x$sinai
    mtx.mayo   <- x$mayo
    }

#------------------------------------------------------------------------------------------------------------------------
split.rna.matrix.by.rsid <- function(rsid="rs28834970")
{
   etx <- EndophenotypeExplorer$new(targetGene, "hg38")
   tbl.map <- etx$getIdMap()
   dim(tbl.map)  # 4026 4

   mtx.geno <- etx$getGenoMatrixByRSID(rsid)
   dim(mtx.geno)
   table(mtx.geno)     # 794  861 239
   samples.hom <- names(which(mtx.geno[1,] == "1/1"))
   samples.het <- names(which(mtx.geno[1,] == "0/1"))
   samples.wt  <- names(which(mtx.geno[1,] == "0/0"))
   stopifnot(length(unique(c(samples.hom, samples.het, samples.wt))) == ncol(mtx.geno))

   length(samples.hom)  # 239
   length(samples.het)  # 861
   length(samples.wt)   # 794

       # todo: use this only 137/235 to track down mapping failure
       # failures.hom <- head(setdiff(patients.hom, colnames(mtx.all)))
       #    "14556" "18198" "18200" "18212" "18224" "18228"
       #     sample patient study assay
       #  205  14556   14556  mayo   vcf
       #  211  18198   18198  mayo   vcf
       #  213  18200   18200  mayo   vcf
       #  222  18212   18212  mayo   vcf
       #  233  18224   18224  mayo   vcf
       #  236  18228   18228  mayo   vcf

   patients.hom <- subset(tbl.map, sample %in% samples.hom)$patient
   length(patients.hom)   # 235
   length(intersect(patients.hom, colnames(mtx.all))) # [1] 137

   patients.het <- subset(tbl.map, sample %in% samples.het)$patient
   length(patients.het)   # 834
   length(intersect(patients.het, colnames(mtx.all))) # [1] 489

   patients.wt <- subset(tbl.map, sample %in% samples.wt)$patient
   length(patients.wt)   # 776
   length(intersect(patients.wt, colnames(mtx.all))) # [1] 453

   patients.hom.in.mtx <- intersect(patients.hom, colnames(mtx.all))
   length(patients.hom)
   length(patients.hom.in.mtx)

   patients.het.in.mtx <- intersect(patients.het, colnames(mtx.all))
   length(patients.het)
   length(patients.het.in.mtx)

   patients.wt.in.mtx <- intersect(patients.wt, colnames(mtx.all))
   length(patients.wt)
   length(patients.wt.in.mtx)

   mtx.all.hom <- mtx.all[, patients.hom.in.mtx]
   mtx.all.het <- mtx.all[, patients.het.in.mtx]
   mtx.all.wt <- mtx.all[, patients.wt.in.mtx]
   mtx.all.mut <- mtx.all[, c(patients.hom.in.mtx, patients.het.in.mtx)]
   return(list(het=mtx.all.het, hom=mtx.all.hom, mut=mtx.all.mut, wt=mtx.all.wt))

} # split.rna.matrix.by.rsid
#------------------------------------------------------------------------------------------------------------------------
run.wt.trena <- function(mtx.rna, tbl.tms.sub, extra.tfs=c())
{
    tfs <- unique(tbl.tms.sub$tf)
    solver <- EnsembleSolver(mtx.rna,
                             targetGene=targetGene,
                             candidateRegulators=c(tfs, extra.tfs),
                             solverNames=c("lasso", "Ridge", "Spearman", "Pearson", "RandomForest", "xgboost"),
                             geneCutoff=0.9)
    tbl.trena <- run(solver)
    tbl.trena <- tbl.trena[order(tbl.trena$rfScore, decreasing=TRUE),]
    tbl.trena$tfbs <- unlist(lapply(tbl.trena$gene, function(gene) nrow(subset(tbl.tms.sub, tf==gene))))
    rownames(tbl.trena) <- NULL
    numeric.colnames <- names(which(lapply(tbl.trena, class) == "numeric"))
    for(colname in numeric.colnames)
        tbl.trena[,colname] <- round(tbl.trena[,colname], digits=3)
    tbl.trena

} # run.wt.trena
#------------------------------------------------------------------------------------------------------------------------
run.trenas <- function()
{
    tbl.tms.sub <- subset(tbl.tms, (oc | chip) & gh > 600)
    tbl.tms.sub <- subset(tbl.tms, (oc | chip) & gh > 600 & fimo_pvalue < 0.001)
    tfs <- unique(tbl.tms.sub$tf)
    "RARA" %in% tfs
    "CUX2" %in% tfs
    length(tfs)  #483

    models <- list()
    models[["all"]] <- run.wt.trena(mtx.all, tbl.tms.sub)

    mtx.div <- split.rna.matrix.by.rsid()
    names(mtx.div)
    lapply(mtx.div, dim)

    models[["mut"]] <- run.wt.trena(mtx.div$mut, tbl.tms.sub)
    models[["wt"]] <- run.wt.trena(mtx.div$wt, tbl.tms.sub)

} # run.trenas
#----------------------------------------------------------------------------------------------------
