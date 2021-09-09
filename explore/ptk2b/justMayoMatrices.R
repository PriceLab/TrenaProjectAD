library(trena)
library(EndophenotypeExplorer)
library(RUnit)
targetGene <- "PTK2B"
#----------------------------------------------------------------------------------------------------
#load.matrices <- function()
#{

if(!exists("mtx.tcx")){
   dir <- "~/github/TrenaProjectAD/prep/rna-seq-counts-from-synapse/eqtl"
   file.tcx <- "mtx.mayo.tcx.eqtl-optimized-geneSymbols-sampleIDs-with-vcf17009x257.RData"
   file.cer <- "mtx.mayo.cer.eqtl-optimized-geneSymbols-sampleIDs-with-vcf17009x255.RData"
   file.exists(file.path(dir, file.cer))

   mtx.tcx <- get(load(file.path(dir, file.tcx)))
   mtx.cer <- get(load(file.path(dir, file.cer)))

   mtx.tcx[which(is.na(mtx.tcx))] <- 0
   mtx.cer[which(is.na(mtx.cer))] <- 0

   dir <- "~/github/TrenaProjectAD/inst/extdata/expression"
   file.covariateCorrected <- "mayo.tcx.16969x262.covariateCorrection.log+scale.RData"
   file.windsorized <- "Scaled_Winsorized_MayoRNAseq_TCX-ENSG.RData"
   mtx.old.tcx.cov <- get(load(file.path(dir, file.covariateCorrected)))
   mtx.old.tcx.win <- get(load(file.path(dir, file.windsorized)))
   mtx.rosmap <- get(load("~/github/TrenaProjectAD/inst/extdata/expression/rosmap.14235x632.RData"))
   }

#----------------------------------------------------------------------------------------------------
split.rna.matrix.by.rsid <- function(mtx, rsid="rs28834970", suffix="_CER")
{
   etx <- EndophenotypeExplorer$new(targetGene, "hg38")
   tbl.map <- etx$getIdMap()
   dim(tbl.map)  # 4026 4

   mtx.geno <- etx$getGenoMatrixByRSID(rsid)
   dim(mtx.geno)
   table(mtx.geno)     # 794  861 239   for rs28834970

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


   patients.hom <- subset(tbl.map, sample %in% samples.hom & study=="mayo" & assay=="vcf")$patient
   length(patients.hom)   # 53
   rna.samples.hom <- subset(tbl.map, patient %in% patients.hom & study=="mayo" & assay=="rnaseq")$sample
   rna.samples.hom.tcx <- grep(suffix, rna.samples.hom, value=TRUE)
   rna.samples.hom.tcx <- sub(suffix, "", rna.samples.hom.tcx)
   length(rna.samples.hom.tcx)   # 41
   length(intersect(rna.samples.hom.tcx, colnames(mtx))) # 41

   patients.het <- subset(tbl.map, sample %in% samples.het & study=="mayo" & assay=="vcf")$patient
   length(patients.het)   # 167
   rna.samples.het <- subset(tbl.map, patient %in% patients.het & study=="mayo" & assay=="rnaseq")$sample
   rna.samples.het.tcx <- grep(suffix, rna.samples.het, value=TRUE)
   rna.samples.het.tcx <- sub(suffix, "", rna.samples.het.tcx)
   length(rna.samples.het.tcx)   # 122
   length(intersect(rna.samples.het.tcx, colnames(mtx))) # 122

   patients.wt <- subset(tbl.map, sample %in% samples.wt & study=="mayo" & assay=="vcf")$patient
   length(patients.wt)   # 129
   rna.samples.wt <- subset(tbl.map, patient %in% patients.wt & study=="mayo" & assay=="rnaseq")$sample
   rna.samples.wt.tcx <- grep(suffix, rna.samples.wt, value=TRUE)
   rna.samples.wt.tcx <- sub(suffix, "", rna.samples.wt.tcx)
   length(rna.samples.wt.tcx)   # 94
   length(intersect(rna.samples.wt.tcx, colnames(mtx))) # 94


   mtx.hom <- mtx[, rna.samples.hom.tcx]
   mtx.het <- mtx[, rna.samples.het.tcx]
   mtx.wt <-  mtx[, rna.samples.wt.tcx]
   mtx.mut <- mtx[, c(rna.samples.hom.tcx, rna.samples.het.tcx)]
   return(list(het=mtx.het, hom=mtx.hom, mut=mtx.mut, wt=mtx.wt))

} # split.rna.matrix.by.rsid
#----------------------------------------------------------------------------------------------------
test_split.rna.matrix.by.rsid <- function()
{
   x <- split.rna.matrix.by.rsid(mtx.tcx, rsid="rs28834970", suffix="_TCX")
   dims <- lapply(x, dim)
   checkEquals(names(dims), c("het", "hom", "mut", "wt"))
   checkEquals(as.numeric(lapply(dims, "[", 2)), c(122, 41, 163, 94))

} # test_split.rna.matrix.by.rsid
#----------------------------------------------------------------------------------------------------
run.wt.trena <- function(mtx.rna, tbl.tms.sub, extra.tfs=c())
{
    tfs <- unique(tbl.tms.sub$tf)
    solver <- EnsembleSolver(mtx.rna,
                             targetGene=targetGene,
                             candidateRegulators=c(tfs, extra.tfs),
                             solverNames=c("lasso", "Ridge", "Spearman", "Pearson", "bicor", "RandomForest", "xgboost"),
                             geneCutoff=0.9)
    tbl.trena <- run(solver)
    tbl.trena <- tbl.trena[order(abs(tbl.trena$spearman), decreasing=TRUE),]
    #tbl.trena <- tbl.trena[order(tbl.trena$rfScore, decreasing=TRUE),]
    tbl.trena$tfbs <- unlist(lapply(tbl.trena$gene, function(gene) nrow(subset(tbl.tms.sub, tf==gene))))
    rownames(tbl.trena) <- NULL
    numeric.colnames <- names(which(lapply(tbl.trena, class) == "numeric"))
    for(colname in numeric.colnames)
        tbl.trena[,colname] <- round(tbl.trena[,colname], digits=3)
    tbl.trena

} # run.wt.trena
#------------------------------------------------------------------------------------------------------------------------
show.models <- function(models, order.by="spearmanCoeff", n=10)
{
    f <- function(model){
       new.order <- order(model[[order.by]], decreasing=TRUE)
       head(model[new.order,], n=n)
       }

   lapply(models, f)


} # show.models
#----------------------------------------------------------------------------------------------------
run.trenas <- function()
{
   if(!exists("tbl.tms"))
      load("tbl.tms.RData")
    dim(tbl.tms) # 129446     21

    tbl.tms.sub <- subset(tbl.tms, (oc | chip) & gh > 600)
    dim(tbl.tms.sub)
    tbl.tms.sub <- subset(tbl.tms, (oc | chip) & gh > 600 & fimo_pvalue < 0.001)
    dim(tbl.tms.sub)
    tfs <- unique(tbl.tms.sub$tf)
    "RARA" %in% tfs
    "CUX2" %in% tfs
    length(tfs)  #483

    models <- list()
    models[["all.cer"]] <- run.wt.trena(mtx.cer, tbl.tms.sub)
    #models[["all.tcx"]] <- run.wt.trena(mtx.tcx, tbl.tms.sub)

    #mtx.split <- split.rna.matrix.by.rsid(mtx.tcx, suffix="_TCX")
    mtx.split <- split.rna.matrix.by.rsid(mtx.cer, rsid="rs28834970", suffix="_CER")
    names(mtx.split)
    lapply(mtx.split, ncol)   # het 122, hom 41, mut 163, wt 94

    models[["mut"]] <- run.wt.trena(mtx.split$mut, tbl.tms.sub)
    models[["wt"]] <- run.wt.trena(mtx.split$wt, tbl.tms.sub)
    models[["hom"]] <- run.wt.trena(mtx.split$hom, tbl.tms.sub)
    models[["het"]] <- run.wt.trena(mtx.split$het, tbl.tms.sub)
    models[["old.cov"]] <- run.wt.trena(mtx.old.tcx.cov, tbl.tms.sub)
    models[["rosmap"]] <- run.wt.trena(mtx.rosmap, tbl.tms.sub)

} # run.trenas
#----------------------------------------------------------------------------------------------------
explore.wt.mut.deltas <- function()
{
    # get models from "run.trenas" above

    tbl.wt <- models[["wt"]]
    tbl.mut <- models[["mut"]]   # 94 wt, 163 mut in mtx.cer

    dim(tbl.wt)      # 245
    dim(tbl.mut)     # 211

    tbl.bicor <- data.frame(tf=sort(unique(c(tbl.wt$gene, tbl.mut$gene))),  method="bicor", stringsAsFactors=FALSE)
    tbl.lasso <- data.frame(tf=sort(unique(c(tbl.wt$gene, tbl.mut$gene))),  method="lasso", stringsAsFactors=FALSE)
    tbl.ridge <- data.frame(tf=sort(unique(c(tbl.wt$gene, tbl.mut$gene))),  method="ridge", stringsAsFactors=FALSE)
    tbl.rf    <- data.frame(tf=sort(unique(c(tbl.wt$gene, tbl.mut$gene))),  method="rf", stringsAsFactors=FALSE)
    tbl.boost  <- data.frame(tf=sort(unique(c(tbl.wt$gene, tbl.mut$gene))), method="boost", stringsAsFactors=FALSE)

    tbl.bicor$wt <- 0
    tbl.bicor$mut <- 0

    tbl.lasso$wt <- 0
    tbl.lasso$mut <- 0

    tbl.ridge$wt <- 0
    tbl.ridge$mut <- 0

    tbl.rf$wt <- 0
    tbl.rf$mut <- 0

    tbl.boost$wt <- 0
    tbl.boost$mut <- 0

    coi <- "bicor"
    indices <- match(tbl.wt$gene, tbl.bicor$tf)
    tbl.bicor$wt[indices] <- tbl.wt[, coi]
    indices <- match(tbl.mut$gene, tbl.bicor$tf)
    tbl.bicor$mut[indices] <- tbl.mut[, coi]
    tbl.bicor[, "delta"] <- tbl.bicor$mut - tbl.bicor$wt
    new.order <- order(abs(tbl.bicor$delta), decreasing=TRUE)
    tbl.bicor <- tbl.bicor[new.order,]
    rownames(tbl.bicor) <- NULL
    head(tbl.bicor, n=10)

    coi <- "betaLasso"
    indices <- match(tbl.wt$gene, tbl.lasso$tf)
    tbl.lasso$wt[indices] <- tbl.wt[, coi]
    indices <- match(tbl.mut$gene, tbl.lasso$tf)
    tbl.lasso$mut[indices] <- tbl.mut[, coi]
    tbl.lasso[, "delta"] <- tbl.lasso$mut - tbl.lasso$wt
    new.order <- order(abs(tbl.lasso$delta), decreasing=TRUE)
    tbl.lasso <- tbl.lasso[new.order,]
    rownames(tbl.lasso) <- NULL
    head(tbl.lasso, n=10)

    coi <- "betaRidge"
    indices <- match(tbl.wt$gene, tbl.ridge$tf)
    tbl.ridge$wt[indices] <- tbl.wt[, coi]
    indices <- match(tbl.mut$gene, tbl.ridge$tf)
    tbl.ridge$mut[indices] <- tbl.mut[, coi]
    tbl.ridge[, "delta"] <- tbl.ridge$mut - tbl.ridge$wt
    new.order <- order(abs(tbl.ridge$delta), decreasing=TRUE)
    tbl.ridge <- tbl.ridge[new.order,]
    rownames(tbl.ridge) <- NULL
    head(tbl.ridge, n=10)

    coi <- "rfScore"
    indices <- match(tbl.wt$gene, tbl.rf$tf)
    tbl.rf$wt[indices] <- tbl.wt[, coi]
    indices <- match(tbl.mut$gene, tbl.rf$tf)
    tbl.rf$mut[indices] <- tbl.mut[, coi]
    tbl.rf[, "delta"] <- tbl.rf$mut - tbl.rf$wt
    new.order <- order(abs(tbl.rf$delta), decreasing=TRUE)
    tbl.rf <- tbl.rf[new.order,]
    rownames(tbl.rf) <- NULL
    head(tbl.rf, n=10)

    coi <- "xgboost"
    indices <- match(tbl.wt$gene, tbl.boost$tf)
    tbl.boost$wt[indices] <- tbl.wt[, coi]
    indices <- match(tbl.mut$gene, tbl.boost$tf)
    tbl.boost$mut[indices] <- tbl.mut[, coi]
    tbl.boost[, "delta"] <- tbl.boost$mut - tbl.boost$wt
    new.order <- order(abs(tbl.boost$delta), decreasing=TRUE)
    tbl.boost <- tbl.boost[new.order,]
    rownames(tbl.boost) <- NULL
    head(tbl.boost, n=10)

    max <- 5
    max <- 10;

    tbl.joint <- do.call(rbind, list(tbl.bicor[1:max,], tbl.lasso[1:max,], tbl.boost[1:max,], tbl.ridge[1:max,], tbl.rf[1:max,]))
    tbl.counts <- as.data.frame(sort(table(tbl.joint$tf), decreasing=TRUE))
    candidates <- as.character(subset(tbl.counts, Freq > 2)$Var1)
    tbl.joint.strong <- subset(tbl.joint, tf %in% candidates)
    new.order <- order(tbl.joint.strong$tf)
    tbl.joint.strong <- tbl.joint.strong[new.order,]
    rownames(tbl.joint.strong) <- NULL

    common <- Reduce(intersect, list(
                                     tbl.bicor$tf[1:max],
                                     #tbl.rf$tf[1:max],
                                     #tbl.boost$tf[1:max],
                                     #tbl.rf$tf[1:max],
                                     tbl.lasso$tf[1:max],
                                     #tbl.ridge$tf[1:max],
                                     tbl.mut$gene
                                     ))
    print(common)
    tbls <- list(tbl.bicor, tbl.lasso, tbl.rf, tbl.ridge, tbl.boost)
    for(tbl in tbls)
        print(subset(tbl, tf %in% common))
    tbl.counts <- as.data.frame(sort(table(tbl$tf), decreasing=FALSE))
    subset(tbl.counts, Freq > 1)

    fivenum(abs(tbl.deltas[, coi]))
    fivenum(tbl.deltas[, coi])
    subset(tbl.deltas, abs(bicor) > 0.3)

    coi <- "spearmanCoeff"
    indices <- match(tbl.wt$gene, tbl.deltas$tf)
    tbl.deltas$wt[indices] <- tbl.wt[, coi]


    indices <- match(tbl.mut$gene, tbl.deltas$tf)
    tbl.deltas$mut[indices] <- tbl.mut[, coi]
    tbl.deltas[, coi] <- tbl.deltas$mut - tbl.deltas$wt
    fivenum(abs(tbl.deltas[, coi]))

    coi <- "betaLasso"
    indices <- match(tbl.wt$gene, tbl.deltas$tf)
    tbl.deltas$wt[indices] <- tbl.wt[, coi]

    indices <- match(tbl.mut$gene, tbl.deltas$tf)
    tbl.deltas$mut[indices] <- tbl.mut[, coi]
    tbl.deltas[, coi] <- tbl.deltas$mut - tbl.deltas$wt
    fivenum(abs(tbl.deltas[, coi]))  #  0.000 0.000 0.000 0.000 0.578
    subset(tbl.deltas, abs(betaLasso) > 0.5)

    coi <- "betaRidge"
    indices <- match(tbl.wt$gene, tbl.deltas$tf)
    tbl.deltas$wt[indices] <- tbl.wt[, coi]

    indices <- match(tbl.mut$gene, tbl.deltas$tf)
    tbl.deltas$mut[indices] <- tbl.mut[, coi]
    tbl.deltas[, coi] <- tbl.deltas$mut - tbl.deltas$wt
    fivenum(abs(tbl.deltas[, coi]))  #  0.000 0.005 0.012 0.022 0.189
    subset(tbl.deltas, abs(betaRidge) > 0.10)

    coi <- "rfScore"
    indices <- match(tbl.wt$gene, tbl.deltas$tf)
    tbl.deltas$wt[indices] <- tbl.wt[, coi]

    indices <- match(tbl.mut$gene, tbl.deltas$tf)
    tbl.deltas$mut[indices] <- tbl.mut[, coi]
    tbl.deltas[, coi] <- tbl.deltas$mut - tbl.deltas$wt
    fivenum(abs(tbl.deltas[, coi]))  #  0.000 0.005 0.012 0.022 0.189
    subset(tbl.deltas, abs(rfScore) > 2.4)


} # explore.wt.mut.deltas
#----------------------------------------------------------------------------------------------------
plot.pbx2 <- function()
{
    plot(mtx.split$wt["PBX2",], mtx.split$wt["PTK2B",], pch=18, col="darkgreen", xlab="PBX2", ylab="PTK2B",
         main="spearman & bicor not confused by the outlier", xlim=c(0,6), ylim=c(0,6))
    points(mtx.split$mut["PBX2",], mtx.split$mut["PTK2B",], pch=18, col="red")
    legend(5,6,  c("wt", "mut"), c("black", "red"))

    plot(mtx.split$wt["EGR4",], mtx.split$wt["PTK2B",], pch=18, col="darkgreen", xlab="EGR4", ylab="PTK2B",
         main="spearman & bicor not confused by the outlier", xlim=c(-2,6), ylim=c(0,6))
    points(mtx.split$mut["EGR4",], mtx.split$mut["PTK2B",], pch=18, col="red")
    legend(5,6,  c("wt", "mut"), c("black", "red"))

    t.test(mtx.split$wt["PBX2",], mtx.split$wt["PTK2B",])
    t.test(mtx.split$wt["EGR4",], mtx.split$wt["PTK2B",])






} # plot.pbx2
#----------------------------------------------------------------------------------------------------
