library(RUnit)
source("~/github/TrenaProjectAD/explore/rs4575098/tmsCore.R")
trenaProject <- TrenaProjectAD()

print(load("~/github/TrenaProjectAD/explore/mayo-epigenetics/atac/earlyResultMayoWholeGenome.RData"))
colnames(tbl.dba)[1] <- "chrom"
tbl.dba$chrom <- as.character(tbl.dba$chrom)
tbl.atac.merged <- tbl.dba


runTMS <- function(targetGene, chrom, start, end, mtx.rna, fimo.threshold=1e-4, verbose=TRUE)
{

    x <- TMS$new(targetGene, trenaProject)
    x$addGeneHancer()
    tbl.gh <- x$getGeneHancer()
    #checkEquals(dim(tbl.gh), c(25, 16))
    tbl.ghRegion <- x$getGeneHancerRegion()
    #checkEquals(tbl.ghRegion$start, min(tbl.gh$start) - 1000)
    #checkEquals(tbl.ghRegion$end, max(tbl.gh$end) + 1000)
    # classic promoter for DOT1L, plus a few kb: chr19:2,158,665-2,177,305
    # intersect.with.genehancer T/F get either 4 or 7 areas of open chromatin
    # all gh+atac-merged found in this 2 MB region: chr19:917,289-3,303,350
    #chrom <- "chr19"
    #start <- 2158665
    #end   <- 2177305
    #start <- 917289
    #end   <- 3303350
    printf("retrieving open chromatin for %s, %5.1f kb", targetGene, (end-start)/1000)
    x$addOpenChromatin(chrom, start, end,
                       intersect.with.genehancer=FALSE,
                       promoter.only=FALSE)
    tbl.oc <- x$getOpenChromatin()

    # x$viewOpenChromatin()
    #checkEquals(dim(tbl.oc), c(88, 7))
    printf("querying MotifDb")
    motifs <- query(MotifDb, c("sapiens"), c("jaspar2018", "hocomocov11-core-a"))
    printf("running FIMO, motifs: %d, threshold: %f", length(motifs), fimo.threshold)
    x$addAndScoreFimoTFBS(fimo.threshold=fimo.threshold, motifs=motifs)
    tbl.tms <- x$getTfTable()
    dim(tbl.tms)  # 165 16
    length(unique(tbl.tms$tf))  # [1] 44
    printf("adding RNA binding proteins")
    #x$addRBP()
    #tbl.rbp <- x$getRbpTable()
    #dim(tbl.rbp)  # 7652  12
    printf("adding mRNA tf correlations")
    x$add.tf.mrna.correlations(mtx.rna, featureName="cor.all")
    print(1)
    tbl.tms <- x$getTfTable()
    print(2)
    print(dim(tbl.tms))
    print(3)
    browser()
    print(fivenum(tbl.tms$cor.all))
    print(4)
    #printf("adding mRNA rbp correlations")
    #x$add.rbp.mrna.correlations(mtx.rna, featureName="cor.all")
    #tbl.rbp <- x$getRbpTable()
    #dim(tbl.rbp)
    #fivenum(tbl.rbp$cor.all)
    #dim(tbl.rbp)  # 7652   13

    tfs <- unique(subset(tbl.tms, abs(cor.all) > 0.5)$tf)
    printf("candidate tfs: %d", length(tfs))
    #rbps <- unique(subset(tbl.rbp, celltype=="K562" & abs(cor.all) > 0.5)$gene)
    #printf("candidate rbps in K562 cells: %d", length(rbps))
    #tbl.trena.tf <- x$build.trena.model(tfs, list(), mtx.rna)
    #dim(tbl.trena.tf)  # 8 7
    rbps <- c()
    printf("running trena")
    tbl.trena.both <- x$build.trena.model(tfs, rbps, mtx.rna)
    linear.models <- list()
    linear.models[["spearman"]] <- x$build.linear.model(mtx.rna, sort.by.column="spearmanCoeff", candidate.regulator.max=15)
    linear.models[["pearson"]] <- x$build.linear.model(mtx.rna, sort.by.column="spearmanCoeff", candidate.regulator.max=15)
    linear.models[["betaLasso"]] <- x$build.linear.model(mtx.rna, sort.by.column="betaLasso", candidate.regulator.max=15)
    linear.models[["betaRidge"]] <- x$build.linear.model(mtx.rna, sort.by.column="betaRidge", candidate.regulator.max=15)
    linear.models[["rfScore"]] <- x$build.linear.model(mtx.rna, sort.by.column="rfScore", candidate.regulator.max=15)
    linear.models[["xgboose"]] <- x$build.linear.model(mtx.rna, sort.by.column="xgboost", candidate.regulator.max=15)
    printf("--- modeling %s", targetGene)
    #print(subset(tbl.model, p.value <= 0.2))
    #printf("adjusted r-squared: %5.2f", x$get.lm.Rsquared())
    return(x)

} # runTMS
