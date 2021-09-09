library(TrenaMultiScore)
library(TrenaProjectAD)
library(trena)
library(EndophenotypeExplorer)
library(RUnit)

targetGene <- "PTK2B"
etx <- EndophenotypeExplorer$new(targetGene, "hg38")

#------------------------------------------------------------------------------------------------------------------------
run.wt.trena <- function(mtx.rna, tbl.sub, extra.tfs=c())
{
    tfs <- unique(tbl.sub$tf)
    solver <- EnsembleSolver(mtx.rna,
                             targetGene=targetGene,
                             candidateRegulators=c(tfs, extra.tfs),
                             solverNames=c("lasso", "Ridge", "Spearman", "Pearson", "RandomForest", "xgboost"),
                             geneCutoff=0.9)
    tbl.trena <- run(solver)
    tbl.trena <- tbl.trena[order(tbl.trena$rfScore, decreasing=TRUE),]
    tbl.trena$tfbs <- unlist(lapply(tbl.trena$gene, function(gene) nrow(subset(tbl.sub, tf==gene))))
    rownames(tbl.trena) <- NULL
    numeric.colnames <- names(which(lapply(tbl.trena, class) == "numeric"))
    for(colname in numeric.colnames)
        tbl.trena[,colname] <- round(tbl.trena[,colname], digits=3)
    tbl.trena

} # run.wt.trena
#------------------------------------------------------------------------------------------------------------------------
# this sort of works - but the variant filtered mtx.rna are too small for sensible trena
#
run.variant.trena <- function(mtx.rna, tbl.breaks, rsid)
{
   # identify gain of tfbs TFs
  tbl.rsid <- subset(tbl.breaks, SNP_id==rsid & effect == "strong" & pctAlt > 0.80)
  mtx.geno <- etx$getGenoMatrixByRSID(rsid)
  new.names <- etx$locsToRSID(rownames(mtx.geno), "hg19")
  dim(mtx.geno)
  table(mtx.geno)   #  0/0  0/1  1/1
                    # 1621  265    8
  hets <- colnames(mtx.geno)[which(mtx.geno[1,] == "0/1")]
  homs <- colnames(mtx.geno)[which(mtx.geno[1,] == "1/1")]
  wt   <- colnames(mtx.geno)[which(mtx.geno[1,] == "0/0")]
  checkEquals(sum(length(wt), length(homs), length(hets)), ncol(mtx.geno))
  dir <- "~/github/EndophenotypeExplorer/inst/extdata/idMapping"
  file <- "tbl.sampleToPatientMap-rosmap-mayo-sinai.RData"

  tbl.ids <- get(load(file.path(dir, file)))

  homs.patients <- subset(tbl.ids, sample %in% homs)$patient
  rna.seq.sample.ids <- subset(tbl.ids, patient %in% homs.patients & assay=="rnaseq" & study=="mayo")$sample
  colnames(mtx.rna) <- sub("^X", "", colnames(mtx.rna))
  rna.seq.sample.ids %in% colnames(mtx.rna)   # just 2/4
  mtx.rna.homs <- mtx.rna[, intersect(rna.seq.sample.ids, colnames(mtx.rna))]
  dim(mtx.rna.homs)

  hets.patients <- subset(tbl.ids, sample %in% hets)$patient
  rna.seq.sample.ids <- subset(tbl.ids, patient %in% hets.patients & assay=="rnaseq" & study=="mayo")$sample # 68
  colnames(mtx.rna) <- sub("^X", "", colnames(mtx.rna))
  rna.seq.sample.ids %in% colnames(mtx.rna)   # just 34/34
  mtx.rna.hets <- mtx.rna[, intersect(rna.seq.sample.ids, colnames(mtx.rna))]
  dim(mtx.rna.hets) # [1] 16969    34

  tfs <- unique(c(tbl.tms.sub$tf, tbl.rsid$geneSymbol))
  solver <- EnsembleSolver(mtx.rna.homs,
                           targetGene=targetGene,
                           candidateRegulators=tfs,
                           solverNames=c("lasso", "Ridge", "Spearman", "Pearson", "RandomForest", "xgboost"),
                           geneCutoff=0.9)
  tbl.trena <- run(solver)
  tbl.trena <- tbl.trena[order(tbl.trena$rfScore, decreasing=TRUE),]
  tbl.trena$tfbs <- unlist(lapply(tbl.trena$gene, function(gene) nrow(subset(tbl.sub, tf==gene))))
  rownames(tbl.trena) <- NULL


} # run.variant.trena
#------------------------------------------------------------------------------------------------------------------------
run.tms <- function()
{
    tpad <- TrenaProjectAD()
    fimo.file <- "tbl.fimo.PTK2B-atac.RData"
    file.exists(fimo.file)
    tbl.fimo <- get(load(fimo.file))
    dim(tbl.fimo) # 45208 9

    f <- "~/github/TrenaProjectAD/explore/ptk2b/tbl.oc.ptk2b.atac.eqtl.merged.RData"
    file.exists(f)
    tbl.oc <- get(load(f))
    dim(tbl.oc)

# loc.chrom <- tbl.fimo$chrom[1]
# shoulder <- 1000
# loc.start <- min(tbl.fimo$start) - shoulder
# loc.end   <- max(tbl.fimo$end)   + shoulder

    tms <- TrenaMultiScore(tpad, targetGene, tbl.fimo, tbl.oc, quiet=FALSE)
    getGeneHancerRegion(tms)

    #  mtx.rna <- get(load("~/github/TrenaProjectAD/inst/extdata/expression/mayo.tcx.16969x262.covariateCorrection.log+scale.RData"
    f <- "mtx.mayo.tcx.eqtl-optimized-geneSymbols-sampleIDs-17009x262.RData"
    dir <- "~/github/TrenaProjectAD/prep/rna-seq-counts-from-synapse/eqtl"
    mtx.rna <- get(load(file.path(dir, f)))
    dim(mtx.rna) # 17009   262
    mtx.rna[is.na(mtx.rna)] <- 0

    addGeneExpressionCorrelations(tms, mtx.rna)

                                        # with spearman: fivenum(tbl.tms$cor) [1] -0.21 -0.02  0.00  0.03  0.31
    tbl.tms <- getMultiScoreTable(tms)
    fivenum(tbl.tms$cor)

    scoreMotifHitsForGeneHancer(tms)
    scoreMotifHitsForConservation(tms)
    scoreMotifHitsForOpenChromatin(tms)
                                        # scoreMotifHitsForOpenChromatin(tms)
    addDistanceToTSS(tms)
    addGenicAnnotations(tms)
    addChIP(tms)

    tbl.tms <- getMultiScoreTable(tms)
    dim(tbl.tms)   # 45208  17
    save(tbl.tms, file="tbl.tms.merged-eqtl-atac.RData")


} # run.tms
#------------------------------------------------------------------------------------------------------------------------
#
#     tbl.rich.01 <- subset(tbl.tms, abs(cor) > 0.9 & gh > 100 & chip & oc)
#     dim(tbl.rich.01)
#     as.data.frame(sort(table(tbl.rich.01$tf)))
#
#     tbl.rich.02 <- subset(tbl.tms, abs(cor) > 0.9 & gh > 600)
#     dim(tbl.rich.02)
#     as.data.frame(sort(table(tbl.rich.02$tf)))
#
#     tbl.rich.03 <- subset(tbl.tms, abs(cor) > 0.7 & gh > 600 &
#                                    ((fimo_pvalue < 1e-6 & oc) | (chip & oc)))
#     dim(tbl.rich.03)
#     as.data.frame(sort(table(tbl.rich.03$tf)))
#
#     tfs <- unique(tbl.rich.03$tf)
# length(tfs)
#
# tbl.trena.03 <- run.wt.trena(mtx.rna, tbl.rich.03, extra.tfs=c())
# head(tbl.trena.03, n=10)
#      gene betaLasso betaRidge spearmanCoeff pearsonCoeff  rfScore xgboost tfbs
# 1   MEF2C     0.228     0.104         0.850        0.871 2215.261   0.572    2
# 2     SP1     0.502     0.075         0.872        0.935 1461.969   0.313  160
# 3    NFYA     0.000     0.022         0.841        0.874  696.742   0.000   11
# 4   CEBPB     0.000     0.023         0.867        0.843  672.027   0.000    5
# 5    RFX5     0.000     0.018         0.847        0.924  616.054   0.001   10
# 6  NFE2L2     0.323     0.233         0.848        0.867  613.361   0.003    1
# 7    TAF1     0.000     0.030         0.848        0.897  485.185   0.001   34
# 8   TFDP1     0.000     0.074         0.727        0.822  300.302   0.000    7
# 9   STAT3     0.000     0.017         0.850        0.931  291.661   0.000   16
# 10   MXI1     0.024     0.062         0.859        0.911  270.636   0.004   20

# mtx.rna.rosmap <- mtx.rosmap
# mtx.rna.rosmap[is.na(mtx.rna.rosmap)] <- 0
# tbl.trena.04 <- run.wt.trena(mtx.rna.rosmap, tbl.rich.03, extra.tfs=c())
# head(tbl.trena.04, n=10)
#      gene betaLasso betaRidge spearmanCoeff pearsonCoeff rfScore xgboost tfbs
# 1    MXI1     0.000     0.053        -0.057        0.172  28.312   0.094   20
# 2   MEF2C     0.198     0.109         0.039        0.349  21.982   0.053    2
# 3  NFE2L2     0.364     0.248        -0.030        0.264  21.491   0.042    1
# 4   CREB1     1.090     0.602         0.006        0.345  20.893   0.159    3
# 5    LYL1     0.000     0.002         0.017       -0.004  17.450   0.109    2
# 6     ZFX     0.000     0.054        -0.024        0.007  17.424   0.100    7
# 7   CEBPB     0.134     0.148         0.019        0.299  15.441   0.001    5
# 8   NR2F1     0.000    -0.035        -0.132       -0.063  14.695   0.013    9
# 9    ELF1     0.000     0.082         0.017        0.127  14.091   0.088   58
# 10   ETS1     0.000     0.109        -0.006        0.220  11.725   0.001   24

# mtx.rna.sinai <- mtx.sinai
# mtx.rna.sinai[is.na(mtx.rna.sinai)] <- 0
# tbl.trena.05 <- run.wt.trena(mtx.rna.sinai, tbl.rich.03, extra.tfs=c())
# head(tbl.trena.05, n=10)
# #      gene betaLasso betaRidge spearmanCoeff pearsonCoeff rfScore xgboost tfbs
# # 1   MEF2C     0.245     0.117         0.016        0.584 104.685   0.198    2
# # 2     SP1     0.819     0.355        -0.016        0.536  99.448   0.222  160
# # 3  NFE2L2     0.327     0.216        -0.027        0.352  41.326   0.045    1
# # 4    ELK1     0.296     0.207         0.226        0.381  36.868   0.104   14
# # 5    ZEB1     0.000    -0.044        -0.247       -0.047  34.171   0.001   69
# # 6   FOXP1     0.087     0.198         0.082        0.365  29.059   0.000    3
# # 7  HMBOX1     0.000     0.018        -0.229        0.120  28.007   0.079    1
# # 8    NFYA     0.054     0.134        -0.085        0.377  24.259   0.001   11
# # 9   FOXK2     0.000    -0.016         0.114       -0.006  23.183   0.007    1
# # 10  CEBPB     0.000     0.121         0.222        0.228  21.194   0.000    5
# mtx.rna.mayo <- mtx.mayo
# mtx.rna.mayo[is.na(mtx.rna.mayo)] <- 0
# tbl.trena.06 <- run.wt.trena(mtx.rna.mayo, tbl.rich.03, extra.tfs=c())
# head(tbl.trena.06, n=10)
# #     gene betaLasso betaRidge spearmanCoeff pearsonCoeff rfScore xgboost tfbs
# # 1   ELK1     0.205     0.132         0.188        0.343  12.180   0.078   14
# # 2   MAFK    -0.047    -0.022        -0.564       -0.421  12.166   0.189    2
# # 3  NR2F1     0.000    -0.016        -0.272       -0.037   8.764   0.047    9
# # 4  TEAD1    -0.114    -0.029        -0.577       -0.378   8.518   0.031    1
# # 5   ETV5     0.000    -0.012        -0.342       -0.262   7.994   0.105    1
# # 6    ERF     0.000    -0.009        -0.446       -0.333   6.477   0.056    1
# # 7   KLF9     0.000     0.014         0.540        0.349   6.016   0.077   24
# # 8   ELK4    -0.194    -0.227        -0.158       -0.213   5.693   0.003   16
# # 9   KLF4     0.000    -0.027        -0.288       -0.241   3.911   0.047   57
# # 10  REST     0.000    -0.008        -0.536       -0.375   3.634   0.033   26
#
#
# # change the tms thresholds to include CUX2
#
# "CUX2" %in%    subset(tbl.tms, abs(cor) > 0.7 & gh > 600 &((fimo_pvalue < 1e-3 & oc) | (chip & oc)))$tf
# tbl.rich.04 <- subset(tbl.tms, abs(cor) > 0.7 & gh > 600 &((fimo_pvalue < 1e-3 & oc) | (chip & oc)))
#
# dim(tbl.rich.04)
# as.data.frame(sort(table(tbl.rich.04$tf)))
#
# tfs <- unique(tbl.rich.04$tf)
# length(tfs)
#
# tbl.trena.04 <- run.wt.trena(mtx.rna, tbl.rich.04, extra.tfs=c())
# head(tbl.trena.04, n=10)
# #     gene betaLasso betaRidge spearmanCoeff pearsonCoeff  rfScore xgboost tfbs
# # 1  MEF2C     0.159     0.083         0.850        0.871 2077.938   0.444   13
# # 2    SP1     0.349     0.064         0.872        0.935  892.147   0.279  316
# # 3  MEIS3     0.045     0.055         0.805        0.860  880.139   0.060    7
# # 4  SMAD2     0.146     0.067         0.860        0.911  416.603   0.023   19
# # 5   NFYA     0.000     0.007         0.841        0.874  397.745   0.000   33
# # 6  ZNF24     0.000     0.047         0.843        0.889  355.195   0.001    2
# # 7  STAT4     0.000     0.007         0.869        0.861  349.824   0.090    9
# # 8  CEBPB     0.000     0.036         0.867        0.843  325.900   0.000    7
# # 9   RORA     0.000     0.025         0.856        0.927  323.007   0.000    2
# # 10   MGA     0.000     0.012         0.776        0.858  298.444   0.000   10
#
#
# tbl.trena.05 <- run.wt.trena(mtx.rosmap, tbl.rich.04, extra.tfs=c())
# head(tbl.trena.05, n=10)
#
# #       gene betaLasso betaRidge spearmanCoeff pearsonCoeff rfScore xgboost tfbs
# # 1   MEF2C     0.245     0.117         0.016        0.584 104.685   0.198    2
# # 2     SP1     0.819     0.355        -0.016        0.536  99.448   0.222  160
# # 3  NFE2L2     0.327     0.216        -0.027        0.352  41.326   0.045    1
# # 4    ELK1     0.296     0.207         0.226        0.381  36.868   0.104   14
# # 5    ZEB1     0.000    -0.044        -0.247       -0.047  34.171   0.001   69
# # 6   FOXP1     0.087     0.198         0.082        0.365  29.059   0.000    3
# # 7  HMBOX1     0.000     0.018        -0.229        0.120  28.007   0.079    1
# # 8    NFYA     0.054     0.134        -0.085        0.377  24.259   0.001   11
# # 9   FOXK2     0.000    -0.016         0.114       -0.006  23.183   0.007    1
# # 10  CEBPB     0.000     0.121         0.222        0.228  21.194   0.000    5
# >
#
# mtx.sinai[is.na(mtx.sinai)] <- 0
#
# tbl.trena.06 <- run.wt.trena(mtx.sinai, tbl.rich.04, extra.tfs=c())
# head(tbl.trena.06, n=10)
#
# #     gene betaLasso betaRidge spearmanCoeff pearsonCoeff rfScore xgboost tfbs
# # 1   ELK1     0.205     0.132         0.188        0.343  12.180   0.078   14
# # 2   MAFK    -0.047    -0.022        -0.564       -0.421  12.166   0.189    2
# # 3  NR2F1     0.000    -0.016        -0.272       -0.037   8.764   0.047    9
# # 4  TEAD1    -0.114    -0.029        -0.577       -0.378   8.518   0.031    1
# # 5   ETV5     0.000    -0.012        -0.342       -0.262   7.994   0.105    1
# # 6    ERF     0.000    -0.009        -0.446       -0.333   6.477   0.056    1
# # 7   KLF9     0.000     0.014         0.540        0.349   6.016   0.077   24
# # 8   ELK4    -0.194    -0.227        -0.158       -0.213   5.693   0.003   16
# # 9   KLF4     0.000    -0.027        -0.288       -0.241   3.911   0.047   57
# # 10  REST     0.000    -0.008        -0.536       -0.375   3.634   0.033   26
#
# mtx.mayo[is.na(mtx.mayo)] <- 0
# tbl.trena.07 <- run.wt.trena(mtx.mayo, tbl.rich.04, extra.tfs=c())
#
# head(tbl.trena.07, n=10)
# #      gene betaLasso betaRidge spearmanCoeff pearsonCoeff rfScore xgboost tfbs
# # 1    ELK1     0.150     0.100         0.188        0.343  13.160   0.065   26
# # 2   MEF2D     0.256     0.119         0.493        0.533   8.456   0.095    8
# # 3  MLXIPL     0.096     0.044         0.598        0.469   7.759   0.040   11
# # 4    MAFK     0.000    -0.014        -0.564       -0.421   5.780   0.140    9
# # 5   TEAD1     0.000    -0.010        -0.577       -0.378   5.643   0.015    9
# # 6    ETV5     0.000    -0.003        -0.342       -0.262   5.205   0.037   13
# # 7   FOXP2     0.000    -0.011        -0.429       -0.319   5.021   0.001    2
# # 8   FOXO3     0.000     0.048         0.600        0.436   4.012   0.012    5
# # 9     ERF     0.000    -0.011        -0.446       -0.333   3.827   0.065   11
# # 10  KLF13     0.000     0.058         0.522        0.363   3.779   0.035   11
#
#
#------------------------------------------------------------------------------------------------------------------------
test_run.variant.trena <- function()
{
   rsid <- "rs12754503"  # very strong GATA2 motif creator, from 65% to 99%
   tbl.trena <- run.variant.trena(mtx.rna, tbl.breaks, rsid)

} # test_run.variant.trena
#------------------------------------------------------------------------------------------------------------------------
# tbl.tms.sub <-
#     subset(tbl.tms, abs(cor) > 0.8 & gh > 10 & fimo_pvalue < 1e-6)
# dim(tbl.tms.sub)
# length(unique(tbl.tms.sub$tf))
# tbl.trena <- run.wt.trena(mtx.rna, tbl.tms.sub)
# print(tbl.trena)
#
# for(TF in tbl.trena$gene){   # have any existing motifs been broken?
#     printf("--- %s", TF)
#     tbl.tms.tf <- subset(tbl.tms.sub, tf==TF)
#     tbl.breaks.tf <- subset(tbl.breaks2, geneSymbol==TF & effect=="strong")
#     tbl.ov <- as.data.frame(findOverlaps(GRanges(tbl.breaks.tf), GRanges(tbl.tms.tf)))
#     if(nrow(tbl.ov) > 0){
#        print(tbl.ov)
#        tms.rows <- unique(tbl.ov$subjectHits)
#        break.rows <- unique(tbl.ov$queryHits)
#        print(tbl.tms.tf[tms.rows,])
#        print(tbl.breaks.tf[break.rows,])
#       }
#     }
#
# fimo.tfs <- unique(tbl.tms$tf)
# tfs.with.gain <- unique(intersect(subset(tbl.breaks2, alleleDiff > 1.0)$geneSymbol, fimo.tfs)) # 51
# length(tfs.with.gain)
# new.tfs.with.gain <- setdiff(tfs.with.gain, tbl.trena$gene)
# length(new.tfs.with.gain)
# tbl.trena <- run.wt.trena(mtx.rna, tbl.tms.sub, tfs.with.gain)
#
