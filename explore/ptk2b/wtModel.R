library(TrenaMultiScore)
library(TrenaProjectAD)
library(trena)
library(EndophenotypeExplorer)
library(RUnit)
#----------------------------------------------------------------------------------------------------
mtx.rna <- get(load("~/github/TrenaProjectAD/inst/extdata/expression/temporalCortex.15167x264.RData"))
fivenum(mtx.rna)
#----------------------------------------------------------------------------------------------------
targetGene <- "PTK2B"
tpad <- TrenaProjectAD()

fimo.file <- "tbl.fimo.PTK2B-atac.RData"
file.exists(fimo.file)
tbl.fimo <- get(load(fimo.file))
dim(tbl.fimo) # 45208 9

f <- "~/github/TrenaProjectAD/explore/ptk2b/tbl.oc.ptk2b.atac.eqtl.merged.RData"
file.exists(f)
tbl.oc <- get(load(f))
dim(tbl.oc)

#tbl.oc <- get(load("~/github/TrenaProjectAD/explore/mayo-epigenetics/atac/mayoAllPeaks.1052789x4.RData"))
#dim(tbl.oc)
tms <- TrenaMultiScore(tpad, targetGene, tbl.fimo, tbl.oc, quiet=FALSE)

dim(mtx.rna)  # 15167 264
addGeneExpressionCorrelations(tms, mtx.rna)

  # with spearman: fivenum(tbl.tms$cor) [1] -0.21 -0.02  0.00  0.03  0.31

scoreMotifHitsForGeneHancer(tms)
scoreMotifHitsForConservation(tms)
scoreMotifHitsForOpenChromatin(tms)
# scoreMotifHitsForOpenChromatin(tms)
addDistanceToTSS(tms)
addGenicAnnotations(tms)
addChIP(tms)

tbl.tms <- getMultiScoreTable(tms)
fivenum(tbl.tms$cor)

tbl.rich.03 <- subset(tbl.tms, abs(cor) > 0.4 & gh > 10 & ((fimo_pvalue < 1e-4 & oc) | (chip & oc)) & phast7 > 0.2)
tfs <- unique(tbl.rich.03$tf)
length(tfs)
tbl.tfs <- as.data.frame(sort(table(tbl.rich.03$tf), decreasing=TRUE))
extra.tfs <- c()
solver <- EnsembleSolver(mtx.rna,
                         targetGene=targetGene,
                         candidateRegulators=c(tfs, extra.tfs),
                         solverNames=c("lasso", "Ridge", "Spearman", "Pearson", "RandomForest", "xgboost", "bicor"),
                         geneCutoff=0.9)
tbl.trena <- run(solver)
tbl.trena <- tbl.trena[order(tbl.trena$rfScore, decreasing=TRUE),]
tbl.trena <- tbl.trena[order(abs(tbl.trena$bicor), decreasing=TRUE),]
tbl.trena$tfbs <- unlist(lapply(tbl.trena$gene, function(gene) nrow(subset(tbl.rich.03, tf==gene))))
rownames(tbl.trena) <- NULL
numeric.colnames <- names(which(lapply(tbl.trena, class) == "numeric"))
for(colname in numeric.colnames)
    tbl.trena[,colname] <- round(tbl.trena[,colname], digits=3)
tbl.trena

rf.order <- order(tbl.trena$rfScore, decreasing=TRUE)
head(tbl.trena[rf.order,])
length(tbl.trena[rf.order,]$gene) # 50

lasso.order <- order(abs(tbl.trena$betaLasso), decreasing=TRUE)
head(tbl.trena[lasso.order,], n=20)

lasso.keepers <- subset(tbl.trena, abs(betaLasso) > 0)$gene
rf.keepers    <- subset(tbl.trena, rfScore > 4)$gene
bicor.keepers <- subset(tbl.trena, abs(bicor) > 0.5)$gene
tf.keepers <- sort(unique(c(lasso.keepers, rf.keepers, bicor.keepers)))
length(tf.keepers)

tfs.2 <- c("GABPA", "PKNOX2")
tfs.50 <- c("GABPA", "PKNOX2", "MEF2C", "STAT4", "ELK4", "NFATC3", "SP3",
            "EGR3", "BACH1", "REST", "TBR1", "BHLHE41", "NEUROD1", "MGA",
            "E2F3", "NFYC", "SP1", "EGR4", "HEY2", "RARA", "FOXO1", "ELF1",
            "PPARG", "SRF", "FOSB", "TP73", "MEIS3", "NEUROD2", "IRF2",
            "ZFX", "NFYA", "MTF1", "ATF3", "NFE2L2", "TAF1", "TCF3", "RBPJ",
            "ZBTB7B", "MEF2D", "TFEB", "CTCF", "PRRX1", "RELA", "MXI1",
            "SP4", "TGIF2", "MAFK", "SMAD2", "NFIA", "NFKB1")

motif.ids <- unique(subset(tbl.fimo, tf %in% tfs.2)$motif_id)
mdb.2 <- MotifDb[motif.ids]
save(mdb.2, file="motifs.2tfs.3motifs.pkt2b.RData")


amotif.ids <- unique(subset(tbl.fimo, tf %in% tf.keepers)$motif_id)
length(motif.ids)  # 60
mdb.38 <- MotifDb[motif.ids]
save(mdb.38, file="motifs.38tfs.60motifs.pkt2b.RData")


mdb.50 <- MotifDb[motif.ids] # MotifDb object of length 80
save(mdb.50, file="motifs.50.80.pkt2b.RData")
