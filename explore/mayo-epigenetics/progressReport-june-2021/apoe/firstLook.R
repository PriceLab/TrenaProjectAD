library(ghdb)
ghdb <- GeneHancerDB()
targetGene <- "MAPT"
targetGene <- "APOE"
targetGene <- "TREM2"
igv <- start.igv(targetGene, "hg38")
showGenomicRegion(igv, targetGene)

tbl.gh <- retrieveEnhancersFromDatabase(ghdb, targetGene, tissues="all")

tbl.track <- tbl.gh[, c("chrom", "start", "end", "combinedscore")]
colnames(tbl.track)[4] <- "score"
tbl.track$score <- asinh(tbl.track$score)
track <- DataFrameQuantitativeTrack("gh", tbl.track, autoscale=TRUE, color="red")
displayTrack(igv, track)

shoulder <- 10000
roi <- sprintf("%s:%d-%d",
               tbl.gh$chrom[1],
               min(tbl.gh$start) - shoulder,
               max(tbl.gh$end) + shoulder)
showGenomicRegion(igv, roi)

# zoomOut(igv)
# zoomOut(igv)
# zoomOut(igv)

tbl.atac <- get(load("~/github/TrenaProjectAD/explore/mayo-epigenetics/atac/dbaConsensusRegionsScored.74273x30.RData"))
roi <- getGenomicRegion(igv)
tbl.atac.sub <- subset(tbl.atac, chrom==roi$chrom & start >= roi$start & end <= roi$end)
dim(tbl.atac.sub)

for(col in c("AD", "PSP", "CTL")){
   tbl.track <- tbl.atac.sub[, c("chrom", "start", "end", col)]
   fivenum(tbl.track$AD)
   track <- DataFrameQuantitativeTrack(col, tbl.track, color="random", autoscale=FALSE, min=0, max=0.025)
   displayTrack(igv, track)
   }



source("~/github/fimoService/batchMode/fimoBatchTools.R")
meme.file <- "jaspar2018-hocomocoCore.meme"
motifs <- query(MotifDb, c("sapiens"), c("jaspar2018", "HOCOMOCOv11-core-A"))
length(motifs) # 718

export(motifs, con=meme.file, format="meme")

tbl.fimo <- fimoBatch(tbl.atac.sub[, c("chrom", "start", "end")],
                      matchThreshold=1e-3,
                      genomeName="hg38",
                      pwmFile=meme.file)
dim(tbl.fimo)
head(tbl.fimo)
tail(as.data.frame(sort(table(tbl.fimo$tf))))

library(trena)
candidate.tfs <- unique(tbl.fimo$tf)
length(candidate.tfs)   # 535
mtx <- get(load("~/github/TrenaProjectAD/inst/extdata/expression/mayo.tcx.16969x262.covariateCorrection.log+scale.RData"))
dim(mtx)   # 16969 262
length(setdiff(c(candidate.tfs, targetGene), rownames(mtx))) # 17/80
candidate.tfs <- intersect(candidate.tfs, rownames(mtx))
length(candidate.tfs)  # 370

solver <- EnsembleSolver(mtx,
                         targetGene=targetGene,
                         candidateRegulators=candidate.tfs,
                         solverNames=c("lasso", "Ridge", "Spearman", "Pearson", "RandomForest", "xgboost"),
                         geneCutoff=1.0)
tbl.out <- run(solver)
tbl.out <- tbl.out[order(abs(tbl.out$spearmanCoeff), decreasing=TRUE),]
dim(tbl.out)
head(tbl.out, n=10)

unlist(lapply(candidate.tfs, function(tf) length(grep(tf, tbl.fimo$tf))))
dim(tbl.out)

tbl.out$bs.3 <- unlist(lapply(tbl.out$gene,
                              function(TF)
                                  nrow(subset(tbl.fimo, tf==TF & p.value <= 1e-3))))
tbl.out$bs.4 <- unlist(lapply(tbl.out$gene,
                              function(TF)
                                  nrow(subset(tbl.fimo, tf==TF & p.value <= 1e-4))))

tbl.out$bs.5 <- unlist(lapply(tbl.out$gene,
                              function(TF)
                                  nrow(subset(tbl.fimo, tf==TF & p.value < 1e-5))))
tbl.out$bs.6 <- unlist(lapply(tbl.out$gene,
                              function(TF)
                                  nrow(subset(tbl.fimo, tf==TF & p.value < 1e-6))))

tbl.out.12 <- head(tbl.out, n=12)

tfs <- tbl.out.12$gene

for(TF in tfs){
    tbl.fimo.sub <- subset(tbl.fimo, tf==TF)
    track <- DataFrameAnnotationTrack(TF, tbl.fimo.sub[, c("chrom", "start", "end")], color="random",
                                      trackHeight=25)
    displayTrack(igv, track)
    }



