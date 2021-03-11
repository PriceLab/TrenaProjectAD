library(R6)
library(trena)
library(TrenaProjectErythropoiesis)
library(RPostgreSQL)
library(ghdb)
# library(MotifDb)
# motifs <- query(MotifDb, "sapiens", c("jaspar2018", "hocomocoV11A"))
# length(motifs)
# export(motifs, con="jasparHocomoco.meme", format="meme")

#----------------------------------------------------------------------------------------------------
Predictor = R6Class("Predictor",

    #--------------------------------------------------------------------------------
    private = list(targetGene = NULL,
                   geneRegDB = NULL,
                   tbl.rbp = NULL,
                   mtx.rna = NULL,
                   meme.filename = NULL,
                   tbl.gh = NULL,
                   tbl.fimo = NULL
                   ),

    #--------------------------------------------------------------------------------
    public = list(

        igv = NULL,

        initialize = function(targetGene){
            private$targetGene <- targetGene
            private$geneRegDB <- dbConnect(PostgreSQL(), user= "trena", password="trena",
                                           dbname="genereg2021", host="khaleesi")

            dir <- "~/github/TrenaProjectErythropoiesis/inst/extdata/expression"
            file <- "brandLabDifferentiationTimeCourse-27171x28-namesCorrected.RData"
            tbl.rna <- get(load(file.path(dir, file)))
            private$mtx.rna <- as.matrix(tbl.rna)
            private$meme.filename = "jasparHocomoco.meme"
            },

        getTargetGene = function(){
           private$targetGene
           },

        getRBP = function(){
           query <- sprintf("select * from rbp where target='%s'", private$targetGene)
           private$tbl.rbp <- dbGetQuery(private$geneRegDB, query)
           private$tbl.rbp
           },

        getExpressionMatrix = function(){
           invisible(private$mtx.rna)
           },

        vizRBP = function(min.score=13){
           stopifnot(nrow(private$tbl.rbp) > 0)
           if(is.null(self$igv))
              self$igv <- start.igv(private$targetGene, "hg38")
           tbl.sub <- subset(private$tbl.rbp, score >= min.score)
           rbps <- sort(unique(tbl.sub$gene))
           for(rbp in rbps){
               tbl.rbp <- subset(tbl.sub, gene==rbp)[, c("chrom", "start", "endpos", "celltype", "score")]
               track <- DataFrameAnnotationTrack(rbp, tbl.rbp, color="random")
               displayTrack(self$igv, track)
               }
           },

        discoverTfCandidates = function(p.value=1e-4, elite.gh.only=TRUE, gh.min.score=20){
           ghdb <- GeneHancerDB()
           tbl.gh <- retrieveEnhancersFromDatabase(ghdb, private$targetGene, tissues="all")
           if(elite.gh.only)
              tbl.gh <- subset(tbl.gh, elite)
           tbl.gh$chrom <- sprintf("chr%s", tbl.gh$chrom)
           tbl.gh$score <- asinh(tbl.gh$combinedscore)
           private$tbl.gh <- tbl.gh
           source("~/github/fimoService/batchMode/fimoBatchTools.R")
           private$tbl.fimo <- fimoBatch(tbl.gh[, c("chrom", "start", "end")],
                                         matchThreshold=p.value, genomeName="hg38",
                                         pwmFile=private$meme.filename)
           },

        getGH = function(){
           invisible(private$tbl.gh)
           },

        getFimo = function(){
           invisible(private$tbl.fimo)
           },

        buildModel = function(candidateSources=c("tf", "rbp", "both")){
            }

       ) # public
    ) # class
#--------------------------------------------------------------------------------

