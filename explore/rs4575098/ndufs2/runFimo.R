library(ghdb)
source("~/github/fimoService/batchMode/fimoBatchTools.R")

targetGene <- "NDUFS2"

meme.file <- "jaspar2018-hocomocoCoreA.meme"
motifs <- query(MotifDb, c("sapiens"), c("jaspar2018", "HOCOMOCOv11-core-A"))

length(motifs)
export(motifs, con=meme.file, format="meme")

# NDUFS2 promoter + open atac
smallRegion <- function(){
    data.frame(chrom="chr1", start=161201350, end=161203624, stringsAsFactors=FALSE)
    }

bigRegion <- function(){
    ghdb <- GeneHancerDB()
    tbl.gh <- retrieveEnhancersFromDatabase(ghdb, targetGene, tissues="all")
    chromosome <- tbl.gh$chrom[1]
    start <- min(tbl.gh$start) - 1000
    end <- max(tbl.gh$end) + 1000
    data.frame(chrom=chromosome, start=start, end=end, stringsAsFactors=FALSE)
    }

tbl.regions <- bigRegion()
with(tbl.regions, printf("span: %dkb", round((end-start)/1000)))

system.time(tbl.fimo <- fimoBatch(tbl.regions,
                                  matchThreshold=1e-4,
                                  genomeName="hg38", pwmFile=meme.file))
save(tbl.fimo, file="ndufs2-fimo4-378kb-continuous.RData")

