source("~/github/fimoService/batchMode/fimoBatchTools.R")

targetGene <- "NDUFS2"
mtx.rna <- get(load("~/github/TrenaProjectAD/inst/extdata/expression/rosmap.14235x632.RData"))
motifs.oi <- query(MotifDb, "sapiens", c("jaspar2018", "hocomocov11-core-a"))
tfs <- sort(unique(mcols(motifs.oi)$geneSymbol))
tfs <- intersect(tfs, rownames(mtx.rna))
length(tfs) # 280
tf.cors <- lapply(tfs, function(tf) cor(mtx.rna[targetGene,], mtx.rna[tf,]))
names(tf.cors) <- tfs
fivenum(abs(as.numeric(tf.cors)))
keepers <- which(abs(as.numeric(tf.cors)) > 0.20) # 81
tfs.cor <- tfs[keepers]
length(tfs.cor)   # 81
tbl.match <- as.data.frame(findMatches(tfs.cor, mcols(motifs.oi)$geneSymbol))
keepers <- unique(tbl.match[,2]) # 122
length(keepers)
motifs.oi.cor <- motifs.oi[keepers]
length(motifs.oi.cor)  # 122
meme.file <- "ndufs2-cor-jaspar2018-hocomocoCoreA.meme"
export(motifs.oi.cor, con=meme.file, format="meme")

# NDUFS2 promoter + open atac
smallRegion <- function(){
    data.frame(chrom="chr1", start=161201350, end=161203624, stringsAsFactors=FALSE)
    }

nearbyRegion <- function(){
    data.frame(chrom="chr1", start=161180064, end=161230686, stringsAsFactors=FALSE)
    }

regRegion <- function(){
    data.frame(chrom="chr1", start=161179568, end=161243370, stringsAsFactors=FALSE)
    }

bigRegion <- function(){
    ghdb <- GeneHancerDB()
    tbl.gh <- retrieveEnhancersFromDatabase(ghdb, targetGene, tissues="all")
    chromosome <- tbl.gh$chrom[1]
    start <- min(tbl.gh$start) - 1000
    end <- max(tbl.gh$end) + 1000
    data.frame(chrom=chromosome, start=start, end=end, stringsAsFactors=FALSE)
    }

tbl.regions <- nearbyRegion()
tbl.regions <- smallRegion()
tbl.regions <- regRegion()
with(tbl.regions, printf("span: %dkb", round((end-start)/1000)))

system.time(tbl.fimo <- fimoBatch(tbl.regions,
                                  matchThreshold=1e-3,
                                  genomeName="hg38", pwmFile=meme.file))
save(tbl.fimo, file="ndufs2-fimo4-64kb-continuous.RData")


