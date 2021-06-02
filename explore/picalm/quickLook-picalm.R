igv <- start.igv("PICALM", "hg38")
snp.primary <- "rs10792832"
  # 11:86156833 (GRCh38)
  # 11:85867875 (GRCh37)
track <- DataFrameAnnotationTrack(snp.primary, data.frame(chrom="chr11", start=86156832, end=86156833,
                                                          stringsAsFactors=FALSE), color="darkred",
                                  trackHeight=23)
displayTrack(igv, track)

   # some table rows are ill-formed, but only in later columns.  use just the first 12:

tbl.snps <- read.table("haploreg-r2gt0.4-results.txt", sep="\t", as.is=TRUE,
                       header=TRUE, quote="", nrow=-1, fill=TRUE)[, 1:12]
dim(tbl.snps)
head(tbl.snps, n=4)
   #   chr pos_hg38   r2 Dprime is_query_snp       rsID ref       alt  AFR  AMR  ASN  EUR
   # 1  11 85959342 0.44   0.72            0 rs68172244 CAA         C 0.82 0.61 0.58 0.65
   # 2  11 85959902 0.45   0.74            0   rs673751   A         C 0.91 0.62 0.58 0.66
   # 3  11 85960209 0.42   0.67            0 rs67911082   A AATCT,ATC 0.88 0.62 0.56 0.63
   # 4  11 85970540 0.47   0.78            0   rs639012   A         G 0.90 0.62 0.60 0.67
threshold <- 0.75
tbl.snps.LD <- subset(tbl.snps, r2 >= threshold)
tbl.track <- tbl.snps.LD[, c("chr", "pos_hg38", "pos_hg38", "r2")]
colnames(tbl.track) <- c("chrom", "start", "end", "score")
tbl.track$chrom <- sprintf("chr%d", tbl.track$chrom)
tbl.track$start <- tbl.track$start - 1
track <- DataFrameQuantitativeTrack("LD.snps", tbl.track, autoscale=TRUE, color="maroon")
displayTrack(igv, track)

library(SNPlocs.Hsapiens.dbSNP151.GRCh38)
library("BSgenome.Hsapiens.UCSC.hg38")
library(motifbreakR)

snps.gr <- snps.from.rsid(rsid = tbl.snps.LD$rsID,
                          dbSNP=SNPlocs.Hsapiens.dbSNP151.GRCh38,
                          search.genome=BSgenome.Hsapiens.UCSC.hg38)

motifs <- query(MotifDb, c("sapiens"), c("jaspar2018", "HOCOMOCOv11-core"))

results <- motifbreakR(snpList = snps.gr,
                       filterp = TRUE,
                       pwmList = motifs,
                       show.neutral=FALSE,
                       method = c("ic", "log", "notrans")[2],
                       bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                       BPPARAM = BiocParallel::bpparam(),
                       verbose=TRUE)
tbl.results <- as.data.frame(results, row.names=NULL)
tbl.2 <- calculatePvalue(results)
tbl.results <- tbl.2
save(tbl.results, file="tbl.motifbreakR.results.r2.75.RData")
