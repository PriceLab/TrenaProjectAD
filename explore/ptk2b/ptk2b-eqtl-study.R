library(RUnit)
library(EndophenotypeExplorer)
library(TrenaMultiScore)
library(TrenaProjectAD)
source("~/github/TrenaMultiScore/tools/runner/v2/tmsCore.R")

tbl.fimo <- get(load("tbl.fimo.PTK2B.RData"))
tbl.atac <- get(load("~/github/TrenaProjectAD/explore/mayo-epigenetics/atac/mayoAllPeaks.1052789x4.RData"))
head(tbl.atac)
tpad <- TrenaProjectAD()
targetGene <- "PTK2B"

tbl.posthuma <- get(load("../../inst/extdata/gwasLoci/tbl.posthuma-38-loci-curated.RData"))
tbl.wms <- get(load("../../inst/extdata/gwasLoci/williams-natureNeuroscience2020.RData"))
subset(tbl.wms, locusOrGene=="PTK2B")
   #   locusOrGene       rsid
   # 24       PTK2B rs28834970

etx <- EndophenotypeExplorer$new(targetGene, "hg38")
tbl.eQTL <- etx$getEQTLsForGene()
tbl.eQTL.sig <- subset(tbl.eQTL, pvalue < 0.001 & study=="ampad-mayo")
new.order <- order(tbl.eQTL.sig$pvalue, decreasing=FALSE)
tbl.eQTL.sig <- tbl.eQTL.sig[new.order,]
rsids <- tbl.eQTL.sig$rsid
length(rsids)


subset(tbl.eQTL, rsid=="rs28834970")
#       chrom     hg19     hg38       rsid    pvalue            ensg genesymbol        study tissue
# 1617   chr8 27195121 27337604 rs28834970 0.7766134 ENSG00000120899      PTK2B   ampad-mayo    tcx
# 7207   chr8 27195121 27337604 rs28834970 0.1501057 ENSG00000120899      PTK2B   ampad-mayo    cer
# 12792  chr8 27195121 27337604 rs28834970 0.5057638 ENSG00000120899      PTK2B ampad-rosmap  dlpfc
#
# haploreg LD >= 0.8
# 8	27337604	1	1	rs28834970
# 8	27347529	0.98	0.99	rs57735330
# 8	27350609	0.99	1	rs6987305
# 8	27354393	0.97	0.99	rs2322599
# 8	27362470	0.97	0.99	rs73223431
# 8	27362793	0.95	0.99	rs17057043

tag.snp <- "rs28834970"
ld.snps <- c("rs57735330", "rs6987305", "rs2322599", "rs73223431", "rs17057043")
coi <- c("chrom", "hg38", "rsid", "pvalue")
subset(tbl.eQTL, rsid %in% c(tag.snp, ld.snps) & study=="ampad-mayo" & tissue == "cer")[, coi]
  # all at about 0.16, hardly signficiant.


tms <- TMS$new(tpad, targetGene, tbl.fimo=tbl.fimo, tbl.oc=tbl.atac)
tms$addGeneHancer()
tms$scoreFimoTFBS()
tbl.tms <- tms$getTfTable()
dim(tbl.tms)   # 51491 17
table(tbl.tms$chip)   # 46393 5098
table(tbl.tms$oc)     # 38081 13410

tbl.strong.01 <- subset(tbl.tms, oc & gh > 600 & (chip | fimo_pvalue < 1e-6))
tfs.01 <- unique(tbl.strong.01$tf)
length(tfs.01)

tbl.strong.02 <- subset(tbl.tms, oc & (chip)) # | fimo_pvalue < 1e-3))
tfs.02 <- unique(tbl.strong.02$tf)
length(tfs.02)
length(intersect(tfs.01, tfs.02))

tbl.strong.03  <- subset(tbl.tms, oc  & (chip | fimo_pvalue < 1e-3))
tfs.03 <- unique(tbl.strong.03$tf)
length(tfs.03)
length(intersect(tfs.01, tfs.03))

if(exists("igv")){
   track <- DataFrameAnnotationTrack("tf", tbl.strong.03[, c("chrom", "start", "end")], color="purple")
   displayTrack(igv, track)
   }


tbl.tfs <- as.data.frame(sort(table(tbl.strong.03$tf), decreasing=TRUE))
subset(tbl.tfs, Freq > 1)

tfs <- subset(tbl.tfs, Freq > 4)$Var1
length(tfs)   # 70
dir <- "~/github/TrenaProjectAD/prep/rna-seq-counts-from-synapse/eqtl"
file.cer <- "mtx.mayo.cer.eqtl-optimized-geneSymbols-sampleIDs-with-vcf17009x255.RData"
mtx.cer <- get(load(file.path(dir, file.cer)))

#------------------------------------------------------------------------------------------------------------------
select.rsids <- function()
{
     #------------------------------------------------------------
     # intersect atac with eqtls
     #------------------------------------------------------------

   gr.eQTL <- GRanges(seqnames=tbl.eQTL$chrom, IRanges(start=tbl.eQTL$hg38-1, end=tbl.eQTL$hg38))
   tbl.ov <- as.data.frame(findOverlaps(GRanges(tbl.oc), gr.eQTL))
   dim(tbl.ov)
   tbl.eQTL.atac <- tbl.eQTL[tbl.ov$subjectHits,]
   tbl.eQTL.atac.sig <- subset(tbl.eQTL.atac, pvalue < 0.05 & study=="ampad-mayo" & tissue=="cer")
   dim(tbl.eQTL.atac.sig)
   eQTLS <- unique(tbl.eQTL.atac.sig$rsid)
   length(eQTLS)
   eQTLS
}
#------------------------------------------------------------------------------------------------------------------
tbl.oc <- tms$getOpenChromatin()
tbl.ghr <- tms$getGeneHancerRegion()

tbl.oc <- subset(tbl.oc, chrom==tbl.ghr$chrom & start >= tbl.ghr$start & end <= tbl.ghr$end)
dim(tbl.oc)

rsids <- select.rsids()
length(rsids)
length(tfs.03)

result <- list()
rsids <- "rs11780653"
for(rsid in rsids){
   printf("--- %s", rsid)
   #rsid <- "rs28834970"  # the tag.snp
   #rsid <- "rs10087271"  # mayo cer eQTL pval = 0.806
   #rsid <- "rs10503792"  # pval 0.50
   x <- etx$splitExpressionMatrixByMutationStatusAtRSID(mtx.cer, rsid, study.name="mayo")
   adequate.sample.balance <- x$genotypes$wt/ncol(mtx.cer) > 0.2 & x$genotypes$mut/ncol(mtx.cer) > 0.2
   if(!adequate.sample.balance){
      printf("rsid %s poorly distributed in rna, %d wt, %d mut", rsid, x$genotypes$wt, x$genotypes$mut)
      next;
      }
   xx <- etx$trenaScoreGenotypeStratifiedExpression(x$wt, x$mut, targetGene, tfs.03)
   print(head(xx$rf.delta, n=10))
   print(head(xx$bicor.delta, n=10))
   trena.wt.rank <- lapply(xx$bicor.delta$tf, function(tf) grep(sprintf("^%s$", tf), xx$trena.1$gene))
   failures <- which(unlist(lapply(trena.wt.rank, length))==0)
   trena.wt.rank[failures] <- -1
   trena.wt.rank <- unlist(trena.wt.rank)

   trena.mut.rank <- lapply(xx$bicor.delta$tf, function(tf) grep(sprintf("^%s$", tf), xx$trena.2$gene))
   failures <- which(unlist(lapply(trena.mut.rank, length))==0)
   trena.mut.rank[failures] <- -1
   trena.mut.rank <- unlist(trena.mut.rank)

   xx$bicor.delta$trena.wt.rank <- trena.wt.rank
   xx$bicor.delta$trena.mut.rank <- trena.mut.rank
   tbl.delta.summary <- subset(xx$bicor.delta,
                              ((trena.wt.rank > 0 & trena.wt.rank <= 10) | (trena.mut.rank > 0 & trena.mut.rank <= 10)) &
                              abs(delta) > 0.2)
   if(nrow(tbl.delta.summary) > 0)
       tbl.delta.summary$rsid <- rsid
   xx$delta.summary <- tbl.delta.summary
   result[[rsid]] <- list(matrices=x, models=xx)
   }

lapply(result, function(item) item$models$delta.summary)

#------------------------------------------------------------------------------------------------------------------------
# $rs2292974
#        tf method         wt        mut      delta trena.wt.rank trena.mut.rank      rsid
# 4    RARA  bicor  0.1815432  0.5210778  0.3395345            79              1 rs2292974
# 7   PROX1  bicor -0.4131147 -0.1014805  0.3116342             9            177 rs2292974
# 9   MYBL1  bicor -0.5283956 -0.2267122  0.3016834             1             49 rs2292974
# 30 ZNF263  bicor -0.4755963 -0.2200591  0.2555371             7             77 rs2292974
# 36  SCRT1  bicor  0.5154509  0.2730314 -0.2424195             3             38 rs2292974
#
#    rs10105298:  wt mut het hom
#                   85 170 132  38

forMax <- function()
{

   rsids <- c("rs35298960",   #0.4875633
              "rs9797", #0.3364394
              "rs17426801", #0.2117080
              "rs328071", #0.1673497
              "rs17404686", #0.5927874
              "rs78945847", #0.6459213
              "rs11987934", #0.4034604
              "rs12676581", #0.3352473
              "rs11780653", #0.9567037
              "rs939267",
              "rs1560344",
              "rs11780653",
) #0.4086222


   for(rsid in rsids){
      x <- etx$splitExpressionMatrixByMutationStatusAtRSID(mtx.cer, rsid, study.name="mayo")
      printf("--- %s", rsid)
      print(unlist(x$genotypes))
      }

   #rsid.oi <- "rs10105298"
                                          #  141 114 102  12

    # $rs11780653
    #      tf method         wt        mut      delta trena.wt.rank trena.mut.rank       rsid
    #2  FOXO3  bicor  0.4555155  0.1745318 -0.2809837             2             55 rs11780653
    #15  E2F4  bicor -0.1569254 -0.3590907 -0.2021653           108              8 rs11780653


   mtx.wt <- result[[rsid.oi]]$matrices$wt
   mtx.mut <- result[[rsid.oi]]$matrices$mut
   mtx.hom <- result[[rsid.oi]]$matrices$hom
   mtx.het <- result[[rsid.oi]]$matrices$het
   dim(mtx.hom)   # 12 for rs11780653
   dim(mtx.hom)   # 59 for rs1560344

   rara.wt <- mtx.wt["RARA",]
   ptk2b.wt <- mtx.wt["PTK2B",]
   rara.mut <- mtx.mut["RARA",]
   ptk2b.mut <- mtx.mut["PTK2B",]

   motifs <- query(MotifDb, c("sapiens"), c("jaspar2018", "HOCOMOCOv11-core"))
   geneSymbols <- unique(mcols(motifs)$geneSymbol)  # 692
   goi <- c("PTK2B", intersect(rownames(mtx.wt), geneSymbols))
   length(goi) # 346

   mtx.wt.sub <- mtx.wt[goi,]
   mtx.mut.sub <- mtx.mut[goi,]
   mtx.hom.sub <- mtx.hom[goi,]
   mtx.het.sub <- mtx.het[goi,]

   cor(mtx.mut.sub["FOXO3",], mtx.mut.sub["PTK2B",], use="pairwise.complete", method="spearman")  # 0.408
   cor(mtx.wt.sub["FOXO3",], mtx.wt.sub["PTK2B",], use="pairwise.complete", method="spearman")  # 0.408
   save(mtx.wt.sub, mtx.mut.sub, mtx.hom.sub, mtx.het.sub, file="mayo.cer.rnaseq.rs11780653-separated-346genes.RData")
   save(mtx.wt.sub, mtx.mut.sub, mtx.hom.sub, mtx.het.sub, file="mayo.cer.rnaseq.rs1560344-separated-346genes.RData")

   cor(rara.mut, ptk2b.mut, use="pairwise.complete", method="spearman")  # 0.408
   cor(rara.wt, ptk2b.wt, use="pairwise.complete", method="spearman")  # 0.408
   t.test(ptk2b.wt, ptk2b.mut)$p.value   # 0.0288
   plot(rara.mut, ptk2b.mut)

} # forMax
#------------------------------------------------------------------------------------------------------------------------

checkEquals(sort(names(x)), c("genotypes", "het", "hom", "mut", "wt"))
checkEquals(x$genotypes, list(wt=95, mut=160, het=124, hom=36))

print(load("~/github/TrenaProjectAD/explore/ptk2b/tbl.tms.RData"))

x <- etx$trenaScoreGenotypeStratifiedExpression(x$wt, x$mut, targetGene, tfs)
head(x$rf.delta, n=10)
head(x$bicor.delta, n=10)

#----------------------------------------------------------------------------------------------------
viz <- function()
{
   igv <- start.igv("PTK2B")
   zoomOut(igv);   zoomOut(igv);
   tbl.gh <- tms$getGeneHancer()
   tbl.gh$combinedscore <- asinh(tbl.gh$combinedscore)
   track <- DataFrameQuantitativeTrack("GH", tbl.gh[, c("chrom", "start", "end", "combinedscore")],
                                       color="brown", autoscale=TRUE)
   displayTrack(igv, track)

   tbl.eQTL <- etx$getEQTLsForGene()
   tbl.eQTL.sig <- subset(tbl.eQTL, pvalue < 0.01 & study=="ampad-mayo")
   tbl.eQTL.sig <- subset(tbl.eQTL, pvalue < 0.001 & study=="ampad-mayo" & tissue=="cer")
   dim(tbl.eQTL.sig)
   setdiff(c(tag.snp, ld.snps), tbl.eQTL.sig$rsid)   # "rs57735330"
   tbl.track <- tbl.eQTL.sig[, c("chrom", "hg38", "hg38", "pvalue", "rsid")]
   colnames(tbl.track) <- c("chrom", "start", "end", "score", "rsid")
   tbl.track$start <- tbl.track$start - 1
   tbl.track$score <- -log10(tbl.track$score)
   track <- DataFrameQuantitativeTrack("eQTL < 0.001", tbl.track, color="red", autoscale=TRUE)
   displayTrack(igv, track)
   tbl.tag.snp <- data.frame(chrom="chr8", start=27337604, end=27337604, name="rs28834970", stringsAsFactors=FALSE)
   track <- DataFrameAnnotationTrack("GWAS tag snp", tbl.tag.snp, color="black")
   displayTrack(igv, track)

   tbl.ghr <- tms$getGeneHancerRegion()
   tbl.oc <- tms$getOpenChromatin()
   tbl.oc <- subset(tbl.oc, chrom==tbl.ghr$chrom & start >= tbl.ghr$start & end <= tbl.ghr$end)
   dim(tbl.oc)
   track <- DataFrameQuantitativeTrack("ATAC", tbl.oc, autoscale=TRUE, color="darkblue")
   displayTrack(igv, track)


   tbl.track <- tbl.eQTL.atac.sig[, c("chrom", "hg38", "hg38", "pvalue", "rsid")]
   colnames(tbl.track) <- c("chrom", "start", "end", "score", "rsid")
   tbl.track$start <- tbl.track$start - 1
   tbl.track$score <- -log10(tbl.track$score)
   dups <- which(duplicated(tbl.track$rsid))
   tbl.track <- tbl.track[-dups,]
   track <- DataFrameQuantitativeTrack("eQTL atac < 0.05", tbl.track, color="red", autoscale=TRUE)
   displayTrack(igv, track)
   dim(tbl.eQTL.atac)

   tbl.tead1 <- subset(tbl.strong.03, tf=="TEAD1")[, c("chrom", "start", "end", "fimo_pvalue")]
   colnames(tbl.tead1)[4] <- "score"
   tbl.tead1$score <- -log10(tbl.tead1$score)
   track <- DataFrameQuantitativeTrack("TEAD1", tbl.tead1, color="blue", autoscale=TRUE)
   displayTrack(igv, track)

   tbl.neurod1 <- subset(tbl.strong.03, tf=="NEUROD1")[, c("chrom", "start", "end", "fimo_pvalue")]
   colnames(tbl.neurod1)[4] <- "score"
   tbl.neurod1$score <- -log10(tbl.neurod1$score)
   track <- DataFrameQuantitativeTrack("NEUROD1", tbl.neurod1, color="blue", autoscale=TRUE)
   displayTrack(igv, track)



} # viz
#----------------------------------------------------------------------------------------------------
