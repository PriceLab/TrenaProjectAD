tbl.genes <-read.table("posthuma-38-loci.csv", sep=",", header=TRUE, as.is=TRUE)
dim(tbl.genes)
lapply(tbl.genes, class)
head(tbl.genes)
chrom <- paste0("chr", unlist(lapply(strsplit(tbl.genes$hg19, ":"), "[", 1)))
tbl.genes$chrom <- chrom
hg19 <- as.numeric(sub("*.:", "", tbl.genes$hg19))
tbl.genes$hg19 <- hg19
# library(rtracklayer)
# gr.hg19 <- GRanges(seqnames=tbl.genes$chrom, IRanges(hg19))
#
# system("curl -O http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz")
# system("gunzip hg19ToHg38.over.chain.gz")
#
# chain <- import.chain("hg19ToHg38.over.chain")
# x <- liftOver(gr.hg19, chain)
# gr.hg38 <- unlist(x)
# tbl.hg38 <- as.data.frame(gr.hg38)
# tbl.genes$hg38 <- tbl.hg38$start
#
# colnames(tbl.genes)[grep("Gene", colnames(tbl.genes))] <- "sym.orig"
# source("~/github/projects/utils/geneIdMapping/symToGeneID.R"); test_assignGeneIDs()
# x <- assignGeneIDs(tbl.genes$sym.orig)
#
# mapIds(EnsDb.Hsapiens.v79, c("ENSG00000176022", "ENSG00000164008", "ENSG00000164008"), "SYMBOL", "GENEID")
# mapIds(EnsDb.Hsapiens.v79, c("B3GALT6", "C1orf50"), "GENEID", "SYMBOL")


gr.snps <- snpsById(SNPlocs.Hsapiens.dbSNP151.GRCh38, tbl.genes$leadVariant)
tbl.snps <- as.data.frame(gr.snps)

tbl <- unique(merge(tbl.genes, tbl.snps, by.x="leadVariant", by.y="RefSNP_id"))
tbl$seqnames <- as.integer(tbl$seqnames)
tbl <- tbl[order(tbl$seqnames, decreasing=FALSE),]
colnames(tbl)[grep("pos", colnames(tbl))] <- "hg38"
strand.column <- grep("strand", colnames(tbl))
tbl <- tbl[, -(strand.column)]

tbl$alleles <- unlist(lapply(tbl$alleles_as_ambig, function(code) IUPAC_CODE_MAP[[code]]))

ambig.column <- grep("alleles_as_ambig", colnames(tbl))
tbl <- tbl[, -(ambig.column)]
stopifnot(all(tbl$chrom == paste0("chr", tbl$seqnames)))

seqnames.column <- grep("seqnames", colnames(tbl))
tbl <- tbl[, -(seqnames.column)]

colnames(tbl)[grep("Gene", colnames(tbl))] <- "geneSymbol"

library(EnsDb.Hsapiens.v79)

ensg <- mapIds(EnsDb.Hsapiens.v79, tbl$geneSymbol, "GENEID", "SYMBOL")
stopifnot(!any(is.na(ensg)))
tbl$ensg <- as.character(ensg)
rownames(tbl) <- NULL
colnames(tbl)

coi <- c(
    "geneSymbol",
    "leadVariant",
    "chrom",
    "hg19",
    "hg38",
    "ensg",
    "alleles",
    "A1",
    "A1.frequency",
    "BETA",
    "SE",
    "P",
    "N")

tbl <- tbl[, coi]
dim(tbl)
save(tbl, file="tbl.posthuma-38-loci-curated.RData")
