library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
library(SNPlocs.Hsapiens.dbSNP151.GRCh38)
tbl.loci <- read.table("supplementaryTable8-fineMap-3823-variants.tsv", sep="\t", header=TRUE, as.is=TRUE)
dim(tbl.loci)
lapply(tbl.loci, class)
tokens <- strsplit(tbl.loci$Variant, "\\:")
chrom <- paste0("chr", unlist(lapply(tokens, "[", 1)))
length(chrom)
head(chrom)
hg19 <- as.integer(unlist(lapply(tokens, "[", 2)))
tbl.loci$chrom <- chrom
tbl.loci$hg19 <- hg19
colnames(tbl.loci)
new.names <- c("associatedLocus", "variant", "finemap.pip", "susie.pip", "susie.cs", "chrom", "hg19")
colnames(tbl.loci) <- new.names
coi <- c("associatedLocus", "variant", "chrom", "hg19", "finemap.pip",  "susie.pip",  "susie.cs")
tbl.loci <- tbl.loci[, coi]

gr.variants <- GRanges(paste(sub("chr", "", tbl.loci$chrom), tbl.loci$hg19, sep=":"))
gr.snps <- snpsByOverlaps(SNPlocs.Hsapiens.dbSNP144.GRCh37, gr.variants)
tbl.rsids <- as.data.frame(gr.snps)
dim(tbl.rsids)   # 3603/3822
dim(tbl.loci)    # 3822
colnames(tbl.rsids) <- c("chrom", "hg19", "strand", "rsid", "iupac")
tbl.rsids$chrom <- as.character(tbl.rsids$chrom)
tbl.rsids$chrom <- paste0("chr", tbl.rsids$chrom)
head(tbl.rsids)

tbl.rsids$signature <- paste(tbl.rsids$chrom, tbl.rsids$hg19, sep=":")
tbl.loci$signature <- paste(tbl.loci$chrom, tbl.loci$hg19, sep=":")
dim(tbl.rsids)
dim(tbl.loci)

tbl.assoc <- unique(merge(tbl.loci, tbl.rsids, by="signature", all.x=TRUE))
dim(tbl.assoc)  # 3828 13   6 variant locations have
colnames(tbl.assoc)
coi <- c("associatedLocus", "chrom.x", "hg19.x", "variant", "rsid", "iupac", "finemap.pip", "susie.pip", "susie.cs")
tbl.assoc <- tbl.assoc[, coi]
colnames(tbl.assoc) <- c("associatedLocus", "chrom", "hg19", "variant", "rsid", "iupac", "finemap.pip",
                         "susie.pip", "susie.cs")
length(which(is.na(tbl.assoc$iupac)))  # [1] 225


variants.wo.rsid <- which(is.na(tbl.assoc$iupac))

tbl.with <- tbl.assoc[-variants.wo.rsid,]
tbl.wo   <- tbl.assoc[variants.wo.rsid,]
nrow(tbl.with) + nrow(tbl.wo) == nrow(tbl.assoc)

    #----------------------------------------
    #  get hg38 coords for tbl.with
    #----------------------------------------

gr.hg38 <- snpsById(SNPlocs.Hsapiens.dbSNP151.GRCh38, tbl.with$rsid, ifnotfound="drop")
tbl.hg38 <- as.data.frame(gr.hg38)
dim(tbl.with)
dim(tbl.hg38)   # dropped 7 rsids
colnames(tbl.hg38)[2] <- "hg38"
length(which(duplicated(tbl.hg38$hg38)))
tbl.with.2 <- merge(tbl.with, tbl.hg38, by.x="rsid", by.y="RefSNP_id", all.x=TRUE)
f <- function(code){
    if(code %in% names(IUPAC_CODE_MAP)){
        return(IUPAC_CODE_MAP[[code]])
    } else {
        return(code)
    }
}

alleles <- unlist(lapply(tbl.with.2$alleles_as_ambig, f))
tbl.with.2$alleles <- alleles
coi <- c("associatedLocus", "chrom", "hg19", "hg38", "variant", "rsid",
         "iupac", "finemap.pip", "susie.pip", "susie.cs", "alleles")
setdiff(coi, colnames(tbl.with.2))
tbl.with.3 <- tbl.with.2[, coi]

    #----------------------------------------
    # try liftOver
    #----------------------------------------

system("curl -O http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz")
system("gunzip hg19ToHg38.over.chain.gz")

chain <- import.chain("hg19ToHg38.over.chain")
gr.wo <- GRanges(paste(tbl.wo$chrom, tbl.wo$hg19, sep=":"))
gr.wo$variant <- tbl.wo$variant

x <- liftOver(gr.wo, chain)
gr.hg38 <- unlist(x)
tbl.wo.hg38 <- as.data.frame(gr.hg38)
colnames(tbl.wo.hg38)[2] <- "hg38"
length(gr.wo)
nrow(tbl.wo.hg38)

tbl.wo.2 <- merge(tbl.wo, tbl.wo.hg38, by="variant", all.x=TRUE)
colnames(tbl.wo.2)

coi <-  c("associatedLocus", "chrom", "hg19", "hg38", "variant", "rsid", "iupac", "finemap.pip",
          "susie.pip", "susie.cs")
tbl.wo.3 <- tbl.wo.2[, coi]

tokens <- strsplit(tbl.wo.3$variant, ":")
alleles.raw <- unlist(lapply(tokens, "[", 3))
alleles <- gsub("_", "", alleles.raw)
tbl.wo.3$alleles <- alleles

tbl.wo$rsid <- tbl.wo$variant

length(variants.wo.rsid)   # 225

tbl.$alleles <- unlist(lapply(tbl$alleles_as_ambig, function(code) IUPAC_CODE_MAP[[code]]))

tbl.assoc.1 <- subset(tbl.assoc, is.character(rsid))
dim(tbl.assoc.1)

    #----------------------------------------------
    # now combine & check
    #---------------------------------------------

tbl.assoc <- rbind(tbl.with.3, tbl.wo.3)
dim(tbl.assoc)
coi <- c("chrom", "hg19", "hg38", "variant", "rsid", "iupac", "finemap.pip", "susie.pip", "susie.cs", "alleles", "associatedLocus")
tbl.assoc <- tbl.assoc[, coi]

associated.gene <- tbl.genes$Gene[tbl.assoc$associatedLocus]
tbl.assoc$gene <- associated.gene

tbl.assoc <- tbl.assoc[order(tbl.assoc$associatedLocus, decreasing=FALSE),]
rownames(tbl.assoc) <- NULL


save(tbl.assoc, file="tbl.posthuma-38-geneAssociations-curated-3828x12.RData")




