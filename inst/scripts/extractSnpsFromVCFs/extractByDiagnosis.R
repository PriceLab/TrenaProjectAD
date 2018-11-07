library(TrenaProjectIGAP)
library(RUnit)
library(GenomicRanges)
library(VariantAnnotation)
#----------------------------------------------------------------------------------------------------
igap <- TrenaProjectIGAP()
goi <- getSupportedGenes(igap)
load(system.file(package="TrenaProjectIGAP", "extdata", "misc", "tbl.enhancerRegions.RData"))
dim(tbl.enhancerRegions)
length(goi)

tbl.covar <- getCovariatesTable(igap)
#ad.samples <- sub("_TCX", "", subset(tbl.covar, Diagnosis=="AD")$ID)
#ctl.samples <- sub("_TCX", "", subset(tbl.covar, Diagnosis=="Control")$ID)
load(system.file(package="TrenaProjectIGAP", "extdata", "misc", "tbl.enhancerRegions.RData"))
dim(tbl.enhancerRegions)
     
vcf.dir <- "~/s/data/sage/ad-mayo"
vcf.file.template <- "NIA_JG_1898_samples_GRM_WGS_b37_JointAnalysis01_2017-12-08_%s.recalibrated_variants.Mayo.vcf.gz"

#----------------------------------------------------------------------------------------------------
getSnpTable <- function(gene, diagnosis)
{
   stopifnot(gene %in% rownames(tbl.enhancerRegions))
   chromosome <- tbl.enhancerRegions[gene,"hg19.chrom"]
   chromosome.stripped <- sub("chr", "", chromosome)
   loc.start  <- tbl.enhancerRegions[gene,"hg19.start"]
   loc.end  <- tbl.enhancerRegions[gene,"hg19.end"]
   gr <- GRanges(chromosome.stripped, IRanges(loc.start, loc.end))
   sampleIDs <-  switch(diagnosis,
                        AD = sub("_TCX", "", subset(tbl.covar, Diagnosis=="AD")$ID),
                        Control = sub("_TCX", "", subset(tbl.covar, Diagnosis=="Control")$ID))
   filename <- sprintf(vcf.file.template, chromosome.stripped)
   full.path <- file.path(vcf.dir, filename)
   stopifnot(file.exists(full.path))
   params <- ScanVcfParam(which=gr, samples=sampleIDs)
   suppressWarnings(
      vcf <- expand(readVcf(full.path,  "hg19", params))  # samples not in file reported
      )
   dim(vcf)
   tbl <- vcfToTable(vcf)
   tbl$gene <- gene

} # getSnpTable
#----------------------------------------------------------------------------------------------------
test_getSnpTable <- function()
{
    printf("--- test_getSnpTable")
    x <- getSnpTable("INPP5D", "AD")

} # test_getSnpTable
#----------------------------------------------------------------------------------------------------
vcfToTable <- function(vcf)
{
    # only want the locations with non-NA AF (allele frequency), AC (allele count)
  # browser()
  deleters <- which(is.na(unlist(info(vcf)[["AF"]])))
  keepers <- which(!is.na(unlist(info(vcf)[["AF"]])))
  length(deleters); range(deleters)
  length(keepers); range(keepers)
  
  tbl.info <- as.data.frame(info(vcf))[keepers,]
  dim(tbl.info)
 
  mtx.geno <- geno(vcf)$GT[keepers,]
  dim(mtx.geno)
  mtx.012 <- matrix(0, nrow=nrow(mtx.geno), ncol=ncol(mtx.geno), dimnames=list(rownames(mtx.geno), colnames(mtx.geno)))
  mtx.012[which(mtx.geno=="0/1")] <- 1
  mtx.012[which(mtx.geno=="1/1")] <- 2
  mtx.geno <- mtx.012
  tbl.geno <- as.data.frame(mtx.geno)
  dim(tbl.geno)
  head(lapply(tbl.geno, class))

  tbl.loc <- as.data.frame(rowRanges(vcf[keepers]))
  dim(tbl.loc)
  tbl <- cbind(cbind(tbl.loc, tbl.info), tbl.geno)
  dim(tbl)  #  18405 383      (349 + 24 + 10)
  colnames(tbl)[1] <- "chr"
  tbl$chr <- paste("chr", tbl$chr, sep="")
  invisible(tbl)    

} # vcfToTable
#----------------------------------------------------------------------------------------------------
test_vcfToTable <- function()
{
   if(!exists("vcf.demo"))
      load("inpp5d.AD.vcf.RData")

   tbl <- vcfToTable(expand(vcf.demo))
   checkEquals(dim(tbl), c(8215, 113))
   
} # test_vcfToTable
#----------------------------------------------------------------------------------------------------
getAllTables <- function()
{
  x.ad <- lapply(goi, function(gene) getSnpTable(gene, "AD"))
  x.ctl <- lapply(goi, function(gene) getSnpTable(gene, "Control"))

} # getAllTables
#----------------------------------------------------------------------------------------------------
