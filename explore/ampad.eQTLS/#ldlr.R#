library(VariantAnnotation)
#----------------------------------------------------------------------------------------------------
get.eqtls <- function()
{
   dir <- "/proj/price4/cory/AMP-AD_eQTLs"
   files <- sort(dir(dir))
   names(files) <- c("mayo.cer", "rosmap.dlpfc", "mayo.tcx")
   tbls <- list()
   for(i in seq_len(length(files))){
      file <- as.character(files)[i]
      short.name <- names(files)[i]
      full.path <- file.path(dir, file)
      output.filename <- sprintf("%s.ldlr.csv", short.name)
      if(!file.exists(output.filename)){
         cmd <- sprintf("grep ',LDLR,' %s > %s", full.path, output.filename)
         system(cmd)
         }
      tbl <- read.table(output.filename, sep=",", as.is=TRUE, nrow=-1)
      colnames(tbl) <- c("chromosome","snpLocation","snpid","gene","geneSymbol","statistic","pvalue","FDR",
                         "beta","A1","A2","A2freq","expressionIncreasingAllele","strand","geneBiotype",
                         "geneStartPosition","geneEndPosition")
      tbl$source <- short.name
      tbls[[short.name]] <- tbl
      }

   tbl.eqtl <- do.call(rbind, tbls)
   rownames(tbl.eqtl) <- NULL
   filename <- "ldlr-eqtls.RData"
   printf("saving %d eqtls to %s", nrow(tbl.eqtl), filename)
   save(tbl.eqtl, file=filename)

} # get.eqtls
#----------------------------------------------------------------------------------------------------
getVariants <- function()
{
   regions <- GRanges(seqnames="19", IRanges(start=10200760, end=12242459))
   vcf.file <- "../ampad-1898-samples/vcf/NIA_JG_1898_samples_GRM_WGS_b37_JointAnalysis01_2017-12-08_19.recalibrated_variants.vcf.gz"
   file.exists(vcf.file)
   file.exists(sprintf("%s.tbi", vcf.file))
    
   vcf <- readVcf(vcf.file, "hg19", regions)
   length(vcf)
   dim(geno(vcf)$GT)  # 5 1894
    
   mtx.geno <- geno(vcf)$GT
   save(mtx.geno, file="ldlr-geno.RData")
   
    
    mtx.geno[1:10, 1:10]
    wdth(100)
    length(unique(colnames(mtx.geno)))
    vcf.sample.names <- colnames(mtx.geno)
    save(vcf.sample.names, file="vcf.1894.sample.names.RData")

} # getVariants
#----------------------------------------------------------------------------------------------------
intersectVariantsAndEQTLs <- function()
{
   if(!exists("mtx.geno"))
      mtx.geno <- get(load("ldlr-geno.RData"))

   if(!exists("tbl.eqtl"))
      tbl.eqtl <- get(load("ldlr-eqtls.RData"))

    tbl.geno <- as.data.frame(mtx.geno)

    geno.locs <- as.integer(sub("_.*$", "", sub("^.*:", "", rownames(tbl.geno))))
    tbl.geno$loc <- geno.locs
    dim(tbl.eqtl) # 9675

       #----------------------------------------------------------
       # one eqtl reported (by mayo?, rosmap?) not in the vcf file
       #----------------------------------------------------------

    setdiff(tbl.eqtl$snpLocation, tbl.geno$loc) # [1] 12242460
    subset(tbl.eqtl, snpLocation==12242460)
     #      chromosome snpLocation     snpid            gene geneSymbol   statistic    pvalue       FDR         beta A1 A2    A2freq expressionIncreasingAllele strand    geneBiotype geneStartPosition geneEndPosition       source
     # 3226         19    12242460 rs1036593 ENSG00000130164       LDLR -1.09826536 0.2731245 0.8975289 -0.143770746  A  T 0.1186700                          A      1 protein_coding          11089362        11133816     mayo.cer
     # 6449         19    12242460 rs1036593 ENSG00000130164       LDLR  0.79758140 0.4254479 0.9049935  0.069192134  A  T 0.1203762                          T      1 protein_coding          11089362        11133816 rosmap.dlpfc
     # 9675         19    12242460 rs1036593 ENSG00000130164       LDLR -0.05739204 0.9542777 0.9975668 -0.007479608  A  T 0.1186700                          A      1 protein_coding          11089362        11133816     mayo.tcx

       #----------------------------------------------------------
       # just 3299 distinct locations in eQTL & VCF
       #----------------------------------------------------------
      
    length(intersect(tbl.geno$loc, tbl.eqtl$snpLocation))   # 3299
    length(unique(tbl.eqtl$snpLocation))  # [1] 3300
    length(setdiff(tbl.eqtl$snpLocation, tbl.geno$loc))     # 

    tbl.geno.ldlr <- subset(tbl.geno, loc %in% tbl.eqtl$snpLocation)
    dim(tbl.geno.ldlr)
    save(tbl.geno.ldlr, file="tbl.geno.eqtl.ldlr.RData")

} # intersectVariantsAndEQTLs
#----------------------------------------------------------------------------------------------------
mapToPatientID <- function()
{
   if(!exists("tbl.geno.ldlr"))
      tbl.geno.ldlr <- get(load("tbl.geno.eqtl.ldlr.RData"))
   dim(tbl.geno.ldlr)

   ids <- colnames(tbl.geno.ldlr)
   length(ids)

   dir <- "~/github/TrenaProjectAD/explore/WGS-Harmonization/metadata-all"
   tbl.rosmap <- read.table(file.path(dir, "ROSMAP_biospecimen_metadata.csv"), sep=",", as.is=TRUE, header=TRUE)
   tbl.sinai <- read.table(file.path(dir, "MSBB_biospecimen_metadata.csv"), sep=",", as.is=TRUE, header=TRUE)
   tbl.mayo <- read.table(file.path(dir, "MayoRNAseq_biospecimen_metadata.csv"), sep=",", as.is=TRUE, header=TRUE)

   length(intersect(ids, tbl.rosmap$specimenID)) # 1151
   length(intersect(ids, tbl.sinai$specimenID))  #  349
   length(intersect(ids, tbl.mayo$specimenID))   #  349
                                                 # 1849

   vcf.ids.rosmap <- intersect(ids, tbl.rosmap$specimenID)
   ids.sinai <- intersect(ids, tbl.sinai$specimenID)
   ids.mayo <- intersect(ids, tbl.mayo$specimenID)

   length(intersect(ids.sinai, ids.mayo))       # though both have a 349 count, they do not overlap

   

   

} # mapToPatientID
#----------------------------------------------------------------------------------------------------


