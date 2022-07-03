library(VariantAnnotation)
library(RUnit)
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
   file.exists(dir)
   tbl.rosmap <- read.table(file.path(dir, "ROSMAP_biospecimen_metadata.csv"), sep=",", as.is=TRUE, header=TRUE)
   tbl.sinai <- read.table(file.path(dir, "MSBB_biospecimen_metadata.csv"), sep=",", as.is=TRUE, header=TRUE)
   tbl.mayo <- read.table(file.path(dir, "MayoRNAseq_biospecimen_metadata.csv"), sep=",", as.is=TRUE, header=TRUE)

   length(intersect(ids, tbl.rosmap$specimenID)) # 1151
   length(intersect(ids, tbl.sinai$specimenID))  #  349
   length(intersect(ids, tbl.mayo$specimenID))   #  349
                                                 # 1849

   vcf.ids.rosmap <- intersect(ids, tbl.rosmap$specimenID)
   vcf.ids.sinai <- intersect(ids, tbl.sinai$specimenID)
   vcf.ids.mayo <- intersect(ids, tbl.mayo$specimenID)

   length(intersect(ids.sinai, ids.mayo))       # though both have a 349 count, they do not overlap

   map.rosmap <- lapply(vcf.ids.rosmap, function(id) subset(tbl.rosmap, specimenID==id)$individualID[1])
   names(map.rosmap) <- vcf.ids.rosmap

   map.sinai <- lapply(vcf.ids.sinai, function(id) subset(tbl.sinai, specimenID==id)$individualID[1])
   names(map.sinai) <- vcf.ids.sinai

   map.mayo <- lapply(vcf.ids.mayo, function(id) subset(tbl.mayo, specimenID==id)$individualID[1])
   names(map.mayo) <- vcf.ids.mayo

   length(which(is.na(map.rosmap)))   # 0
   length(which(is.na(map.sinai)))    # 4
   sinai.failures <- as.integer(which(is.na(map.sinai)))
   map.sinai <- map.sinai[-sinai.failures]
   length(map.sinai)
   length(which(is.na(map.mayo)))    # 4

   length(map.rosmap) +    # 1141
   length(map.sinai) +     # 345
   length(map.mayo)        # 349
                           # 1845

      #------------------------------------------------------------------------
      # now use these maps to construct an id table: vcf, patient, project
      #------------------------------------------------------------------------

   tbl.id <- data.frame(vcf=ids, patient=NA, cohort=NA, stringsAsFactors=FALSE)
   rosmap.rows <- match(names(map.rosmap), tbl.id$vcf)
   sinai.rows  <- match(names(map.sinai),  tbl.id$vcf)
   mayo.rows   <- match(names(map.mayo),   tbl.id$vcf)

   stopifnot(length(intersect(rosmap.rows, sinai.rows)) == 0)
   stopifnot(length(intersect(rosmap.rows, mayo.rows)) == 0)
   stopifnot(length(intersect(mayo.rows, sinai.rows)) == 0)

   tbl.id$cohort[rosmap.rows] <- "rosmap"
   tbl.id$cohort[sinai.rows]  <- "sinai"
   tbl.id$cohort[mayo.rows]   <- "mayo"

   tbl.id$patient[rosmap.rows] <- as.character(map.rosmap)
   tbl.id$patient[sinai.rows] <- as.character(map.sinai)
   tbl.id$patient[mayo.rows] <- as.character(map.mayo)


      #------------------------------------------------------------------------
      # 50 vcf sample ids, plus my own artifact "loc", failed to map
      #------------------------------------------------------------------------
   length(which(is.na(tbl.id$cohort)))  # 50
   failures <- which(is.na(tbl.id$cohort))
   tbl.id <- tbl.id[-failures,]
   save(tbl.id, file="tbl.vcfToPatientIDs.RData")


   tbl.eqtl.geno.ldlr <- tbl.geno.ldlr[, -failures]
   dim(tbl.eqtl.geno.ldlr)   # 3299 1845

   stopifnot(all(colnames(tbl.eqtl.geno.ldlr) ==  tbl.id$vcf))
   colnames(tbl.eqtl.geno.ldlr) <- tbl.id$patient
   save(tbl.eqtl.geno.ldlr, file="tbl.eqtl.geno.ldlr-3299x1845-24may2021.RData")


      #------------------------------------------------------------------------
      # now use these maps to relabel the columns of (a copy of) tbl.geno.ldr
      #------------------------------------------------------------------------

   new.colnames <- colnames(tbl.geno.ldlr)

   new.colnames[match(names(map.rosmap), new.colnames)] <- as.character(map.rosmap)
   new.colnames[match(names(map.sinai), new.colnames)] <- as.character(map.sinai)
   new.colnames[match(names(map.mayo), new.colnames)] <- sprintf("mayo.%s", as.character(map.mayo))


   tbl.eqtl.geno.ldlr <- tbl.geno.ldlr
   colnames(tbl.eqtl.geno.ldlr) <- new.colnames
   length(grep("^R",    colnames(x))) # [1] 1171
   length(grep("^mayo", colnames(x))) # 349
   length(grep("^AMPAD_MSSM_", colnames(x))) # 345



   save(tbl.eqtl.geno.ldlr, file="tbl.eqtl.geno.ldlr.RData")

} # mapToPatientID
#----------------------------------------------------------------------------------------------------
test_vcfToPatientMapping <- function()
{
   dir <- "~/github/TrenaProjectAD/explore/WGS-Harmonization/metadata-all"
   file.exists(dir)
   tbl.rosmap <- read.table(file.path(dir, "ROSMAP_biospecimen_metadata.csv"), sep=",", as.is=TRUE, header=TRUE)
   tbl.sinai <- read.table(file.path(dir, "MSBB_biospecimen_metadata.csv"), sep=",", as.is=TRUE, header=TRUE)
   tbl.mayo <- read.table(file.path(dir, "MayoRNAseq_biospecimen_metadata.csv"), sep=",", as.is=TRUE, header=TRUE)

   tbl.id <- get(load("~/github/TrenaProjectAD/explore/ampad.eQTLS/tbl.vcfToPatientIDs.RData"))

   set.seed(17)
   for(i in seq_len(1000)){
     index <- sample(seq_len(nrow(tbl.id)), 1)
     cohort <- tbl.id$cohort[index]
     vcf.id <- tbl.id$vcf[index]
     patient.id <- tbl.id$patient[index]
     printf("%s: %s  %s", cohort, vcf.id, patient.id)
     if(cohort == "rosmap"){
        checkEquals(subset(tbl.rosmap, specimenID == vcf.id)$individualID, patient.id)
        printf("good rosmap: %s", vcf.id)
        } # if rosmap
     if(cohort == "mayo"){
        checkEquals(as.character(subset(tbl.mayo, specimenID == vcf.id)$individualID)[1], patient.id)
        printf("good mayo: %s", vcf.id)
        }
     if(cohort == "sinai"){
        checkEquals(subset(tbl.sinai, specimenID == vcf.id)$individualID, patient.id)
        printf("good mayo: %s", vcf.id)
        }
     } # for i

} # test_vcfToPatientMapping
#----------------------------------------------------------------------------------------------------
test_finalGenoTable <- function()
{
  tbl.0 <- get(load("ldlr-geno.RData"))
  tbl.1 <- get(load("tbl.eqtl.geno.ldlr.RData"))

     # tbl.1 got a loc column added at the end, remove it for these tests
  tbl.1 <- tbl.1[,-ncol(tbl.1)]

     # tbl.1 colnames got "fixed" somewhere today.  unfix them:
  colnames(tbl.1) <- colnames(tbl.0)
  tbl.2 <- get(load("tbl.eqtl.geno.ldlr-3299x1845-24may2021.RData"))
  tbl.id <- get(load("~/github/TrenaProjectAD/explore/ampad.eQTLS/tbl.vcfToPatientIDs.RData"))

    #------------------------------------------------------------
    # pick a patient from the mapped tbl.2, find the corresponding
    # column in tbl.0
    # check that all of the variant calls are identical
    #------------------------------------------------------------

  for(mapped.sample in colnames(tbl.2)){
    printf("mapped.sample: %s", mapped.sample)
    vcf.id <- subset(tbl.id, patient==mapped.sample)$vcf[1]
    if(!vcf.id %in% colnames(tbl.1)){
        printf("vcf.id missing: %s", vcf.id)
        next;
        }
    variant.calls.1 <- tbl.1[, vcf.id]
    variant.calls.2 <- tbl.2[, mapped.sample]
    checkEquals(length(variant.calls.1), length(variant.calls.2))
    checkEquals(variant.calls.1, variant.calls.2)
    } # for mapped.sample

} # test_finalGenoTable
#----------------------------------------------------------------------------------------------------
use.rsid.rowNames <- function()
{
   tbl.2 <- get(load("tbl.eqtl.geno.ldlr-3299x1845-24may2021.RData"))
   tbl.eqtls <- get(load("ldlr-eqtls.RData"))

   head(rownames(tbl.2))

   locs <- as.integer(sub("_.*$", "", sub("^.*:", "", rownames(tbl.2))))
   matches <- match(locs, tbl.eqtls$snpLocation)
   length(which(is.na(matches)))  # 0
   rsids <- tbl.eqtls$snpid[matches]
   stopifnot(length(rsids) == nrow(tbl.2))
   rownames(tbl.2) <- rsids

   tbl <- tbl.2
   save(tbl, file="tbl.eqtl.geno.ldlr-3299x1845-rsids-24may2021.RData")


} # use.rsid.rowNames
#---------------------------------------------------------------------------------------------------
