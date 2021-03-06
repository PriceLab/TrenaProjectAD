
file <- "NIA_JG_1898_samples_GRM_WGS_b37_JointAnalysis01_2017-12-08_7.recalibrated_variants.Broad_Rush.vcf.gz"
vcf.path <- file.path("/proj/price4/cory/AMP-AD_VCF/ROSMAP", file)
file.exists(vcf.path)

samples.vcf <- samples(scanVcfHeader(vcf.path))
length(samples.vcf)   # 1196

tbl.biospecimen <- read.table("metadata-all/ROSMAP_biospecimen_metadata.csv", sep=",", header=TRUE, as.is=TRUE)
tbl.clinical <-  read.table("metadata-all/ROSMAP_clinical.csv", sep=",", as.is=TRUE, header=TRUE, nrow=-1)

mtx.rna <- get(load("/users/pshannon/github/TrenaProjectAD/inst/extdata/expression/rosmap.14235x632.RData"))
samples.rna <- colnames(mtx.rna)
length(samples.rna)  # 632

samples.vcf.rna <- c(samples.vcf, samples.rna)
length(samples.vcf.rna) # 1828
length(unique(samples.vcf.rna)) # 1828
length(intersect(samples.vcf.rna, tbl.biospecimen$specimenID)) # 1783
length(setdiff(samples.vcf.rna, tbl.biospecimen$specimenID))   # 25 start with MAP, 20 start with ROS, all from samples.vcf

samples.all <- intersect(samples.vcf.rna, tbl.biospecimen$specimenID) # 1783

coi <- c("individualID", "specimenID", "organ", "tissue", "cellType")
tbl.xref <- subset(tbl.biospecimen, specimenID %in% samples.all)[, coi]
tbl.xref <- unique(tbl.xref)
rownames(tbl.xref) <- NULL
dim(tbl.xref)  # 1785 5
length(unique(tbl.xref$individualID))  # 1224

length(intersect(tbl.xref$specimenID, samples.rna))
length(intersect(tbl.xref$specimenID, samples.vcf))


ids.vcf <- subset(tbl.xref, specimenID %in% samples.vcf)$individualID
ids.rna <- subset(tbl.xref, specimenID %in% samples.rna)$individualID

ids.both <- intersect(ids.vcf, ids.rna)  
length(ids.both)  # 553
length(unique(ids.both)) # 553

dim(tbl.xref)  # 1785 5
tbl.xref <- subset(tbl.xref, individualID %in% ids.both)
dim(tbl.xref)  # 1111 5
assay <- rep("rna", nrow(tbl.xref))
assay[grep("^SM", tbl.xref$specimenID)] <- "vcf"
tbl.xref$assay <- assay


    # five individuals have two vcf entries.  keep them for now
    # "R3257830" "R4990009" "R7483736" "R7874995" "R8444624"
triples <- names(which(table(tbl.xref$individualID) == 3))
tbl.triples <- subset(tbl.xref, individualID %in% triples)
tbl.triples <- tbl.triples[order(tbl.triples$individualID),]

tbl.xref2 <- merge(tbl.xref, tbl.clinical, by="individualID", all.x=TRUE)
dim(tbl.xref)    # 1111 6
dim(tbl.xref2)   # 1111 23

wdth(20)
lapply(tbl.xref2, class)

    # three age columns, all character, all use "90+".
    # change those to "95", then cast to numeric, and round to 1 decimal digit

grep("age_", colnames(tbl.xref2), value=TRUE)
tbl.xref2$age_at_visit_max[tbl.xref2$age_at_visit_max == "90+"] <- "95"
tbl.xref2$age_at_visit_max <- as.numeric(tbl.xref2$age_at_visit_max)
tbl.xref2$age_at_visit_max <- round(tbl.xref2$age_at_visit_max, digits=1)
fivenum(tbl.xref2$age_at_visit_max) # [1] 67.0 83.8 88.2 95.0 95.0

tbl.xref2$age_first_ad_dx[tbl.xref2$age_first_ad_dx == "90+"] <- "95"
tbl.xref2$age_first_ad_dx <- as.numeric(tbl.xref2$age_first_ad_dx)
tbl.xref2$age_first_ad_dx <- round(tbl.xref2$age_first_ad_dx, digits=1)
fivenum(tbl.xref2$age_first_ad_dx) # 71.4 83.4 87.7 95.0 95.0

tbl.xref2$age_death[tbl.xref2$age_death == "90+"] <- "95"
tbl.xref2$age_death <- as.numeric(tbl.xref2$age_death)
tbl.xref2$age_death <- round(tbl.xref2$age_death, digits=1)
fivenum(tbl.xref2$age_death) # 67.4 84.7 89.0 95.0 95.0

tbl.cov <- tbl.xref2
colnames(tbl.cov)

     # "individualID" "specimenID" "organ" "tissue" "cellType" "assay"
     # "projid" "Study" "msex" "educ" "race" "spanish" "apoe_genotype"
     # "age_at_visit_max" "age_first_ad_dx" "age_death"
     # "cts_mmse30_first_ad_dx" "cts_mmse30_lv" "pmi" "braaksc"
     # "ceradsc" "cogdx" "dcfdx_lv"

     #  mmse30: 30 item questionare for dementia severity
     #      lv: last valid
     #     pmi: post-mortem interval in hours (between death and tissue preservation)
     # braaksc: semi-quantitative measure of neurofibrillary tangles.
     # ceradsc: assessment of neuritic plaques:
     #          value  coding
     #          -----  ------
     #              1 Definite
     #              2 Probable
     #              3 Possible
     #              4 No AD
     #              9 Missing
     #   cogdx: final clinical consensus diagnosis, blinded to all postmortem data
     #          value coding
     #          1     NCI, No cognitive impairment (No impaired domains)
     #          2     MCI, Mild cognitive impairment (One impaired domain) and NO other cause of CI
     #          3     MCI, Mild cognitive impairment (One impaired domain) AND another cause of CI
     #          4     AD, Alzheimer's disease and NO other cause of CI (NINCDS PROB AD)
     #          5     AD, Alzheimer's disease AND another cause of CI (NINCDS POSS AD)
     #          6     Other dementia. Other primary cause of dementia
     # dcfdx_lv: age where first AD dx was given

     # some modest renaming
cop <- c("patientID", "specimen", "organ", "tissue", "cellType", "assay", "project", "study", "gender",
         "educ", "race", "spanish", "apoe_genotype", "age_at_visit_max", "age_first_ad_dx", "age_death",
         "cts_mmse30_first_ad_dx", "cts_mmse30_lv", "pmi", "braaksc", "ceradsc", "cogdx", "dcfdx_lv")

colnames(tbl.cov) <- cop
save(tbl.cov, file="tbl.covariates.rna.vcf.RData")
fivenum(tbl.cov$apoe_genotype)
