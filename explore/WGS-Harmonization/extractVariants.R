library(VariantAnnotation)
# chr7:99,371,951-100,633,313
region <- GRanges(seqnames="7", IRanges(start=99371951, end=100633313))
#file <- file.path("/proj/price4/cory/AMP-AD_VCF/ROSMAP",
                                        #                  "NIA_JG_1898_samples_GRM_WGS_b37_JointAnalysis01_2017-12-08_7.recalibrated_variants.Broad_Rush.vcf.gz")
vcf.file <- file.path("~/github/TrenaProjectAD/inst/extdata/variants",
                  "NIA_JG_1898_samples_GRM_WGS_b37_JointAnalysis01_2017-12-08_7.recalibrated_variants.annotated.vcf")
file.exists(vcf.file)

tbi.filename <- indexTabix(vcf.file, format="vcf")  #  creates filename with .tbi suffix
file.exists(tbi.filename)
vcf <- readVcf(file, "hg38", region)
length(vcf) # 28966
writeVcf(vcf, filename="pilra-region.vcf")

tbl.info <- as.data.frame(info(vcf));#  [1] 5162   24
dim(tbl.info)  # 28966 24
wdth(100)
colnames(tbl.info)
mtx.geno <- geno(vcf)$GT
dim(mtx.geno)  #28966 1196
head(colnames(mtx.geno))

[1] "AC"                  "AF"                  "AN"                  "BaseQRankSum"
[5] "DP"                  "DS"                  "END"                 "ExcessHet"
[9] "FS"                  "HaplotypeScore"      "InbreedingCoeff"     "MLEAC"
[13] "MLEAF"               "MQ"                  "MQ0"                 "MQRankSum"
[17] "NEGATIVE_TRAIN_SITE" "POSITIVE_TRAIN_SITE" "QD"                  "ReadPosRankSum"
[21] "SOR"                 "VQSLOD"              "VariantType"         "culprit"

--- val answers this identical question in 2016:    https://support.bioconductor.org/p/79389/
ScanVcfParam(samples=c("samp1", "samp2", "samp45"))
head(samples(scanVcfHeader(f))) # [1] "1000" "1005" "1010" "1015" "1019" "1027"

