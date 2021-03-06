library(VariantAnnotation)
library(RUnit)

file <- "NIA_JG_1898_samples_GRM_WGS_b37_JointAnalysis01_2017-12-08_7.recalibrated_variants.Broad_Rush.vcf.gz"
vcf.path <- file.path("/proj/price4/cory/AMP-AD_VCF/ROSMAP", file)
file.exists(vcf.path)
index.path <- file.path("/proj/price4/cory/AMP-AD_VCF/ROSMAP", sprintf("%s.tbi", file))
file.exists(index.path)

   # define a region for readVcf: just the PILRA gene

start <- 100371659
end   <- 100401098
1 + end - start   # 29k

region <- GRanges(seqnames="7", IRanges(start, end))
vcf <- readVcf(vcf.path, "hg38", region)
length(vcf) # 876

mtx.geno <- geno(vcf)[[1]]
dim(mtx.geno)

set.seed(17)
sample(colnames(mtx.geno), size=40)
  # [1] "SM-CJEHU"    "SM-CJGM8"    "SM-CJGHV"    "SM-CJEJA"    "SM-CTEM5"    "SM-CJK3F"   
  # [7] "SM-CJGIY"    "SM-CJK58"    "SM-CTEN2"    "SM-CJGHC"    "SM-CTDTM"    "SM-CJK5P"   
  #[13] "SM-CTDS3"    "SM-CJGIG"    "MAP50302680" "SM-CJIYH"    "SM-CTDSZ"    "SM-CTEMP"   
  #[19] "SM-CJEFZ"    "SM-CTEDW"    "SM-CTEF7"    "SM-CJEIS"    "SM-CJIXK"    "SM-CTEEV"   
  #[25] "SM-CJIX8"    "SM-CTEDG"    "SM-CTDSF"    "SM-CJEIJ"    "SM-CTEDO"    "SM-CTDU6"   
  #[31] "SM-CJIX1"    "MAP33332646" "SM-CJGMY"    "SM-CTDQP"    "SM-CJEJ6"    "SM-CTDS5"   
  #[37] "MAP50106992" "SM-CTEEI"    "SM-CTEEF"    "SM-CJJ25"   
  #     "ROS10524640"
ncol(mtx.geno)                             #  1196
length(grep("^SM-", colnames(mtx.geno)))   #  1151
length(grep("^MAP", colnames(mtx.geno)))   #    25
length(grep("^ROS", colnames(mtx.geno)))   #    20

    # same information available from scanVcfHeader

head(samples(scanVcfHeader(vcf.path)), n=50)
  #  [1] "MAP15387421" "MAP22868024" "MAP26637867" "MAP29629849" "MAP33332646" "MAP34726040"
  #  [7] "MAP46246604" "MAP46251007" "MAP50103679" "MAP50104134" "MAP50104846" "MAP50106442"
  # [13] "MAP50106992" "MAP50108462" "MAP50301099" "MAP50302680" "MAP50402431" "MAP50409406"
  # [19] "MAP61344957" "MAP85980779" "MAP87264456" "MAP89164957" "MAP93787649" "MAP95330358"
  # [25] "MAP95453354" "ROS10524640" "ROS11430815" "ROS11697592" "ROS15114174" "ROS15738428"
  # [31] "ROS20225051" "ROS20251553" "ROS20254452" "ROS20275399" "ROS20327084" "ROS20376029"
  # [37] "ROS20945666" "ROS20946257" "ROS20990085" "ROS20998065" "ROS21001807" "ROS21112011"
  # [43] "ROS21113864" "ROS21274866" "ROS79590778" "SM-CJEFQ"    "SM-CJEFR"    "SM-CJEFS"   
  # [49] "SM-CJEFU"    "SM-CJEFV"   

  # what does the genotype matrix look like?

mtx.geno[5:10, 5:10]
   #                MAP33332646 MAP34726040 MAP46246604 MAP46251007 MAP50103679 MAP50104134
   #7:100371765_C/. "0/0"       "0/0"       "0/0"       "0/0"       "0/0"       "0/0"      
   #7:100371785_T/C "0/1"       "0/1"       "0/1"       "0/1"       "0/1"       "0/0"      
   #7:100371884_C/T "0/0"       "0/0"       "0/0"       "0/0"       "0/0"       "0/0"      
   #7:100371946_C/T "0/0"       "0/0"       "0/0"       "0/0"       "0/0"       "0/0"      
   #7:100371956_C/A "0/0"       "0/0"       "0/0"       "0/0"       "0/0"       "0/0"      
   #7:100372120_C/T "0/0"       "0/0"       "0/0"       "0/0"       "0/0"       "0/0"      


   # save this small vcf for future easy use
   # 28M vs 21G for original

checkEquals(sort(samples(scanVcfHeader(vcf.path))), sort(colnames(mtx.geno)))  # TRUE

writeVcf(vcf, filename="pilra-region.vcf")  

samples.vcf <- sort(colnames(mtx.geno))
length(samples.vcf)  # 1196

tbl.biospecimen <- read.table("metadata-all/ROSMAP_biospecimen_metadata.csv",
                              sep=",", header=TRUE, as.is=TRUE)
dim(tbl.biospecimen)
head(tbl.biospecimen)
missing.ids <- setdiff(samples.vcf, tbl.biospecimen$specimenID)
length(missing.ids)  # 45
all(samples.vcf %in% tbl.biospecimen$specimenID)

length(tbl.biospecimen$individualID)          # 10944
length(unique(tbl.biospecimen$individualID))  # 2481

tbl.clinical <-
    read.table("metadata-all/ROSMAP_clinical.csv", sep=",", as.is=TRUE, header=TRUE, nrow=-1)
dim(tbl.clinical)  # 3584 18
length(unique(subset(tbl.clinical, Study=="ROS")$individualID)) # 1451
table(tbl.clinical$Study)

head(tbl.biospecimen)


