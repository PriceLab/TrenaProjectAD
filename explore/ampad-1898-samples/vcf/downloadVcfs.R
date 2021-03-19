library(synapser) # 0.9.77
synLogin("paul-shannon", password="tril0byt")
files <- list("chr19.vcf" = "syn11714133", # NIA_JG_1898_samples_GRM_WGS_b37_JointAnalysis01_2017-12-08_19.recalibrated_variants.vcf.gz
              "chr19.tbi" = "syn11714174") # NIA_JG_1898_samples_GRM_WGS_b37_JointAnalysis01_2017-12-08_19.recalibrated_variants.vcf.gz.tbi

f <- synGet(files[[1]], downloadLocation=".")
f <- synGet(files[[2]], downloadLocation="./")

wdth(80)
list.files()
