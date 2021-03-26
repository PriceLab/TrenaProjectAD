library(VariantAnnotation)
# --- rs7412
#   C>T
#    19:44908822 (GRCh38)
#    19:45412079 (GRCh37)
# rs429358
# type:SNVAlleles:T>C
#   19:44908684 (GRCh38)
#   19:45411941 (GRCh37)

regions <- GRanges(seqnames="19", IRanges(start=45411941, end=45412079))
vcf.file <- "NIA_JG_1898_samples_GRM_WGS_b37_JointAnalysis01_2017-12-08_19.recalibrated_variants.vcf.gz"
file.exists(vcf.file)
file.exists(sprintf("%s.tbi", vcf.file))
vcf <- readVcf(vcf.file, "hg19", regions)
length(vcf)
dim(geno(vcf)$GT)  # 5 1894

mtx.geno <- geno(vcf)$GT

pathAging.ids <- c("18202_TCX", "18213_TCX", "18214_TCX", "18221_TCX", "18225_TCX", "18232_TCX",
                   "18235_TCX","18237_TCX", "12308_TCX", "18198_TCX", "18199_TCX", "18203_TCX",
                   "18204_TCX", "18209_TCX","18222_TCX", "2517_TCX", "18212_TCX", "18230_TCX",
                   "18231_TCX", "18234_TCX", "18236_TCX", "18201_TCX", "18210_TCX", "18239_TCX",
                   "2210_TCX", "12306_TCX", "18218_TCX", "18223_TCX", "18224_TCX", "18233_TCX",
                   "18238_TCX", "18206_TCX", "18216_TCX", "18219_TCX", "18226_TCX", "18229_TCX",
                   "18200_TCX", "18211_TCX", "18217_TCX", "18220_TCX", "18228_TCX")

pathAging.ids <- sub("_TCX", "", pathAging.ids)
length(intersect(pathAging.ids, colnames(mtx.geno))) # [1] 41  all present

mtx <- mtx.geno[c(1, 5), pathAging.ids]

# rs429358	rs7412	Name
#        C	T	ε1
#        T	T	ε2
#        T	C	ε3
#        C	C	ε4
#
# Common name	Genoset	Magnitude	rs429358	rs7412	Comment
# Apo-ε1/ε1	gs267	6	           (C;C)	(T;T)	the rare missing allele
# Apo-ε1/ε2	gs271	2.5	           (C;T)	(T;T)
# Apo-ε1/ε3	gs270	2.6	           (C;T)	(C;T)	ambiguous ε2/ε4 or ε1/ε3
# Apo-ε2/ε4	gs270	2.6	           (C;T)	(C;T)	ambiguous ε2/ε4 or ε1/ε3
# Apo-ε1/ε4	gs272	2.5	           (C;C)	(C;T)
# Apo-ε2/ε2	gs268	4	           (T;T)	(T;T)	good; lowest risk
# Apo-ε2/ε3	gs269	2	           (T;T)	(C;T)
# Apo-ε3/ε3	gs246	2	           (T;T)	(C;C)	the most common
# Apo-ε3/ε4	gs141	3	           (C;T)	(C;C)
# Apo-ε4/ε4	gs216	6	           (C;C)	(C;C)	~11x increased Alzheimer's risk

tbl <- as.data.frame(t(mtx))
tbl <- tbl[order(as.integer(rownames(tbl))),]
rs429358.variants <- which(tbl[, 1] == "0/1")
rs429358.wt <- setdiff(seq_len(nrow(tbl)), rs429358.variants)
rs7412.variants <- which(tbl[, 2] == "0/1")
rs7412.wt <- setdiff(seq_len(nrow(tbl)), rs7412.variants)

tbl$rs429358 <- "T"
tbl[rs429358.variants, "rs429358"] <- "C"
tbl$rs7412   <- "C"
tbl[rs7412.variants, "rs7412"] <- "T"

tbl$code <- 33
tbl[intersect(rs429358.variants, rs7412.wt), "code"] <- 34          # CC
tbl[intersect(rs429358.variants, rs7412.variants), "code"] <- 24    # CT
tbl[intersect(rs429358.wt, rs7412.variants), "code"] <- 23          # TT

# gustavo notes:
#  Individual 18231 has genotypes 1/1 0/0, which means 4/4. Seems to be miscoded as 3/3 in the table.
#  This cohort has a surprisingly high frequency of 2 alleles...

tbl[match("18231", tbl$individualID), "code"] <- 44

metadata.file <- "../../WGS-Harmonization/metadata-all/MayoRNAseq_individual_metadata.csv"
file.exists(metadata.file)
tbl.md <- read.table(metadata.file, sep=",", header=TRUE, as.is=TRUE)

merge.indices <- match(rownames(tbl), tbl.md$individualID)
tbl.combined <- cbind(tbl, tbl.md[merge.indices,])
rownames(tbl.combined) <- NULL

colanmes(tbl.combined)

coi <- c("individualID", "19:45411941_T/C", "19:45412079_C/T", "rs429358", "rs7412", "code",
         "individualIdSource", "species", "sex", "race", "ethnicity", "yearsEducation",
         "ageDeath", "causeDeath", "mannerDeath", "apoeGenotype", "pmi", "pH", "brainWeight",
         "diagnosis", "diagnosisCriteria", "CERAD", "Braak", "thal")

tbl.combined <- tbl.combined[, coi]

write.table(tbl.combined, file="mayoPathAgingSamplesAnnotated.tsv", sep="\t", quote=FALSE,
            row.names=FALSE, col.names=TRUE)

