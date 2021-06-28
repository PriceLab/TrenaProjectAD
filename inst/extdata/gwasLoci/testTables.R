library(RUnit)

tbl.38 <- get(load("~/github/TrenaProjectAD/inst/extdata/gwasLoci/tbl.posthuma-38-loci-curated.RData"))
tbl.assoc <- get(load("~/github/TrenaProjectAD/inst/extdata/gwasLoci/tbl.posthuma-38-geneAssociations-curated-3828x12.RData"))

loci.genes <- sort(tbl.38$geneSymbol)
assoc.genes <- sort(unique(tbl.assoc$gene))

checkEquals(length(loci.genes), 38)
length(assoc.genes)
setdiff(loci.genes, assoc.genes)   # "APOE"     "HLA-DRB1"

checkEquals(sort(unique(tbl.assoc$gene)), sort(loci.genes))
colnames(tbl.assoc)
checkEquals(sort(unique(tbl.assoc$gene)), sort(loci.genes))

genes <- tbl.38$geneSymbol
for(g in genes){
  associated.loc <- mean(subset(tbl.assoc, gene==g)$hg38, na.rm=TRUE) # 86044126
  gene.loc <-subset(tbl.38, geneSymbol==g)$hg38    # 86089237
  delta <- abs(associated.loc-gene.loc)/gene.loc
  printf("--- %s: %f", g, delta)
  }


