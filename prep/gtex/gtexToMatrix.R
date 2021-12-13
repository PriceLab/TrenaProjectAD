library(EnsDb.Hsapiens.v79)
tbl <- read.table("Brain_Cortex.v8.normalized_expression.bed", sep="\t", as.is=TRUE, nrow=-1,
                  row.names=4)[, -(1:3)]
dim(tbl) # 24849   205

rownames(tbl) <- sub("\\..*$", "", rownames(tbl))

map <- mapIds(EnsDb.Hsapiens.v79, rownames(tbl), "SYMBOL", "GENEID")
  # Unable to map 152 of 24849 requested IDs.
tbl.map <- data.frame(ensg=names(map), symbol=as.character(map), stringsAsFactors=FALSE)
na.indices <- which(is.na(tbl.map$symbol))
length(na.indices)
tbl.map$symbol[na.indices] <- tbl.map$ensg[na.indices]

   # make sure we the two tables are congruent
stopifnot(nrow(tbl) == nrow(tbl.map))

dups <- which(duplicated(tbl.map$symbol))

length(dups) # 28
failures <- names(which(is.na(map)))
map[failures]

tbl.map <- tbl.map[-dups,]
tbl.2 <- tbl[-dups,]

dim(tbl.2)
dim(tbl.map)
rownames(tbl.2) <- tbl.map$symbol

mtx <- as.matrix(tbl.2)
new.colnames <- sprintf("s.%03d", seq_len(ncol(mtx)))
colnames(mtx) <- new.colnames
save(mtx, file="~/github/TrenaProjectAD/inst/extdata/expression/gtex.v8/Brain_Cortex.gtexV8.24821x205.RData")

