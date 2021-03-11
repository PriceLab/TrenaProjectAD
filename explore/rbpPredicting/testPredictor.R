source("Predictor.R")
library(RUnit)
if(!exists("p"))
    p <- Predictor$new("GATA2")
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
    test_ctor()

} # runTests
#----------------------------------------------------------------------------------------------------
test_ctor <- function()
{
    message(sprintf("--- test_ctor"))

    checkEquals(p$getTargetGene(), "GATA2")

} # test_ctor
#----------------------------------------------------------------------------------------------------
test_getRBP <- function()
{
    message(sprintf("--- test_getRBP"))
    tbl.rbp <-  p$getRBP()
    p$vizRBP(min.score=500)

} # test_getRBP
#----------------------------------------------------------------------------------------------------
test_getExpressionMatrix <- function()
{
    message(sprintf("--- test_getExpressionMatrix"))
    mtx.rna <- p$getExpressionMatrix()
    checkEquals(dim(mtx.rna), c(27171, 28))

} # test_getExpressionMatrix
#----------------------------------------------------------------------------------------------------
test_discoverTfCandidates <- function()
{
    message(sprintf("--- test_setDiscoverCandidates"))
    p$discoverTfCandidates(p.value=6e-7, elite.gh.only=TRUE, gh.min.score=500)
    tbl.gh   <- p$getGH()
    dim(tbl.gh)
    tbl.fimo <- p$getFimo()
    dim(tbl.fimo)
    checkTrue(length(unique(tbl.fimo$tf)) > 45)
    expected <- c("CTCF","CTCFL","E2F1","KLF5","SP3","ZNF384","E2F6","EGR1","ZNF263","SP2","SP1")
    checkTrue(all(expected %in% tbl.fimo$tf))

} # test_getExpressionMatrix
#----------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()
