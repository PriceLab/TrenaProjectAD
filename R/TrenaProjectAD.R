#----------------------------------------------------------------------------------------------------
#' @import methods
#' @import TrenaProject
#'
#' @title TrenaProjectAD-class
#'
#' @name TrenaProjectAD-class
#' @rdname TrenaProjectAD-class
#' @aliases TrenaProjectAD
#' @exportClass TrenaProjectAD
#'

.TrenaProjectAD <- setClass("TrenaProjectAD",
                               contains="TrenaProject")

#----------------------------------------------------------------------------------------------------
#' Define an object of class TrenaProjectPlacneta
#'
#' @description
#' Expression, variant and covariate data for the 21 IGAP loci
#'
#' @rdname TrenaProject-class
#'
#' @export
#'
#' @return An object of the TrenaProjectAD class
#'

TrenaProjectAD <- function(quiet=TRUE)

{
   genomeName <- "hg38"

   igap.ad.genes <- sort(c("CR1", "BIN1", "CD2AP", "EPHA1", "CLU", "MS4A6A", "PICALM",
                           "ABCA7", "CD33", "HLA-DRB5", "HLA-DRB1", "PTK2B", "SORL1",
                           "SLC24A4", "RIN3", "DSG2", "INPP5D", "MEF2C", "NME8", "ZCWPW1",
                           "CELF1", "FERMT2", "CASS4", "APOE", "TOMM40", "TREM2"))

   footprintDatabaseNames <- c("brain_wellington_16",
                               "brain_wellington_20",
                               "brain_hint_16",
                               "brain_hint_20")

     # very temporarily
   geneInfoTable.path <- system.file(package="TrenaProject", "extdata", "geneInfoTable.RData")

   expressionDirectory <- system.file(package="TrenaProjectAD", "extdata", "expression")
   variantsDirectory <- system.file(package="TrenaProjectAD", "extdata", "variants")
   footprintDatabaseHost <- "khaleesi.systemsbiology.net"

   covariatesFile <- system.file(package="TrenaProjectAD", "extdata", "covariates", "covariates.RData")
   stopifnot(file.exists(expressionDirectory))
   stopifnot(file.exists(variantsDirectory))
   stopifnot(file.exists(covariatesFile))

   .TrenaProjectAD(TrenaProject(supportedGenes=igap.ad.genes,
                                  genomeName=genomeName,
                                  geneInfoTable.path=geneInfoTable.path,
                                  footprintDatabaseHost=footprintDatabaseHost,
                                  footprintDatabaseNames=footprintDatabaseNames,
                                  expressionDirectory=expressionDirectory,
                                  variantsDirectory=variantsDirectory,
                                  covariatesFile=covariatesFile,
                                  quiet=quiet
                                  ))

} # TrenaProjectAD, the constructor
#----------------------------------------------------------------------------------------------------
