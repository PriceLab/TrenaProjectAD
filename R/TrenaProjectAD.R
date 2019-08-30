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
                               contains="TrenaProjectHG38")

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
   igap.ad.genes <- sort(c("CR1", "BIN1", "CD2AP", "EPHA1", "CLU", "MS4A6A", "PICALM",
                           "ABCA7", "CD33", "HLA-DRB5", "HLA-DRB1", "PTK2B", "SORL1",
                           "SLC24A4", "RIN3", "DSG2", "INPP5D", "MEF2C", "NME8", "ZCWPW1",
                           "CELF1", "FERMT2", "CASS4", "APOE", "TOMM40", "TREM2"))

   footprintDatabaseHost <- "khaleesi.systemsbiology.net"
   footprintDatabasePort <- 5432

   footprintDatabaseNames <- c("brain_wellington_16",
                               "brain_wellington_20",
                               "brain_hint_16",
                               "brain_hint_20")

   dataDirectory <- system.file(package="TrenaProjectAD", "extdata")

   .TrenaProjectAD(TrenaProjectHG38(projectName="TrenaProjectAD",
                                    supportedGenes=igap.ad.genes,
                                    footprintDatabaseHost=footprintDatabaseHost,
                                    footprintDatabasePort=footprintDatabasePort,
                                    footprintDatabaseNames=footprintDatabaseNames,
                                    packageDataDirectory=dataDirectory,
                                    quiet=quiet
                                    ))

} # TrenaProjectAD, the constructor
#----------------------------------------------------------------------------------------------------
