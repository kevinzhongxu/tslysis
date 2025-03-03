% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lysis_rate.R
\name{lysis_rate}
\alias{lysis_rate}
\title{A lysis_rate Function}
\usage{
lysis_rate(otu_table, metadata)
}
\arguments{
\item{otu_table}{'Path/To/Your/OTU_TABLE_FILE'; For example: otu_table="/home/kevin/Desktop/data/10reads-plus/otu_table_qiime2_otu_16S.txt".}

\item{metadata}{'Path/To/Your/METADATA_FILE'; For example: metadata="/home/kevin/Desktop/data/10reads-plus/metadata.txt".}
}
\description{
This function allows calculating both the taxon-specific lysis and lysis rate of microbe, based on the input of a tabular-separated otu table (.txt) and metadata (.txt)
The otu_table contains data from three rDNA/rRNA fractions (Cellular rDNA, Cellular rRNA and Extracellular rRNA) of each water sample.
The metadata contains contains at least columns "sample.id", "Abundance_Total_Cellular.rRNA", "Abundance_Total_Extracellular.rRNA" and "turnover_rate".
}
\examples{
##### If you want to estimate both the taxon-specific lysis and lysis rate for prokaryotes, then run
lysis_rate(otu_table="Path/To/Your/OTU_TABLE_FILE_16S.txt", metadata="Path/To/Your/METADATA.txt")











}
\keyword{lysis_rate}
