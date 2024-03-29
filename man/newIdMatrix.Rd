% Generated by roxygen2 (4.0.2): do not edit by hand
\name{newIdMatrix}
\alias{newIdMatrix}
\title{newIdMatrix}
\usage{
newIdMatrix(x, convertIdTable, genesKept = c("var", "foldchange", "abs",
  "first", "random"))
}
\arguments{
\item{x}{the expression data matrix.}

\item{convertIdTable}{A vector. The names should be the source IDs, and the values should be the target IDs.}

\item{genesKept}{The method to select target gene in more than one targets. "var"/"foldchange"/"abs" means selecting the gene with largest variation/fold change/absolute value. "first" means selecting the first target and "random" means randomly selection.}
}
\description{
A function to convert ID.
}
\details{
A function to convert ID.
}
\examples{
convertIdTable<-paste("New",c(1,2,2,2,1,3,4,4,5,5))
names(convertIdTable)<-paste("Old",1:length(convertIdTable))
temp<-matrix(rnorm(20),ncol=2)
row.names(temp)<-names(convertIdTable)
colnames(temp)<-c("Exp1","Exp2")
newIdMatrix(temp,genesKept="foldchange",convertIdTable)
}

