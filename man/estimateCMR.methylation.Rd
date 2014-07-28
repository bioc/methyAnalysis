\name{estimateCMR.methylation}
\alias{estimateCMR.methylation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate the averaged methylation levels within a chromosome region or transcript promoter
}
\description{
Estimate the averaged methylation levels within a chromosome region defined as a GRanges object or transcript promoter 

}
\usage{
estimateCMR.methylation(cmr, methyGenoSet, tx2probe.corList = NULL, estimateFun = mean, probeAnnotation = NULL, selectGeneElement = c("exon1", "promoter"), mc.cores = min(12, detectCores()))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cmr}{
A GRanges object or transcript ID
}
  \item{methyGenoSet}{
A MethyGenoSet object keeps the DNA methylation data
}
  \item{estimateFun}{
The function used to estimate the methylation levels within the chromosome region
}
  \item{probeAnnotation}{
Pre-calculated probe annotation (a GRanges object)
}
  \item{selectGeneElement}{
Gene elements used to calculate the transcript promoter methylation levels if cmr GRanges object is not provided.
}
  \item{mc.cores}{
Number of cores used to calculate in parallel
}
}
\details{
}
\value{
A numeric matrix (row: cmr, column: samples)
}

\author{
Pan DU
}

\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
