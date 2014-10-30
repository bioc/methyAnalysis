\name{filterBisulfiteVariant}
\alias{filterBisulfiteVariant}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
filtering the variant calls of Bisulfite converted sequencing data
}
\description{
a VRanges object output by NGS pipeline implemented in HTSeqGenie package
}
\usage{
filterBisulfiteVariant(seqVariant, coverage, CpGInfo, cleanVariant = TRUE, minCoverage = 1, convertTh.nonCpG = 0.9)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seqVariant}{
  a VRanges object output by NGS pipeline implemented in HTSeqGenie package
}
  \item{coverage}{
  the genome coverage (a RleList object) output by NGS pipeline
}
  \item{CpGInfo}{
  the pre-calculated CpG-site information (by identifyCpG function)
}
  \item{cleanVariant}{
  whether to filter those variants on non-CpG sites with full CT and GA conversion, or non CT and GA variations
}
  \item{minCoverage}{
  minimum coverage for the variants
}
  \item{convertTh.nonCpG}{
  convert rate threshold used by filtering variants on non-CpG sites (cleanVariant is TRUE)
}
}
\details{
Filtering the variant calls based on following criteria:
Only keeps CG and GA conversion variants.
The variants should have minimum coverage.
When cleanVariant parameter is TRUE, those fully converted non-CpG sites (convert rate higher than convertTh.nonCpG) will be removed.
}
\value{
a filtered VRanges object with two attributes: 'variantStats' and 'filterSettings'
}

\author{
Pan Du
}

\seealso{
\code{\link{estimateMethySeq}}
}
\examples{
}
\keyword{methods}

