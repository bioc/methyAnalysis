\name{estimateMethySeq}
\alias{estimateMethySeq}
\title{
Estimate the methylation level (Beta-value) of Methyl-Seq data 
}
\description{
Estimate the methylation level (Beta-value) of Methyl-Seq data, which is a VRanges object output by NGS pipeline 
}
\usage{
estimateMethySeq(seqVariant, coverage, CpGInfo = NULL, mergeStrand = TRUE, cleanVariant = TRUE, minCoverage = 10)
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
  the precalculated CpG-site information (by identifyCpG function)
}
  \item{mergeStrand}{
  whether to merge the AT and GA conversion on opposite strands
}
  \item{cleanVariant}{
  whether to filter those non-CpG with full CT and GA conversion, or non CT and GA variations
}
  \item{minCoverage}{
  minimum coverage for the variants
}
}
\details{

}
\value{
a GRanges object with the Beta column shows the methylation levels
}

\author{
Pan Du
}

\seealso{

}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
