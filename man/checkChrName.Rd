\name{checkChrName}
\alias{checkChrName}
\title{
check chromosome names
}
\description{
Check chromosome names and make sure chromosome names start with "chr" (or not if addChr is FALSE)
}
\usage{
checkChrName(grange, addChr = TRUE)
}
\arguments{
  \item{grange}{
a GRanges, RangedData object, character or named vector
}
  \item{addChr}{
Whether to add "chr" in front of chromosome names
}
}
\details{
Because some annotation database names the chromosomes without "chr" prefix, while many others do, it causes problems when both types of data exist in the analysis. This function aims to resolve such issues by checking chromosome names and make sure chromosome names start with "chr" (or not if addChr is FALSE).
}
\value{
return the same type of object with chromosome names checked.
}
\author{
Pan Du
}

\examples{

data(exampleMethyGenoSet)
seqlevels(rowRanges(exampleMethyGenoSet))

tt <- checkChrName(exampleMethyGenoSet, addChr = TRUE)
seqlevels(rowRanges(tt))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
