\name{getCoverage}
\alias{getCoverage}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
get the coverage based on a given GRanges object
}
\description{
calculate the average coverage on each element of a given GRanges object
}
\usage{
getCoverage(grange, coverage, startOnly = FALSE, as.GRanges = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{grange}{
  GRanges object specify the genome locations to get coverage
}
  \item{coverage}{
  the genome coverage (a RleList object) output by NGS pipeline
}
  \item{startOnly}{
  whether to calculate the coverage based on the start location or the average of the entire GRanges element
}
  \item{as.GRanges}{
  whether return a GRanges object or a vector of coverage
}
}
\value{
a vector of coverage or a GRanges object (as.GRanges is TRUE)
}
\author{
Pan Du
}

\examples{
}
\keyword{methods}
