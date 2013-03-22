\name{smoothMethyData}
\alias{smoothMethyData}
\title{
Smooth the methylation data
}
\description{
Smooth the methylation data by a sliding window with fixed width in bp unit
}
\usage{
smoothMethyData(methyData, winSize = 250, lib = "IlluminaHumanMethylation450k.db", ...)
}
\arguments{
  \item{methyData}{
A GenoSet or MethyLumiM object
}
  \item{winSize}{
Half sliding window size in bp unit at each side of the probe
}
  \item{lib}{
Methylation annotation library
}
  \item{\dots}{
other parameters
}
}
\details{
The function basically averages the probes within a local window (within winSize bp at each side of the probe).
}
\value{
An object with the methylation values smoothed
}
\author{
Pan Du
}

\seealso{
\code{\link{detectDMR.slideWin}}
}
\examples{
data(exampleMethyGenoSet)
smoothData <- smoothMethyData(exampleMethyGenoSet)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
