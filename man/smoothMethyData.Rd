\name{smoothMethyData}
\alias{smoothMethyData}
\title{
Smooth the methylation data
}
\description{
Smooth the methylation data by a sliding window with fixed width in bp unit
}
\usage{
smoothMethyData(methyData, winSize = 250, lib = "IlluminaHumanMethylation450k.db", bigMatrix=FALSE, 
		dir.bigMatrix='.', savePrefix.bigMatrix=NULL, ...)
}
\arguments{
  \item{methyData}{A GenoSet or MethyLumiM object}
  \item{winSize}{Half sliding window size in bp unit at each side of the probe}
  \item{lib}{Methylation annotation library}
  \item{bigMatrix}{whether to save the data as BigMatrix (designed for very large dataset)}
  \item{dir.bigMatrix}{the parent directory to save the BigMatrix data files}
  \item{savePrefix.bigMatrix}{the folder name prefix of the directory to save the BigMatrix data files. The fold name will be like this: paste(savePrefix.bigMatrix, '_bigmat', sep='') }
  \item{\dots}{other parameters}
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
