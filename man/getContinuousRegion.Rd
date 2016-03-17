\name{getContinuousRegion}
\alias{getContinuousRegion}
\title{
Get continuous chromosome region by merging nearby or overlapping regions 
}
\description{
Get continuous chromosome region by merging nearby or overlapping regions 
}
\usage{
getContinuousRegion(detectResult, scoreColumns = NULL, scoreFuns = c(mean=mean), maxGap = 2000, minGap = 100)
}
\arguments{
  \item{detectResult}{
A GRanges object (with "status" column) or a data.frame with "CHROMOSOME", "POSITION" and "status" columns
}
  \item{scoreColumns}{
The numeric score columns to be summarized in DMR
}
  \item{scoreFuns}{
A named vector of summarizing functions. The vector names will be used in the output columns
}
  \item{maxGap}{
The maximum gap allowed between two nearby probes to be considered within a same DMR
}
  \item{minGap}{
If two nearby DMRs have a gap less than or equal to the minGap, they will be merged as a single DMR
}
}
\details{
The "status" column in the "detectRsult" parameter is required, which is a logical vector indicating the interested probes.
}
\value{
A GRanges object of DMR
}
\author{
Pan Du
}

\seealso{
\code{\link{identifySigDMR}}
}
\examples{

	data(exampleMethyGenoSet)
	## get sample type information
	sampleType <- colData(exampleMethyGenoSet)$SampleType
	
	## Do differential test
	allResult <- detectDMR.slideWin(exampleMethyGenoSet, sampleType=sampleType, testMethod='ttest')

	## Identify the DMR (Differentially Methylated Region) by setting proper parameters.
	## Here we simply using fdr.adjusted p.value cutoff 0.05 to define DMR
	## "status" column is required for getContinuousRegion function.
	values(allResult)$status <- values(allResult)$p.adjust < 0.05
	dmrInfo <- getContinuousRegion(allResult)
	
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{method}

