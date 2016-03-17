\name{export.DMRInfo}
\alias{export.DMRInfo}
\title{
Output the DMR (Differentially Methylated Region) data information
}
\description{
Output the DMR (Differentially Methylated Region) data information
}
\usage{
export.DMRInfo(DMRInfo.ann, methyData = NULL, savePrefix = "")
}
\arguments{
  \item{DMRInfo.ann}{
The annotated DMR information outputted by annotateDMRInfo.
}
  \item{methyData}{
Methylation data information in MethyGenoSet or MethyLumiM class
}
  \item{savePrefix}{
The prefix added to the output file names.
}
}
\details{
This function basically save the annotated DMR information as text .csv files. 
}
\value{
results files.
}
\author{
Pan Du
}

\seealso{
 \code{\link{annotateDMRInfo}}
}
\examples{
	data(exampleMethyGenoSet)
	
	sampleType <- colData(exampleMethyGenoSet)$SampleType
	
	## Do differential test
	allResult <- detectDMR.slideWin(exampleMethyGenoSet, sampleType=sampleType, testMethod='ttest')

	## Identify the DMR (Differentially Methylated Region) by setting proper parameters.
	## Here we just use default ones
	allDMRInfo <- identifySigDMR(allResult)

	## Annotate significant DMR info
	if (require(TxDb.Hsapiens.UCSC.hg19.knownGene)) {
		DMRInfo.ann <- annotateDMRInfo(allDMRInfo, 'TxDb.Hsapiens.UCSC.hg19.knownGene')
		export.DMRInfo(DMRInfo.ann, savePrefix='testExample')
	}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
