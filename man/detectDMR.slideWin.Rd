\name{detectDMR.slideWin}
\alias{detectDMR.slideWin}
\title{
Detect DMR (Differentially Methylated Region) using slide window smoothing
}
\description{
Detect DMR (Differentially Methylated Region) using slide window smoothing
}
\usage{
detectDMR.slideWin(methyGenoSet, sampleType, winSize = 250, testMethod = c("ttest", "wilcox"), 
  p.adjust.method = "fdr", p.value.detection.th = 0.05, ...)
}
\arguments{
  \item{methyGenoSet}{
A GenoSet object includes the methylation data information
}
  \item{sampleType}{
A vector shows the sample type information
}
  \item{winSize}{
Slide window size (half window size, bp at each side of the probe)
}
  \item{testMethod}{
test methods
}
  \item{p.adjust.method}{
p.value FDR adjust method
}
  \item{p.value.detection.th}{
  the threshold of detection p.value used to determine the failed probes, which will be set as NAs.
}
  \item{\dots}{
other paramters
}
}
\details{
The function will check whether the data was previously smoothed. If not, slide-window smoothing will be performed first, and then followed by differential methylation tests 
}
\value{
A GRanges object with additional test information (difference, p.value, p.adjust, and etc.)
}

\author{
Pan Du
}

\seealso{
 \code{\link{identifySigDMR}}
}
\examples{

	data(exampleMethyGenoSet)
	
	sampleType <- pData(exampleMethyGenoSet)$SampleType
	
	## Do differential test
	allResult <- detectDMR.slideWin(exampleMethyGenoSet, sampleType=sampleType, testMethod='ttest')
	head(allResult)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
