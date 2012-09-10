\name{identifySigDMR}
\alias{identifySigDMR}
\title{
Identify significantly DMR (Differentially Methylated Region)
}
\description{
Identify significantly DMR (Differentially Methylated Region)
}
\usage{
identifySigDMR(detectResult, p.adjust.method = "fdr", pValueTh = 0.01, fdrTh = pValueTh, diffTh = 1, minProbeNum = 1, maxGap = 2000, minGap=100, oppositeMethylation = FALSE, topNum = NULL)
}
\arguments{
  \item{detectResult}{
A GRanges object or a data.frame with "CHROMOSOME" and "POSITION" columns
}
  \item{p.adjust.method}{
p.value FDR adjust method
}
  \item{pValueTh}{
The threshold of p.value
}
  \item{fdrTh}{
The threshold of FDR (adjusted p.value)
}
  \item{diffTh}{
The threshold of difference between two conditions
}
  \item{minProbeNum}{
Minimum number of probes in each DMR
}
  \item{maxGap}{
The maximum gap allowed between two nearby probes to be considered within a same DMR
}
  \item{minGap}{
If two nearby DMRs have a gap less than or equal to the minGap, they will be merged as a single DMR
}
  \item{oppositeMethylation}{
Whether require the averaged methylation levels in the DMR are in opposite direction
}
  \item{topNum}{
Whether only returns the top number of probes (ranked by p.value)
}
}
\details{
(with "status" column) 

We define a differentially methylated region (DMR) as a region, in which most measured CpG-sites are differentially methylated. To identify DMRs, the function first determines the differential methylation status of each probe, then merge them as a continuous region. The identifySigDMR function returns a list of two GRanges objects. The sigDMRInfo includes the identified DMRs, and the sigDataInfo includes all differentially methylated probe information.
}
\value{
A list of GRanges objects, sigDMRInfo and sigDataInfo. sigDMRInfo contains DMR information, while sigDataInfo includes the probe level information within the DMRs.
}

\author{
Pan Du
}

\seealso{
\code{\link{getContinuousRegion}}, \code{\link{annotateDMRInfo}}
}
\examples{
	data(exampleMethyGenoSet)
	
	sampleType <- pData(exampleMethyGenoSet)$SampleType
	
	## Do differential test
	allResult <- detectDMR.slideWin(exampleMethyGenoSet, sampleType=sampleType, testMethod='ttest')

	## Identify the DMR (Differentially Methylated Region) by setting proper parameters.
	## Here we just use default ones
	allDMRInfo <- identifySigDMR(allResult)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}

