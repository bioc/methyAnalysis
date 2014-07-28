\name{annotateDMRInfo}
\alias{annotateDMRInfo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Annotate the DMR (Differentially Methylated Region)
}
\description{
Annotate the DMR (Differentially Methylated Region) information
}
\usage{
annotateDMRInfo(DMRInfo, annotationDatabase, CpGInfo = NULL, flankRange = 500, promoterRange = 2000, EntrezDB = "org.Hs.eg.db", as.GRanges = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DMRInfo}{
A GRanges object or a list of GRanges objects (sigDMRInfo and sigDataInfo), which is the return of \code{\link{identifySigDMR}}
}
  \item{annotationDatabase}{
Annotation database: a TxDb package, TxDb object or GRanges object.
}
  \item{CpGInfo}{
A Bed file or GRanges object, which keeps the CpG-island information
}
  \item{flankRange}{
The flank range to be added to the input GRanges object
}
  \item{promoterRange}{
Define the size of promoter range at the upstream of TSS. User can also directly provide the GRanges object
}
  \item{EntrezDB}{
The Entrez database for mapping from Entrez ID to gene symbols
}
  \item{as.GRanges}{
Whether return a GRanges object or a data.frame
}
}
\details{
This function is to annotate the DMRs to the gene promoters or bodies. The annotation information is attached as additional columns of the GRanges object values. 
}
\value{
Return a GRanges object or list of GRanges when the as.GRanges is TRUE. Or else it returns a data.frame or a list of data.frame objects (sigDMRInfo and sigDataInfo). The annotation information is attached as additional columns of the GRanges object values or the data.frame.
}

\author{
Pan Du
}

\seealso{
See Also \code{\link{annotateGRanges}}
}
\examples{
	data(exampleMethyGenoSet)
	## get sample type information
	sampleType <- pData(exampleMethyGenoSet)$SampleType
	
	## Do differential test
	allResult <- detectDMR.slideWin(exampleMethyGenoSet, sampleType=sampleType, testMethod='ttest')

	## Identify the DMR (Differentially Methylated Region) by setting proper parameters.
	## Here we just use default ones
	allDMRInfo <- identifySigDMR(allResult)

	## Annotate significant DMR info
	if (require(TxDb.Hsapiens.UCSC.hg19.knownGene)) {
		DMRInfo.ann <- annotateDMRInfo(allDMRInfo, 'TxDb.Hsapiens.UCSC.hg19.knownGene')
	}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
