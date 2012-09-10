\name{annotateGRanges}
\alias{annotateGRanges}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Annotate a GRanges object
}
\description{
Annotate a GRanges object based on a transcription database
}
\usage{
annotateGRanges(grange, annotationDatabase, CpGInfo = NULL, flankRange = 0, promoterRange = 2000, EntrezDB = "org.Hs.eg.db")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{grange}{
A GRanges object
}
  \item{annotationDatabase}{
Annotation database: a TranscriptDb package, TranscriptDb object or GRanges object.
}
  \item{CpGInfo}{
A Bed file or GRanges object, which keeps the CpG-island information
}
  \item{flankRange}{
The flank range to be added to the input GRanges object
}
  \item{promoterRange}{
Define the size of promoter range at the upstream of TSS. Users can also directly provide the GRanges object
}
  \item{EntrezDB}{
The Entrez database for mapping from Entrez ID to gene symbols
}
}
\details{
This function is to annotate a GRanges object to the gene promoters or bodies. The annotation information is attached as additional columns of the GRanges object values. 
}
\value{
Return an annotated GRanges object with the annotation information attached as additional columns.
}

\author{
Pan Du
}

\seealso{
See Also \code{\link{annotateDMRInfo}}
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
	sigDMRInfo <- allDMRInfo$sigDMRInfo
	class(sigDMRInfo)
	
	## Annotate significant DMR info
	if (require(TxDb.Hsapiens.UCSC.hg19.knownGene)) {
		sigDMRInfo.ann <- annotateDMRInfo(sigDMRInfo, 'TxDb.Hsapiens.UCSC.hg19.knownGene')
	}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}

