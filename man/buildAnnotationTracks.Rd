\name{buildAnnotationTracks}
\alias{buildAnnotationTracks}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Build annotation tracks for visualizing using Gviz package
}
\description{
Build annotation tracks for visualizing using Gviz package
}
\usage{
buildAnnotationTracks(gene, extendRange = c(2000, 2000), includeGeneBody = TRUE, cytobandInfo = NULL, CpGInfo = NULL, genomeAxis = TRUE, lib = "org.Hs.eg.db", genome = "hg19", genomicFeature = "TxDb.Hsapiens.UCSC.hg19.knownGene", ...)
}
\arguments{
  \item{gene}{
An Entrez gene id or a GRanges object with length equals one
}
  \item{extendRange}{
extended range on each side of the  gene
}
  \item{includeGeneBody}{
whether to include genebody of the provided gene
}
  \item{cytobandInfo}{
cytoband information. Set NA to suppress it.
}
  \item{CpGInfo}{
CpG-island information, GRanges or bed file are supported
}
  \item{genomeAxis}{
whether to add genome axis or not
}
  \item{lib}{
gene annotation library
}
  \item{genome}{
genome version
}
  \item{genomicFeature}{
genomic features: "TranscriptDb" library or object, "Mart" object
}
  \item{\dots}{
other parameters used by createTranscriptTrack function
}
}
\details{
This function aims to build annotation tracks to be visualized using Gviz package. If the cytobandInfo and CpGInfo are NULL and internet connection is available, it will download information directly from UCSC website. Set them as NAs if you want suppress this default behavior.
}
\value{
A list of different annotation Tracks
}

\author{
Pan Du
}

\seealso{
\code{\link{help}}
}
\examples{
if (require(TxDb.Hsapiens.UCSC.hg19.knownGene) && require(Gviz)) {
	annotationTracks <- buildAnnotationTracks('1826', includeGeneBody = FALSE, genomicFeature = "TxDb.Hsapiens.UCSC.hg19.knownGene")
	plotTracks(annotationTracks)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
