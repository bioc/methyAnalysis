\name{createTranscriptTrack}
\alias{createTranscriptTrack}
\title{
Create a transcript annotation track
}
\description{
Create a transcript track, which is a GeneRegionTrack object
}
\usage{
createTranscriptTrack(gene, genomicFeature = "TxDb.Hsapiens.UCSC.hg19.knownGene", lib = "org.Hs.eg.db", 
	genome = "hg19", extendRange = c(2000, 2000), includeOtherGene=FALSE, includeGeneBody = TRUE, 
	thinBox_utrOnly = FALSE, background.title = "gray", fill = "#8282d2", ...)
}
\arguments{
  \item{gene}{
An Entrez gene ID or a GRanges object with length equals one
}
  \item{genomicFeature}{
a TranscriptDb library, TranscriptDb object, or Mart object
}
  \item{lib}{
Entrez annotation library
}
  \item{genome}{
The version of genome
}
  \item{extendRange}{
extended range on each side of the  gene
}
  \item{includeOtherGene}{
whether to include other genes in the same chromosome ranges, only useful when "gene" is a gene ID.
}
  \item{includeGeneBody}{
whether to include the whole gene body or not
}
  \item{thinBox_utrOnly}{
whether to only show UTRs as thin boxs in the plot
}
  \item{background.title}{
the background color of the title
}
  \item{fill}{
fill color for transcript track
}
  \item{\dots}{
other parameters
}
}
\details{
This function is to create a GeneRegionTrack object for visualization using Gviz package.
}
\value{
a GeneRegionTrack object
}
\author{
Pan Du
}


\seealso{
\code{\link[Gviz]{plotTracks}}, \code{\link{plotTracksWithDataTrackInfo}}, \code{\link{heatmapByChromosome}}, \code{\link{plotMethylationHeatmapByGene}}
}
\examples{
if (require(TxDb.Hsapiens.UCSC.hg19.knownGene) && require(Gviz)) {
	rangeTrack <- createTranscriptTrack('7157', genomicFeature = "TxDb.Hsapiens.UCSC.hg19.knownGene")
	# plotTracks(rangeTrack)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
