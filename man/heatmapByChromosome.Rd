\name{heatmapByChromosome}
\alias{heatmapByChromosome}
\title{
heatmap with chromosome location as x axis
}
\description{
heatmap with chromosome location as x axis and plot together with other gene annotation information
}
\usage{
heatmapByChromosome(genoSet, gene, annotationTracks = NULL, otherTrackList = NULL, phenoData = NULL,  
	phonoColorMap = NULL, extendRange = c(2000, 2000), includeGeneBody = TRUE, showFullModel = FALSE, sortSample = TRUE,  
	cytobandInfo = NULL, CpGInfo = NULL, genomeAxis = TRUE, dataTrackName = "Methylation Profile", lib = "org.Hs.eg.db", 
	genome = "hg19", genomicFeature = "TxDb.Hsapiens.UCSC.hg19.knownGene", gradient = c("blue", "white", "red"), 
	ncolor = 16, ylim = NULL, th = 0.99, main = "", selSample = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genoSet}{
a GenoSet object keeping the methylation data as the "exprs" numeric matrix in the AssayData
}
  \item{gene}{
a Entrez Gene ID, or a GRanges object to define the chromosome range of the plot.
}
  \item{annotationTracks}{
A annotation tracks list returned by buildAnnotationTracks
}
  \item{otherTrackList}{
A list of other tracks supported by plotTracks function
}
  \item{phenoData}{
a data matrix with the same number of rows or columns as the columns of genoSet.
}
  \item{phonoColorMap}{
the colormap for expression heatmap
}
  \item{extendRange}{
extended range on each side of the  gene
}
  \item{includeGeneBody}{
whether to include genebody of the provided gene
}
  \item{showFullModel}{
whether to show full gene model track when includeGeneBody = FALSE
}
  \item{sortSample}{
whether to sort samples or not
}
  \item{cytobandInfo}{
cytoband information
}
  \item{CpGInfo}{
CpG-island information, GRanges or bed file are supported
}
  \item{genomeAxis}{
whether to add genome axis or not
}
  \item{dataTrackName}{
the title of the data track
}
  \item{lib}{
the Entrez annotation library
}
  \item{genome}{
genome name
}
  \item{genomicFeature}{
genomic features: "TranscriptDb" library or object, "Mart" object
}
  \item{gradient}{
the gradient color used by data track heatmap
}
  \item{ncolor}{
the number of color levels
}
  \item{ylim}{
the range for plotting the data.
}
  \item{th}{
the quantile threshold used to remove outlier, which affects the plot color ranges.
}
  \item{main}{
the title of the plot
}
  \item{selSample}{
subset of samples for plotting. It is designed for BigMatrix, which have to extract the data at the last moment.
}
  \item{\dots}{
other parameters used by \code{\link{plotTracksWithDataTrackInfo}}
}
}
\details{
This function plots heatmap with chromosome location as x axis and together with other gene annotation information. It is adapted based on the \code{\link[Gviz]{plotTracks}} function in Gviz package. Users can also provide a GRanges object to specify a plot region.
}
\value{
returns the grid viewport layout information
}
\author{
Pan Du
}

\seealso{
\code{\link[Gviz]{plotTracks}}, \code{\link{plotTracksWithDataTrackInfo}}, \code{\link{plotMethylationHeatmapByGene}}
}
\examples{

	data(exampleMethyGenoSet)
	if (require(TxDb.Hsapiens.UCSC.hg19.knownGene) && require(Gviz)) {
		## define data track
		exampleMethyGenoSet <- checkChrName(exampleMethyGenoSet)

		## build annotation tracks
		selGene <- '1826'
		annotationTracks <- buildAnnotationTracks(selGene, includeGeneBody = FALSE, genomicFeature = "TxDb.Hsapiens.UCSC.hg19.knownGene", cytobandInfo=NA, CpGInfo=NA)
		heatmapByChromosome(exampleMethyGenoSet, selGene, annotationTracks = annotationTracks)
	}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
\keyword{hplot}
