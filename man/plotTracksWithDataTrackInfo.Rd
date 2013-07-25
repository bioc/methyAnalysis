\name{plotTracksWithDataTrackInfo}
\alias{plotTracksWithDataTrackInfo}
\title{
plot Tracks with additional DataTrack information added to the left of the plot
}
\description{
plot Tracks with additional DataTrack information added to the left of the plot
}
\usage{
plotTracksWithDataTrackInfo(trackList, labels = NULL, grange2show = NULL, dataTrackName = NULL, dataInfo = NULL, dataColorMap = NULL,
	dataInfoRange = NULL, dataBackground = gray(0.9), minHeatmapColumnWidth = 2, labelWidth = 0.1, 
	gradient = c("blue", "white", "red"), ncolor = 16, main = "", newPlot = FALSE, sizes = NULL, ...)
}
\arguments{
  \item{trackList}{
a list of tracks supported by \code{\link[Gviz]{plotTracks}} function
}
  \item{labels}{
the sample labels. By default, it will use the rownames of dataTrack. It can also be a list if there are multiple dataTracks. And the list names should be consistent with dataTrack names. 
Providing a subset of dataTrack labels is allowed. 
}
  \item{grange2show}{
a GRanges to indicate the plot range
}
  \item{dataTrackName}{
the name of the DataTrack
}
  \item{dataInfo}{
a data matrix or data.frame to show the related sample information, e.g. its expression profile
}
  \item{dataColorMap}{
the color map to plot the dataInfo
}
  \item{dataInfoRange}{
the range of dataInfo to control the range of color map
}
  \item{dataBackground}{
the background color for the data tracks
}
  \item{minHeatmapColumnWidth}{
the minimum width (points) of the heatmap data column 
}
  \item{labelWidth}{
the width of the label, which is the ratio of the entire plot width
}
  \item{gradient}{
the gradient color to show the DataTrack
}
  \item{ncolor}{
the number of color levels
}
  \item{main}{
the title of the plot
}
  \item{newPlot}{
whether to create a new plot or add it to previous plot
}
  \item{sizes}{
the track sizes used by plotTracks function
}
  \item{\dots}{
other parameters used by plotTracks
}
}
\details{
This function is adapted based on the \code{\link[Gviz]{plotTracks}} function in Gviz package. It adds sample labels to the heatmap dataTracks.
}
\value{
Grid viewport layout information
}
\author{
Pan Du
}

\seealso{
 See Also \code{\link[Gviz]{plotTracks}}, \code{\link{heatmapByChromosome}}
}
\examples{

	data(exampleMethyGenoSet)
	if (require(TxDb.Hsapiens.UCSC.hg19.knownGene) && require(Gviz)) {
		## define data track
		exampleMethyGenoSet <- checkChrName(exampleMethyGenoSet)
		dTrack <- DataTrack(range=suppressWarnings(as(locData(exampleMethyGenoSet), 'GRanges')), data=t(exprs(exampleMethyGenoSet)), 
			chromosome='chr21', type='heatmap')		

		## build annotation tracks
		annotationTracks <- buildAnnotationTracks('1826', includeGeneBody = FALSE, genomicFeature = "TxDb.Hsapiens.UCSC.hg19.knownGene", cytobandInfo=NA, CpGInfo=NA)
		trackList <- c(annotationTracks, list(dTrack))
		plotTracksWithDataTrackInfo(trackList, labels=colnames(exampleMethyGenoSet), grange2show = attr(annotationTracks, 'grange2show'))
	}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
\keyword{hplot}
