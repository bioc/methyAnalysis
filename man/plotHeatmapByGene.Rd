\name{plotHeatmapByGene}
\alias{plotHeatmapByGene}
\title{
plot methylation heatmap by genes
}
\description{
plot methylation heatmap by genes
}
\usage{
plotHeatmapByGene(selGene, genoSet,  phenoData = NULL, sortBy=c(NA, 'phenoData', 'data'), includeGeneBody = FALSE, 
	sortByTx = FALSE, CpGInfo = NULL, genomicFeature = NULL, phenoColor = list(gradient=c("green", "black", "red")),
	title.suffix = NULL, addLegend = TRUE, genoSetLegendTitle = NULL, gradient = c("blue", "white", "red"), 
	ncolor = 16, main = NULL, newPlot = TRUE, ylim = NULL, ...)
}
\arguments{
  \item{selGene}{
a Entrez Gene ID 
}
  \item{genoSet}{
a GenoSet object or a list of GenoSet objects
}
  \item{phenoData}{
a data.frame for phenotype information
}
  \item{sortBy}{
whether to sort samples based on the phenoData, cluster of genoSet data or NA (no sorting)
}
  \item{includeGeneBody}{
if FALSE, then only shows the promoter region
}
  \item{sortByTx}{
if TRUE, sort the genoset columns based on Gene Model track. (only valid when the genoset column names are matching transcript IDs.)
}
  \item{CpGInfo}{
a bed file or GRanges for CpG island information
}
  \item{genomicFeature}{
used by buildAnnotationTracks function
}
  \item{phenoColor}{
a list of colors corresponding to phenotype
}
  \item{title.suffix}{
a string attached to the end of the title
}
  \item{addLegend}{
whether to add a legend or not
}
  \item{genoSetLegendTitle}{
title for methylation colorbar legend
}
  \item{gradient}{
the gradient color to show the DataTrack
}
  \item{ncolor}{
the number of color levels
}
  \item{main}{
title of the plot. If it is null, then the Gene Symbol will be the plot title
}
  \item{newPlot}{
whether to create a new plot or add it to previous plot
}
  \item{ylim}{
ylim for the genoSet data, which is also used for plotting the legend.
}
  \item{\dots}{
other parameters used by heatmapByChromosome
}
}
\details{
Function, plotHeatmapByGene, is specifically designed for the methylation data. It plots one gene or genomic range each time. Users can add phenotypes or matched gene expression data to the right panel of the plot. Figure legends can be also added. By default, the plotHeatmapByGene plots methylation Beta-values (in the range of 0 to 1) instead of M-values. Users can set useBetaValue as FALSE if they want to change to M-values.
}
\value{
returns the grid viewport information

}

\author{
Pan Du
}

\seealso{
See also \code{\link{heatmapByChromosome}}
}
\examples{

data('exampleMethyGenoSet') 
if (require(TxDb.Hsapiens.UCSC.hg19.knownGene)) {
  genomicFeature <- 'TxDb.Hsapiens.UCSC.hg19.knownGene'
  selGene <- '1826'
  plotHeatmapByGene(selGene, genoSet=exampleMethyGenoSet, phenoData=pData(exampleMethyGenoSet), genomicFeature=genomicFeature)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
\keyword{hplot}
