\name{plotMethylationHeatmapByGene}
\alias{plotMethylationHeatmapByGene}
\title{
plot methylation heatmap by genes
}
\description{
plot methylation heatmap by genes
}
\usage{
plotMethylationHeatmapByGene(selGene, methyGenoSet, gene2tx = NULL, tx2exon = NULL, expression.tx = NULL, expression.exon = NULL, phenotype = NULL, sortBy=c('expression', 'methylation', 'NA'), renameExon = FALSE, showAllTx = TRUE, useBetaValue = TRUE, includeGeneBody = F, CpGInfo = NULL, genomicFeature = NULL, phenoColor = NULL, th = 0.99, title.suffix = NULL, addLegend = TRUE, gradient = c("blue", "white", "red"), ncolor = 16, main = NULL, newPlot = TRUE, ...)
}
\arguments{
  \item{selGene}{
a vector of EntrezIDs or a list of gene2tx
}
  \item{methyGenoSet}{
a GenoSet object for methylation data
}
  \item{gene2tx}{
a gene to transcript mapping list, used for retrieving expression.tx data
}
  \item{tx2exon}{
a transcript to exon mapping list, used for retrieving expression.exon data
}
  \item{expression.tx}{
an ExpressionSet or data matrix for transcript expression
}
  \item{expression.exon}{
an ExpressionSet or data matrix for exon expression
}
  \item{phenotype}{
a data.frame for phenotype information
}
  \item{sortBy}{
whether to sort samples based on the mean of expression profiles, methylation cluster or NA (no sorting)
}
  \item{renameExon}{
whether to rename exons as "transcript_exon1"
}
  \item{showAllTx}{
whether to show all transcript in gene2tx or just those provided in selGene
}
  \item{useBetaValue}{
whether to use methylation Beta-value in the plot. 
}
  \item{includeGeneBody}{
if FALSE, then only shows the promoter region
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
  \item{th}{
the quantile threshold used to remove outlier, which affects the plot color ranges.
}
  \item{title.suffix}{
a string attached to the end of the title
}
  \item{addLegend}{
whether to add a legend or not
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
  \item{\dots}{
other parameters used by heatmapByChromosome
}
}
\details{
Function, plotMethylationHeatmapByGene, is specifically designed for the methylation data. It plots one gene or genomic range each time. Users can add phenotypes or matched gene expression data to the right panel of the plot. Figure legends can be also added. By default, the plotMethylationHeatmapByGene plots methylation Beta-values (in the range of 0 to 1) instead of M-values. Users can set useBetaValue as FALSE if they want to change to M-values.
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
  plotMethylationHeatmapByGene(selGene, methyGenoSet=exampleMethyGenoSet, phenotype=pData(exampleMethyGenoSet), genomicFeature=genomicFeature)

  ## use different color map for expression data
  es.example <- matrix(runif(ncol(exampleMethyGenoSet), max=10), nrow=1)
  rownames(es.example) <- selGene
  plotMethylationHeatmapByGene(selGene, methyGenoSet=exampleMethyGenoSet, expression.tx=es.example, genomicFeature=genomicFeature, phenoColor=list(gradient=c("green", "black", "red")))
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
\keyword{hplot}
