\name{identifyCpG}
\alias{identifyCpG}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Identify the CpG-site locations from a genome library
}
\description{
Identify the CpG-site locations from a genome library
}
\usage{
identifyCpG(bsgenome = "Hsapiens", seqnames, genomeLib = "BSgenome.Hsapiens.UCSC.hg19", pattern = "CG")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bsgenome}{
  a BSgenome object or variant name in the genomeLib
}
  \item{seqnames}{
  chromosome names, if missing all chromosomes will be used.
}
  \item{genomeLib}{
  the BSgenome library in Bioconductor
}
  \item{pattern}{
  the sequence pattern to be matched.
}
}
\value{
a GRanges object with CpG-site locations
}
\author{
Pan Du
}

\seealso{
\code{\link{filterBisulfiteVariant}}
}
\examples{

  # library(GenomicFeatures)
  # library(BSgenome)
  # seqnames <- paste('chr', c(1:22, 'X', 'Y', 'M'), sep='')

}
\keyword{methods}
