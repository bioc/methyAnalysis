\name{export.methyGenoSet}
\alias{export.methyGenoSet}
\title{
Export a MethyGenoSet object to be visualized using external genome browser tools 
}
\description{
Export a MethyGenoSet object to be visualized in IGV, IGB or other tools. Current version supports "gct" or "bw" formats.
}
\usage{
export.methyGenoSet(methyGenoSet, file.format = c("gct", "bw"), exportValue = c("beta", "M", "intensity"), savePrefix = NULL)
}
\arguments{
  \item{methyGenoSet}{
A MethyGenoSet object.
}
  \item{file.format}{
Export file format
}
  \item{exportValue}{
Export methylation values
}
  \item{savePrefix}{
The prefix used in the output filename.
}
}
\details{

An easy way to visualize DNA methylation data is to export the DNA methylation data in certain formats, and visualize these files using some external genome browser tools, like IGV (http://www.broadinstitute.org/igv/) and IGB (http://bioviz.org/igb/index.html). The current implementation of this function supports two output formats: ".gct" and ".bw" files. ".gct" includes all samples in a single file. It is only supported by IGV genome browser. The BigWig format (".bw") is a more general format supported by many visualization tools. Each BigWig file represents one single sample. So it is more flexible for the users only interested in a subset of samples.
}

\value{
Output "gct" (for IGV) or "bw" (BigWig) files

}
\author{
Pan Du
}
\references{
 IGV: http://www.broadinstitute.org/igv/
 IGB: http://bioviz.org/igb/index.html
}
\examples{
	data(exampleMethyGenoSet)
	export.methyGenoSet(exampleMethyGenoSet, file.format='gct', savePrefix='test')
	# export.methyGenoSet(exampleMethyGenoSet, file.format='bw', savePrefix='test')

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
