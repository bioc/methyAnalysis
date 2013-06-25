\name{MethyLumiM2GenoSet}
\alias{MethyLumiM2GenoSet}
\title{
Coerce objects of MethyLumiM-class to MethyGenoSet
}
\description{
    Coerce objects of \code{\link{MethyLumiM-class}} to \code{MethyGenoSet}
}
\usage{
MethyLumiM2GenoSet(methyLumiM, lib = "FDb.InfiniumMethylation.hg19")
}
\arguments{
  \item{methyLumiM}{
a MethyLumiM object
}
  \item{lib}{
lib is a annotation library
}
}
\value{
a MethyGenoSet object
}
\author{
Pan Du
}

\seealso{
 \code{\link{MethyGenoSet}}
}
\examples{

if (require(FDb.InfiniumMethylation.hg19)) {
	data(exampleMethyGenoSet)
	## set as MethyLumiM object
	methyLumiM <- as(exampleMethyGenoSet, 'MethyLumiM')
	## set back as MethyGenoSet object
	methyGenoSet <- MethyLumiM2GenoSet(methyLumiM, lib = "FDb.InfiniumMethylation.hg19")
	class(methyGenoSet)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
