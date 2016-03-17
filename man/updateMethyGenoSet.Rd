\name{updateMethyGenoSet}
\alias{updateMethyGenoSet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Update old version MethyGenoSet objects
}
\description{
Update old version eSet-based MethyGenoSet objects to the SummarizedExperiment-based MethyGenoSet object
}
\usage{
updateMethyGenoSet(methyGenoSet)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{methyGenoSet}{ old version eSet-based MethyGenoSet object}
}
\details{
This function aims to convert the previously saved eSet-based MethyGenoSet objects to the SummarizedExperiment-based MethyGenoSet objects.

Note: Users may see errors when loading the previous saved the eSet-based MethyGenoSet objects. In this case, use "attach" instead of "load". 
}
\value{
A SummarizedExperiment-based MethyGenoSet object
}
\author{
Pan Du
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
