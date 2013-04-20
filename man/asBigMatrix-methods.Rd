\name{asBigMatrix-methods}
\docType{methods}
\alias{asBigMatrix-methods}
\alias{asBigMatrix,GenoSet-method}
\title{convert the data matrix in the assayData of a GenoSet as BigMatrix}
\description{
 convert the data matrix in the assayData of a GenoSet as BigMatrix
}
\usage{
\S4method{asBigMatrix}{GenoSet}(object, rowInd=NULL, colInd=NULL, nCol=NULL, dimNames=NULL, saveDir='.', savePrefix=NULL, ...)
}

\arguments{
\item{object}{ an object of \code{\link[genoset]{GenoSet}} or its inherited class }
\item{rowInd}{ the subset of row index }
\item{colInd}{ the subset of column index }
\item{nCol}{ the number of columns of the data, which can be larger than the real data dimension. It is designed for adding future data.}
\item{dimNames}{ the dimension names, which is a list of two character vectors (rownames and colnames)}
\item{saveDir}{the parent directory to save the BigMatrix data files}
\item{savePrefix}{the folder name prefix of the directory to save the BigMatrix data files. The fold name will be like this: paste(savePrefix, '_bigmat', sep='') }
\item{\dots}{ optional arguments to \code{\link[bigmemoryExtras]{BigMatrix}} }
}
\details{
The BigMatrix data files will be save in the directory file.path(saveDir, paste(savePrefix, '_bigmat', sep=''))
}
\seealso{
\code{\link[bigmemoryExtras]{BigMatrix}}
}

\keyword{methods}
