\name{MethyGenoSet-class}
\docType{class}
\alias{class:MethyGenoSet}
\alias{MethyGenoSet-class}
\alias{MethyGenoSet}

\alias{exprs,MethyGenoSet-method}
\alias{exprs<-,MethyGenoSet-method}
\alias{exprs<-,MethyGenoSet,ANY-method}

\alias{methylated,MethyGenoSet-method}
\alias{methylated<-,MethyGenoSet-method}
\alias{methylated<-,MethyGenoSet,ANY-method}

\alias{unmethylated,MethyGenoSet-method}
\alias{unmethylated<-,MethyGenoSet-method}
\alias{unmethylated<-,MethyGenoSet,ANY-method}

\alias{detection,MethyGenoSet-method}
\alias{detection<-,MethyGenoSet-method}
\alias{detection<-,MethyGenoSet,ANY-method}

\alias{coerce,MethyLumiM,MethyGenoSet-method}
\alias{coerce,MethyGenoSet,MethyLumiM-method}
\alias{coerce,GenoSet,MethyGenoSet-method}

\alias{initialize,MethyGenoSet-method}
\alias{combine,MethyGenoSet,MethyGenoSet-method}
\alias{[,MethyGenoSet-method}
\alias{[,MethyGenoSet,ANY,ANY,ANY-method}
\alias{getHistory,MethyGenoSet-method}

\title{Class MethyGenoSet: contain and describe Illumina Infinium methylation data in GenoSet-class}
\description{This is a class representation for Illumina Infinium methylation microarray data. It directly extends \code{\link[genoset]{GenoSet}}. The purpose of this class is to make the high-density methylation microarray data \code{\link[methylumi]{MethyLumiM-class}} compatible with the Biocoductor infrastructure packages designed for sequencing analysis.}
\section{Extends}{
   Directly extends class \code{\link[genoset]{GenoSet}}.
}

\section{Creating Objects}{
  MethyGenoSet(rowRanges, exprs, methylated, unmethylated, detection = NULL, 
    pData = NULL, annotation = "", universe = NULL, assays=NULL, ...) 
    
  \code{MethyGenoSet} instances are usually created through converting from MethyLumiM object 
  using MethyLumiM2GenoSet function or calling \code{MethyGenoSet} function as shown above. The arguments, 
  rowRanges, exprs, methylated and unmethylated, are required; others can be 
  missing. Please check \code{\link[genoset]{GenoSet}} for more details of other parameters.
    
}
\section{Slots}{
   \describe{
   	  \item{\code{rowRanges}:}{a GRanges object, inherited from \code{\link[SummarizedExperiment]{SummarizedExperiment}}}
      \item{\code{assays}:}{ contains equal dimensional matrices:
 		\code{exprs} (contains the methylation M-value, same as \code{\link[methylumi]{MethyLumiM-class}}),
 		\code{methylated} (contains the methylated probe intensities. Same as \code{\link[methylumi]{MethyLumiM-class}}), 
 		\code{unmethylated} (contains the unmethylated probe intensities. Same as \code{\link[methylumi]{MethyLumiM-class}}),
 		\code{detection} (records the detection p-value of the probe. Same as \code{\link[methylumi]{MethyLumiM-class}}). 
		For more details of \code{assays}, please see \code{\link[SummarizedExperiment]{SummarizedExperiment}} } 
      \item{\code{colData}:}{See \code{\link[SummarizedExperiment]{colData}}}
      \item{\code{annotation}:}{the annotation library of methylation data}
   	  \item{\code{history}:}{ a data.frame recording the operation history of the MethyGenoSet object.}
   }
}

\section{Methods}{

  \bold{Class-specific methods:}
  \describe{
     \item{\code{exprs(MethyGenoSet)}, \code{exprs(MethyGenoSet,matrix)<-}:}{Access and
       set elements named \code{exprs} in the \code{assays} slot.}
     \item{\code{methylated(MethyGenoSet)}, \code{methylated(MethyGenoSet)<-}:}{Access and
       set elements named \code{methylated} in the \code{assays} slot.}
     \item{\code{unmethylated(MethyGenoSet)}, \code{unjmethylated(MethyGenoSet)<-}:}{Access and
       set elements named \code{unmethylated} in the \code{assays} slot.}
     \item{\code{detection(MethyGenoSet)}, \code{detection(MethyGenoSet)<-}:}{Access and
       set elements named \code{detection} in the \code{assays} slot.}

    \item{\code{as(methyGenoSet,"MethyLumiM")}}{Coerce objects of
      \code{\link{MethyGenoSet-class}} to \code{MethyLumiM}}
    \item{\code{as(genoSet,"MethyGenoSet")}}{Coerce objects of
      \code{\link{GenoSet-class}} to \code{MethyGenoSet}}

     \item{\code{getHistory(MethyGenoSet)}:}{Access the operation history of \code{MethyGenoSet} object.}
  }


   \bold{Derived from \code{\link[genoset]{GenoSet}}}:
   \describe{
     \item{\code{rowRanges(MethyGenoSet)}:}{return a GRanges object, which contains the chromosome location information}
   }
   
   \bold{Derived from \code{\link[SummarizedExperiment]{SummarizedExperiment}}} (For the directly inherited methods, 
	 please see \code{\link[genoset]{GenoSet}} and \code{\link[SummarizedExperiment]{SummarizedExperiment}} ):
   \describe{
     \item{\code{combine(MethyGenoSet,missing)}:}{Combine two MethyGenoSet objects, including \code{history} slot}
	 \item{\code{exprs(MethyGenoSet)}, \code{exprs(MethyGenoSet,matrix)<-}:}{Access and
	    set elements named \code{exprs} in the \code{assays} slot.}
	 \item{\code{object[(i,j)}:}{Conduct subsetting of the data in a MethyGenoSet object }
   }

   \bold{Standard generic methods} Please check see \code{\link[genoset]{GenoSet}} and \code{\link[SummarizedExperiment]{SummarizedExperiment}} for other inherited methods, 

}

\author{ Pan Du }
\seealso{\code{\link{MethyLumiM2GenoSet}})
}
\examples{
## load example data
data(exampleMethyGenoSet)
class(exampleMethyGenoSet)

}

\keyword{classes}
