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
  MethyGenoSet(locData, exprs, methylated, unmethylated, detection = NULL, 
    pData = NULL, annotation = "", universe = NULL, ...) 
    
  \code{MethyGenoSet} instances are usually created through converting from MethyLumiM object 
  using MethyLumiM2GenoSet function or calling \code{MethyGenoSet} function as shown above. The arguments, 
  locData, exprs, methylated and unmethylated, are required; others can be 
  missing. Please check \code{\link[genoset]{GenoSet}} for more details of other parameters.
    
}
\section{Slots}{
   \describe{
   	  \item{\code{locData}:}{a RangedData object, inherited from \code{\link[genoset]{GenoSet}}}
      \item{\code{assayData}:}{ contains equal dimensional matrices:
 		\code{exprs} (contains the methylation M-value, same as \code{\link[methylumi]{MethyLumiM-class}}),
 		\code{methylated} (contains the methylated probe intensities. Same as \code{\link[methylumi]{MethyLumiM-class}}), 
 		\code{unmethylated} (contains the unmethylated probe intensities. Same as \code{\link[methylumi]{MethyLumiM-class}}),
 		\code{detection} (records the detection p-value of the probe. Same as \code{\link[methylumi]{MethyLumiM-class}}). 
		For more details of \code{assayData}, please see \code{\link[Biobase]{ExpressionSet}} } 
      \item{\code{featureData}:}{See \code{\link[Biobase]{eSet}}}
      \item{\code{phenoData}:}{See \code{\link[Biobase]{eSet}}}
      \item{\code{experimentData}:}{See \code{\link[Biobase]{eSet}}}
      \item{\code{protocolData}:}{See \code{\link[Biobase]{eSet}}}
      \item{\code{annotation}:}{See \code{\link[Biobase]{eSet}}}
      \item{\code{.__classVersion__}:}{See \code{\link[Biobase]{eSet}}}
   	  \item{\code{history}:}{ a data.frame recording the operation history of the MethyGenoSet object.}
   }
}

\section{Methods}{

  \bold{Class-specific methods:}
  \describe{
     \item{\code{exprs(MethyGenoSet)}, \code{exprs(MethyGenoSet,matrix)<-}:}{Access and
       set elements named \code{exprs} in the \code{AssayData-class} slot.}
     \item{\code{methylated(MethyGenoSet)}, \code{methylated(MethyGenoSet)<-}:}{Access and
       set elements named \code{methylated} in the \code{AssayData-class} slot.}
     \item{\code{unmethylated(MethyGenoSet)}, \code{unjmethylated(MethyGenoSet)<-}:}{Access and
       set elements named \code{unmethylated} in the \code{AssayData-class} slot.}
     \item{\code{detection(MethyGenoSet)}, \code{detection(MethyGenoSet)<-}:}{Access and
       set elements named \code{detection} in the \code{AssayData-class} slot.}

    \item{\code{as(methyGenoSet,"MethyLumiM")}}{Coerce objects of
      \code{\link{MethyGenoSet-class}} to \code{MethyLumiM}}
    \item{\code{as(genoSet,"MethyGenoSet")}}{Coerce objects of
      \code{\link{GenoSet-class}} to \code{MethyGenoSet}}

     \item{\code{getHistory(MethyGenoSet)}:}{Access the operation history of \code{MethyGenoSet} object.}
  }


   \bold{Derived from \code{\link[genoset]{GenoSet}}}:
   \describe{
     \item{\code{locData(MethyGenoSet)}:}{return a RangedData object, which contains the chromosome location information}
   }
   
   \bold{Derived from \code{\link[Biobase]{ExpressionSet}}} (For the directly inherited methods, 
	 please see \code{\link[Biobase]{ExpressionSet}} and \code{\link[Biobase]{eSet}}):
   \describe{
     \item{\code{combine(MethyGenoSet,missing)}:}{Combine two MethyGenoSet objects, including \code{history} slot. See \code{\link[Biobase]{eSet}}}
	 \item{\code{exprs(MethyGenoSet)}, \code{exprs(MethyGenoSet,matrix)<-}:}{Access and
	    set elements named \code{exprs} in the \code{AssayData-class} slot.}
	 \item{\code{object[(i,j)}:}{Conduct subsetting of the data in a MethyGenoSet object }
   }

   \bold{Standard generic methods} Please check \code{\link[Biobase]{ExpressionSet}} and \code{\link[Biobase]{eSet}} for other inherited methods, 

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
