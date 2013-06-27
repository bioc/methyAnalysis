## ---------------------------------------------------------------
## define a new class MethyGenoSet
setClass('MethyGenoSet', 
	representation(history='data.frame'), 
	prototype=list(history=data.frame(
		submitted	 = I(vector()),
		finished	= I(vector()),
		command	 = I(vector()),
		lumiVersion = I(vector())
	)), 
	contains='GenoSet')


setValidity("MethyGenoSet", function(object) 
{
	msg <- Biobase:::validMsg(NULL, Biobase:::isValidVersion(object, "eSet"))
	msg <- Biobase:::validMsg(msg, assayDataValidMembers(assayData(object), c("exprs")))
	msg <- Biobase:::validMsg(msg, assayDataValidMembers(assayData(object), c("methylated")))
	msg <- Biobase:::validMsg(msg, assayDataValidMembers(assayData(object), c("unmethylated")))
	if (is.null(msg)) TRUE else msg
})


## Create MethyGenoSet class
MethyGenoSet <- function(locData, exprs, methylated, unmethylated, detection=NULL, pData=NULL, annotation="", universe=NULL, ...) {
	## convert "RangedData" as "GRanges"
	if (is(locData, 'RangedData')) locData <- as(locData, 'GRanges')
	
	if (is.null(detection)) {
	object <- genoset:::initGenoSet(type="MethyGenoSet", locData=locData, pData=pData, annotation=annotation, universe=universe, exprs=exprs, methylated=methylated, unmethylated=unmethylated, ...)
	} else {
	object <- genoset:::initGenoSet(type="MethyGenoSet", locData=locData, pData=pData, annotation=annotation, universe=universe, exprs=exprs, methylated=methylated, unmethylated=unmethylated, detection=detection, ...)
	}
	return(object)
}

setGeneric("asBigMatrix", function(object, ...) standardGeneric("asBigMatrix"))

setMethod("exprs", signature(object="MethyGenoSet"), function(object) {
	if ('exprs' %in% assayDataElementNames(object)) {
		return(assayDataElement(object,"exprs"))
	} else {
		return(NULL)
	}
})

setReplaceMethod("exprs", signature(object="MethyGenoSet"), function(object, value) {
	if (is.null(value)) {
		assay <- assayData(object)
		if (exists('exprs', envir=assay)) {
			oldMode <- storageMode(assay)
			storageMode(assay) <- 'environment'
			rm(methylated, envir=assay)
			storageMode(assay) <- oldMode
			assayData(object) <- assay
		}
		return(object)
	} else {
		assayDataElementReplace(object, "exprs", value)
	}
})	


setMethod("methylated", signature(object="MethyGenoSet"), function(object) {
	if ('methylated' %in% assayDataElementNames(object)) {
		return(assayDataElement(object,"methylated"))
	} else {
		return(NULL)
	}
})

setReplaceMethod("methylated", signature(object="MethyGenoSet"), function(object, value) {
	if (is.null(value)) {
		assay <- assayData(object)
		if (exists('methylated', envir=assay)) {
			oldMode <- storageMode(assay)
			storageMode(assay) <- 'environment'
			rm(methylated, envir=assay)
			storageMode(assay) <- oldMode
			assayData(object) <- assay
		}
		return(object)
	} else {
		assayDataElementReplace(object, "methylated", value)
	}
})	


setMethod("unmethylated", signature(object="MethyGenoSet"), function(object) {
	if ('unmethylated' %in% assayDataElementNames(object)) {
		return(assayDataElement(object,"unmethylated"))
	} else {
		return(NULL)
	}
})

setReplaceMethod("unmethylated", signature(object="MethyGenoSet"), function(object, value) {
	if (is.null(value)) {
		assay <- assayData(object)
		if (exists('methylated', envir=assay)) {
			oldMode <- storageMode(assay)
			storageMode(assay) <- 'environment'
			rm(unmethylated, envir=assay)
			storageMode(assay) <- oldMode
			assayData(object) <- assay
		}
		return(object)
	} else {
		assayDataElementReplace(object, "unmethylated", value)
	}
})	
	


setMethod("detection", signature(object="MethyGenoSet"), function(object) {
	if ('detection' %in% assayDataElementNames(object)) {
		return(assayDataElement(object,"detection"))
	} else {
		return(NULL)
	}
})

setReplaceMethod("detection", signature(object="MethyGenoSet"), function(object, value) {
	if (is.null(value)) {
		assay <- assayData(object)
		if (exists('detection', envir=assay)) {
			oldMode <- storageMode(assay)
			storageMode(assay) <- 'environment'
			rm(detection, envir=assay)
			storageMode(assay) <- oldMode
			assayData(object) <- assay
		}
		return(object)
	} else {
		assayDataElementReplace(object, "detection", value)
	}
})	


setMethod("getHistory",signature(object="MethyGenoSet"), function(object) object@history)


setMethod("combine", signature=c(x="MethyGenoSet", y="MethyGenoSet"), function(x, y, ...) {

	history.submitted <- as.character(Sys.time())

	if (missing(y)) return(x)
	if (length(list(...)) > 0) 
		return(combine(x, combine(y, ...)))

	## do default processing of 'GenoSet'
	x.comb <- callNextMethod()

	# history tracking
	history.finished <- as.character(Sys.time())
	#history.command <- match.call()
	history.command <- capture.output(print(match.call(combine)))	
	x.comb@history<- rbind(x@history, y@history)
	if (is.null(x.comb@history$lumiVersion) && nrow(x@history) > 0) {
		x.comb@history <- data.frame(x.comb@history, lumiVersion=rep(NA, nrow(x.comb@history)))
	} 
	packageVersion <- paste('methyAnalysis', packageDescription('methyAnalysis')$Version, sep='_')
	x.comb@history<- rbind(x.comb@history, data.frame(submitted=history.submitted,finished=history.finished,command=history.command, lumiVersion=packageVersion))
	return(x.comb)
})


setMethod("[", "MethyGenoSet", function(x, i, j, ..., drop = FALSE)	{

	if (missing(drop)) drop <- FALSE
	history.submitted <- as.character(Sys.time())
	
	## do default processing of 'GenoSet'
	x <- callNextMethod()
	
	ddim <- dim(x)
	if (!missing(i) & !missing(j)) {
			history.command <- paste('Subsetting', ddim[1], 'features and', ddim[2], 'samples.')		
	} else if (!missing(i)) {
			history.command <- paste('Subsetting', ddim[1], 'features.')
	} else if (!missing(j)) {
			history.command <- paste('Subsetting', ddim[2], 'samples.')
	} else {
			return(x)
	}
	
	# history tracking
	history.finished <- as.character(Sys.time())
	if (is.null(x@history$lumiVersion) && nrow(x@history) > 0) {
		x@history <- data.frame(x@history, lumiVersion=rep(NA, nrow(x@history)))
	}
	packageVersion <- paste('methyAnalysis', packageDescription('methyAnalysis')$Version, sep='_')
	x@history<- rbind(x@history, data.frame(submitted=history.submitted,finished=history.finished, command=history.command, lumiVersion=packageVersion))
	
	return(x)
})



## convert MethyLumiM class object to GenoSet class object
setAs("MethyGenoSet", "MethyLumiM", function(from) {
	oldFeatureData <- fData(from)
	locdata <- locData(from)
	chrInfo <- data.frame(CHROMOSOME=as.character(space(locdata)), POSITION=start(locdata))

	methyLumiM <- new('MethyLumiM', phenoData=phenoData(from), annotation=annotation(from), exprs=exprs(from), 
			methylated=methylated(from), unmethylated=unmethylated(from))
	assayData(methyLumiM) <- assayData(from)		
	dataType(methyLumiM) <- 'M'
	if (ncol(oldFeatureData) > 0) {
		ff <- data.frame(chrInfo, oldFeatureData)
	}	else {
		ff <- chrInfo
	}
	fData(methyLumiM) <- ff
	methyLumiM@history <- from@history
	
	## set smoothing attributes if exists
	if (!is.null(attr(from, 'windowIndex')))
		attr(methyLumiM, 'windowIndex') <- attr(from, 'windowIndex')
	if (!is.null(attr(from, 'windowRange')))
		attr(methyLumiM, 'windowRange') <- attr(from, 'windowRange')
	if (!is.null(attr(from, 'windowSize')))
		attr(methyLumiM, 'windowSize') <- attr(from, 'windowSize')

	return(methyLumiM)
})



setAs("GenoSet", "MethyGenoSet", function(from) {
	
	if (!all(c("unmethylated", "methylated") %in% assayDataElementNames(from))) {
		stop("The input should include 'methylated' and 'unmethylated' elements in the assayData slot!\n")
	}
	
	if (is.null(assayDataElement(from,"exprs"))) {
		from <- estimateM(from)
	}
	mm <- assayDataElement(from,"exprs")
	if (!is.null(assayDataElement(from,"detection"))) {
		methyGenoSet <- MethyGenoSet(locData=locData(from), pData=pData(from), annotation=annotation(from), 
				exprs=assayDataElement(from,"exprs"), methylated=assayDataElement(from,"methylated"), unmethylated=assayDataElement(from,"unmethylated"), 
				detection=assayDataElement(from,"detection"))
	} else {
		methyGenoSet <- MethyGenoSet(locData=locData(from), pData=pData(from), annotation=annotation(from), 
				exprs=assayDataElement(from,"exprs"), methylated=assayDataElement(from,"methylated"), unmethylated=assayDataElement(from,"unmethylated"))
	}
	fData(methyGenoSet) <- fData(from)
	return(methyGenoSet)
})


## Example
# library(methyAnalysis)
# data(exampleMethyGenoSet)
# test <- asBigMatrix.genoset(exampleMethyGenoSet, nCol=10)
setMethod('asBigMatrix',
	signature(object='GenoSet'),
	function(object, rowInd=NULL, colInd=NULL, nCol=NULL, dimNames=NULL, saveDir='.', savePrefix=NULL, ...)
{
	## check whether the user just wants a physical copy of the data to a new location
	if (saveDir == '.') saveDir <- getwd()
	if (is.null(savePrefix)) {
		## use the variable name as the bigmatrix directory prefix
		savePrefix <- match.call(asBigMatrix)[['object']]  
	} 	
	saveDir <- file.path(saveDir, paste(savePrefix, 'bigmat', sep='_'))
	
	if (is.null(rowInd) && is.null(colInd) && is.null(nCol) && is.null(dimNames) 
			&& is(assayDataElement(object, assayDataElementNames(object)[1]), 'BigMatrix')) {
		oldDir <- dirname(assayData(object)[[assayDataElementNames(object)[1]]]$datapath)
		if (oldDir != saveDir) {
			if (!file.exists(saveDir)) 	dir.create(saveDir,showWarnings=FALSE)
			sapply(dir(oldDir, full.names=T), file.copy, to=saveDir, overwrite=TRUE, recursive=TRUE)  
		}
		object <- bigmemoryExtras::updateAssayDataElementPaths(object, saveDir)
		return(object)		
	}
	
	if (!is.null(dimNames)) {
		if (!is.list(dimNames) || length(dimNames) != 2) stop("dimNames should be a list with length 2!")
	}
	nRow <- nrow(object)
	if (is.null(nCol)) {
		if (is.null(dimNames)) {
			nCol <- ncol(object)
		} else {
			nCol <- length(dimNames[[2]])
			if (length(nCol) < ncol(object)) stop('The length of dimNames[[2]] should not be less than the columns of the input object!')
		}
	} 
	if (is.null(dimNames)) {
		dimNames <- list(featureNames(object), sampleNames(object))
		## append colnames if nCol is longer than dimNames[[2]]
		if (length(dimNames[[2]]) < nCol) {
			appLen <- nCol - length(dimNames[[2]])
			dimNames[[2]] <- c(dimNames[[2]], paste(rep('unknown', appLen), seq(appLen), sep='.'))
		}
	} else if (length(dimNames[[1]]) < nRow) {
		stop('The length of dimNames[[1]] should be the same as the rows of the input object!')
	}
	subsetMode <- FALSE
	extensionMode <- FALSE 
	if (is.null(rowInd)) {
		rowInd <- 1:nRow
	} else {
		nRow <- length(rowInd)
		dimNames[[1]] <- dimNames[[1]][rowInd]
		subsetMode <- TRUE
	}
	if (is.null(colInd)) {
		colInd <- 1:nCol
	} else {
		nCol <- length(colInd)
		dimNames[[2]] <- dimNames[[2]][colInd]
		subsetMode <- TRUE
	}
	if (nCol > ncol(object)) extensionMode <- TRUE
	
	for (ad.name in assayDataElementNames(object)) {
    matrix.i <- assayDataElement(object, ad.name)
		if (is.null(matrix.i)) next
		backingfile <- file.path(saveDir, ad.name)
		
		if (!is(assayDataElement(object, ad.name), "BigMatrix") && !extensionMode) {
			x.mat <- bigmemoryExtras::BigMatrix(matrix.i[dimNames[[1]], dimNames[[2]]], backingfile=backingfile, nrow=nRow, ncol=nCol, dimnames=dimNames, ...)
		} else {
			x.mat <- bigmemoryExtras::BigMatrix(backingfile=backingfile, nrow=nRow, ncol=nCol, dimnames=dimNames, ...)
			for (i in 1:ncol(matrix.i)) {
				col.i <- colnames(matrix.i)[i]
				x.mat[dimNames[[1]], col.i] <- matrix.i[dimNames[[1]], col.i]
			}
		}
		assayDataElement(object, ad.name) <- x.mat
  }
	object <- bigmemoryExtras::updateAssayDataElementPaths(object, saveDir)

	if (extensionMode || subsetMode) {
		appLen <- nCol - nrow(pData(object))
		if (length(appLen) > 0) {
			pdata <- rbind(as.matrix(pData(object)), matrix(NA, nrow=appLen, ncol=ncol(pData(object))))
			pdata <- as.data.frame(pdata)
			rownames(pdata) <- dimNames[[2]]
		} else {
			pdata <- pData(object)
			if (length(colInd) < nrow(pdata)) pdata <- pdata[colInd,]
		}
		if (class(object) == 'MethyGenoSet') {
			object.new <- MethyGenoSet(locData=locData(object), assayData=assayData(object), pData=pdata)
		} else {
			object.new <- GenoSet(locData=locData(object), assayData=assayData(object), pData=pdata)
		}	
		fData(object.new) <- fData(object)[rowInd,,drop=FALSE]
		object <- object.new
	}
	return(object)
})

