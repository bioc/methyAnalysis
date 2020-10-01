## ---------------------------------------------------------------
## define a new class MethyGenoSet
setClass('MethyGenoSet',
	representation(history='data.frame', annotation='character'),
	prototype=list(history=data.frame(
		submitted	 = I(vector()),
		finished	= I(vector()),
		command	 = I(vector()),
		lumiVersion = I(vector())
	), annotation=''),
	contains='GenoSet')


setValidity("MethyGenoSet", function(object)
{
	msg <- NULL
	if (!all(c('exprs', 'methylated', 'unmethylated') %in% names(assays(object))))
		msg <- 'exprs, methylated and unmethylated data matrix are required in assayData!'
	if (is.null(msg)) TRUE else msg
})


## Create MethyGenoSet class
MethyGenoSet <- function(rowRanges, exprs=NULL, methylated=NULL, unmethylated=NULL, detection=NULL, pData=NULL, annotation="", universe=NULL, assays=NULL, ...) {
	if (!is.null(universe)) genome(rowRanges) <- universe
	if (!is.null(assays)) {
  	  #if (is.null(rownames(rowRanges))) names(rowRanges) <- rownames(assays)
	  if (!all(c('exprs', 'methylated', 'unmethylated') %in% names(assays))) stop("'exprs', 'methylated' and 'unmethylated are required in assayData!")
	  object <- GenoSet(rowRanges=rowRanges, colData=pData, assays=assays, ...)
	} else {
	  assays <- list(exprs=exprs, methylated=methylated, unmethylated=unmethylated, detecton=detection)
	  object <- GenoSet(rowRanges=rowRanges, colData=pData, assays=assays, ...)
	}
	object <- new('MethyGenoSet', object)
	object@annotation <- annotation
	return(object)
}

setGeneric("asBigMatrix", function(object, ...) standardGeneric("asBigMatrix"))
setGeneric("assayElement", function(object, element, ...) standardGeneric("assayElement"))
setGeneric("assayElement<-", function(object, element, value, ...) standardGeneric("assayElement<-"))

setMethod("exprs", signature(object="MethyGenoSet"), function(object) {
	return(assays(object)$exprs)
})

setReplaceMethod("exprs", signature(object="MethyGenoSet"), function(object, value) {
	assays(object)$exprs <- value
	return(object)
})


setMethod("methylated", signature(object="MethyGenoSet"), function(object) {
	return(assays(object)$methylated)
})

setReplaceMethod("methylated", signature(object="MethyGenoSet"), function(object, value) {
	assays(object)$methylated <- value
	return(object)
})


setMethod("unmethylated", signature(object="MethyGenoSet"), function(object) {
	return(assays(object)$unmethylated)
})

setReplaceMethod("unmethylated", signature(object="MethyGenoSet"), function(object, value) {
	assays(object)$unmethylated <- value
	return(object)
})



setMethod("detection", signature(object="MethyGenoSet"), function(object) {
	return(assays(object)$detection)
})

setReplaceMethod("detection", signature(object="MethyGenoSet"), function(object, value) {
	assays(object)$detection <- value
	return(object)
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
	#oldFeatureData <- fData(from)
	locdata <- rowRanges(from)
	oldFeatureData <- mcols(locdata)
	chrInfo <- data.frame(CHROMOSOME=as.character(genoset::chr(locdata)), POSITION=start(locdata))

	methyLumiM <- new('MethyLumiM', phenoData=as(as.data.frame(colData(from)), 'AnnotatedDataFrame'), annotation=from@annotation, exprs=exprs(from),
			methylated=methylated(from), unmethylated=unmethylated(from))
	assayData(methyLumiM) <- as.list(assays(from))
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

	if (!all(c("unmethylated", "methylated") %in% names(assays(from)))) {
		stop("The input should include 'methylated' and 'unmethylated' elements in the assayData slot!\n")
	}

	if (is.null(assayElement(from,"exprs"))) {
		from <- estimateM(from)
	}
	mm <- assayElement(from,"exprs")
	if (!is.null(assayElement(from,"detection"))) {
		methyGenoSet <- MethyGenoSet(rowRanges=rowRanges(from), pData=pData(from), annotation=annotation(from),
				exprs=assayElement(from,"exprs"), methylated=assayElement(from,"methylated"), unmethylated=assayElement(from,"unmethylated"),
				detection=assayElement(from,"detection"))
	} else {
		methyGenoSet <- MethyGenoSet(rowRanges=rowRanges(from), pData=pData(from), annotation=annotation(from),
				exprs=assayElement(from,"exprs"), methylated=assayElement(from,"methylated"), unmethylated=assayElement(from,"unmethylated"))
	}
    #fData(methyGenoSet) <- fData(from)
    rowRanges(methyGenoSet) <- rowRanges(from)
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
	if (all(rowInd == 1:nrow(object))) rowInd <- NULL
	if (all(colInd == 1:ncol(object))) colInd <- NULL

	if (is.null(rowInd) && is.null(colInd) && is.null(nCol) && is.null(dimNames)
			&& is(assayElement(object, assayNames(object)[1]), 'BigMatrix')) {
	  fieldnames <- ls(assayData(object)[[assayNames(object)[1]]])
	  if ('datapath' %in% fieldnames) {
  		oldDir <- dirname(assays(object)[[assayNames(object)[1]]]$datapath)
	  } else if ('backingfile' %in% fieldnames) {
  		oldDir <- dirname(assays(object)[[assayNames(object)[1]]]$backingfile)
	  }

		if (oldDir == saveDir) {
			return(object)
		} else {
			if (!file.exists(saveDir)) 	dir.create(saveDir,showWarnings=FALSE)
			sapply(dir(oldDir, full.names=T), file.copy, to=saveDir, overwrite=TRUE, recursive=TRUE)
		}
	  ###add if
	if(require(bigmemoryExtras))
	{object <- bigmemoryExtras::updateAssayDataElementPaths(object, saveDir)}
		return(object)
	}

	if (!is.null(dimNames)) {
		if (!is.list(dimNames) | length(dimNames) != 2) stop("dimNames should be a list with length 2!")
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
		dimNames <- list(rownames(object), colnames(object))
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

	for (ad.name in assayNames(object)) {
    matrix.i <- assayElement(object, ad.name)
		if (is.null(matrix.i)) next
		backingfile <- file.path(saveDir, ad.name)

		if (!is(assayElement(object, ad.name), "BigMatrix") && !extensionMode) {
			x.mat <- bigmemoryExtras::BigMatrix(matrix.i[dimNames[[1]], dimNames[[2]]], backingfile=backingfile, nrow=nRow, ncol=nCol, dimnames=dimNames, ...)
		} else {
			x.mat <- bigmemoryExtras::BigMatrix(backingfile=backingfile, nrow=nRow, ncol=nCol, dimnames=dimNames, ...)
			for (i in 1:ncol(matrix.i)) {
				col.i <- colnames(matrix.i)[i]
				x.mat[dimNames[[1]], col.i] <- matrix.i[dimNames[[1]], col.i]
			}
		}
		assayElement(object, ad.name) <- x.mat
  }
	object <- bigmemoryExtras::updateAssayDataElementPaths(object, saveDir)

	if (extensionMode | subsetMode) {
		appLen <- nCol - nrow(pData(object))
		if (length(appLen) > 0) {
			pdata <- rbind(as.matrix(pData(object)), matrix(NA, nrow=appLen, ncol=ncol(pData(object))))
			pdata <- as.data.frame(pdata)
			rownames(pdata) <- dimNames[[2]]
		} else {
			pdata <- pData(object)
			if (length(colInd) < nrow(pdata)) pdata <- pdata[colInd,]
		}
        rowRanges(object) <- rowRanges(object)[rowInd]
	}
	return(object)
})


setMethod("assayElement", signature(object="SummarizedExperiment"), function(object, element) {
	return(assays(object)[[element]])
})

setReplaceMethod("assayElement", signature(object="SummarizedExperiment"), function(object, element, value) {
	assays(object)[[element]] <- value
	return(object)
})


updateMethyGenoSet <- function(methyGenoSet) {
	if (.hasSlot(methyGenoSet, 'assayData')) {
		locdata <- methyGenoSet@locData
		fdata <- methyGenoSet@featureData
		if (!is.null(fdata)) mcols(locdata) <- pData(fdata)
		genoset <- GenoSet(rowRanges=locdata, assays=as.list(methyGenoSet@assayData), colData=pData(methyGenoSet@phenoData))
		genoset <- new('MethyGenoSet', genoset)
		genoset@annotation <- methyGenoSet@annotation
		return(genoset)
	} else {
		return(methyGenoSet)
	}
}
