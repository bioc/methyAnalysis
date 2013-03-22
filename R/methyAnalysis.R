


## convert MethyLumiM class object to GenoSet class object
MethyLumiM2GenoSet <- function(methyLumiM, lib="IlluminaHumanMethylation450k.db") {
	oldFeatureData <- fData(methyLumiM)
	methyLumiM <- addAnnotationInfo(methyLumiM, lib=lib)
	ff <- fData(methyLumiM)
	if (is.null(ff$CHROMOSOME))
		stop('Chromosome information is not available!\n')
		
	## remove those without chromosome location
	rmInd <- which(is.na(ff$CHROMOSOME) | ff$CHROMOSOME == '')
	if (length(rmInd) > 0) {
		ff <- ff[-rmInd,]
		methyLumiM <- methyLumiM[-rmInd,]
	}
	ff$CHROMOSOME <- checkChrName(ff$CHROMOSOME)
	
	hgVersion <- ''
	## retrieve the hgVersion information from the annotation library
	if (require(lib, character.only=TRUE)) {
		dbInfo <- do.call(paste(sub('.db$', '', lib), '_dbInfo', sep=''), list())
		rownames(dbInfo) <- dbInfo$name
		hgVersion <- strsplit(dbInfo['GPSOURCEURL', 'value'], '/')[[1]]
		hgVersion <- hgVersion[length(hgVersion)]
	}
	## create RangedData for location information
	locdata <- RangedData(ranges=IRanges(start=ff$POSITION, width=1, names=featureNames(methyLumiM)), space=ff$CHROMOSOME, universe=hgVersion)

	methyGenoSet <- MethyGenoSet(locData=locdata, pData=pData(methyLumiM), annotation=as.character(lib), exprs=exprs(methyLumiM), methylated=methylated(methyLumiM), 
		unmethylated=unmethylated(methyLumiM), detection=detection(methyLumiM), universe=hgVersion)
	fData(methyGenoSet) <- oldFeatureData[featureNames(methyGenoSet), ]
	methyGenoSet@history <- methyLumiM@history
	
	## set smoothing attributes if exists
	if (!is.null(attr(methyLumiM, 'windowIndex')))
		attr(methyGenoSet, 'windowIndex') <- attr(methyLumiM, 'windowIndex')
	if (!is.null(attr(methyLumiM, 'windowRange')))
		attr(methyGenoSet, 'windowRange') <- attr(methyLumiM, 'windowRange')
	if (!is.null(attr(methyLumiM, 'windowSize')))
		attr(methyGenoSet, 'windowSize') <- attr(methyLumiM, 'windowSize')

	return(methyGenoSet)
}


## smooth the methylation data (MethyLumiM or GenoSet objects) using slide window with fixed windowsize (in bp)
smoothMethyData <- function(methyData, winSize=250, lib='IlluminaHumanMethylation450k.db', ...) {

	if (!is(methyData, 'GenoSet') && !is(methyData, 'MethyLumiM')) {
		stop("methyData should be a GenoSet or MethyLumiM object!")
	}
	
	if (is(methyData, 'MethyLumiM')) {
		chrInfo <- getChrInfo(methyData, lib=lib)
	} else {
		chrInfo <- data.frame(PROBEID=featureNames(methyData), CHROMOSOME=space(locData(methyData)), POSITION=start(methyData), END=end(methyData))
	}
	ratioData <- as.data.frame(exprs(methyData))

	if (is.character(chrInfo$POSITION)) chrInfo$POSITION = as.numeric(chrInfo$POSITION)
	# remove those probes lack of position information
	rmInd <- which(is.na(chrInfo$POSITION) | chrInfo$CHROMOSOME == '' )
	if (length(rmInd) > 0)	{
		ratioData <- ratioData[-rmInd,]
		chrInfo <- chrInfo[-rmInd,]
	}
	
	# Sort the ratioData by chromosome and location
	ord <- order(chrInfo$CHROMOSOME, chrInfo$POSITION, decreasing=FALSE)
	ratioData <- ratioData[ord,]
	chrInfo <- chrInfo[ord,]
	
	# split data by Chromosome 
	ratioData.chrList <- split(ratioData, as.character(chrInfo$CHROMOSOME))
	chrInfoList <- split(chrInfo, as.character(chrInfo$CHROMOSOME))
	windowIndex <- vector(mode='list', length=length(chrInfoList))
	windowRange <- vector(mode='list', length=length(chrInfoList))
	smooth.ratioData <- lapply(seq(ratioData.chrList), function(i) {
 
		ratioData.i <- as.matrix(ratioData.chrList[[i]])
		chrInfo.i <- chrInfoList[[i]]
		chr.i <- chrInfo.i$CHROMOSOME[1]
		cat(paste("Smoothing Chromosome", chr.i, "...\n"))

		# return the index of sorted probes in each slide window
		windowIndex.i <- eval(call(".setupSlidingTests", pos_data=chrInfo.i$POSITION, winSize=winSize))
		names(windowIndex.i) <- chrInfoList[[i]]$PROBEID
		smooth.ratio.i <- sapply(windowIndex.i, function(ind) {
			# average the values in slide window and perform test
			smooth.ij <- colMeans(ratioData.i[ind,,drop=FALSE])
			return(smooth.ij)
		})
		smooth.ratio.i <- t(smooth.ratio.i) 
		windowIndex[[i]] <<- windowIndex.i
		windowRange.i <- sapply(windowIndex.i, function(ind) {
			# average the values in slide window and perform test
			range.ij <- range(chrInfo.i$POSITION[ind])
			return(range.ij)
		})
		windowRange[[i]] <<- t(windowRange.i)
		cat("\n")
		return(smooth.ratio.i)	
	})
	smooth.ratioData <- do.call('rbind', smooth.ratioData)
	windowRange <- do.call('rbind', windowRange)
	colnames(windowRange) <- c('startLocation', 'endLocation')

	methyData <- suppressMessages(methyData[rownames(smooth.ratioData),colnames(smooth.ratioData)])
	exprs(methyData) <- smooth.ratioData[featureNames(methyData),]

	attr(methyData, 'windowIndex') <- windowIndex
	attr(methyData, 'windowRange') <- windowRange
	attr(methyData, 'windowSize') <- winSize
	return(methyData)
}


## export a MethyGenoSet object as a 'gct' (for IGV) or 'bw' (big-wig) file
export.methyGenoSet <- function(methyGenoSet, file.format=c('gct', 'bw'), exportValue=c('beta', 'M'), savePrefix=NULL) {
	
	exportValue <- match.arg(exportValue)
	file.format <- match.arg(file.format)
	## get the annotation version
	hgVersion <- universe(methyGenoSet)
	
	chr <- space(locData(methyGenoSet))
	start <- start(methyGenoSet)
	## Sort the rows of ratios.obj
	methyGenoSet <- methyGenoSet[order(chr, start),]
	
	## check overlap probes and average them
	locationID <- paste(chr, start, sep='_')
	dupInd <- which(duplicated(locationID))
	if (length(dupInd) > 0) {
		dupID <- unique(locationID[dupInd])
		for (dupID.i in dupID) {
			selInd.i <- which(locationID == dupID.i)
			exprs(methyGenoSet)[selInd.i[1],] <- colMeans(exprs(methyGenoSet)[selInd.i,], na.rm=TRUE)
		}
		methyGenoSet <- methyGenoSet[-dupInd,]
	}
	## attach chr prefix in the chromosome names if it does not include it
	if (length(grep('^chr', chr)) == 0) {
		names(locData(methyGenoSet)) <- paste('chr', names(locData(methyGenoSet)), sep='')
	} 
	
	methyData <- assayDataElement(methyGenoSet, 'exprs')
	if (exportValue == 'beta') {
		methyData <- m2beta(methyData) # - 0.5
	} 
	
	## only keep 3 digit to save space
	methyData <- signif(methyData, 3)
	
	if (file.format == 'gct') {
		chrInfo <- data.frame(PROBEID=featureNames(methyGenoSet), CHROMOSOME=space(locData(methyGenoSet)), START=start(methyGenoSet), END=end(methyGenoSet),	stringsAsFactors=FALSE)

		# remove those probes lack of position information
		rmInd <- which(is.na(chrInfo$START))
		if (length(rmInd) > 0)	{
			chrInfo <- chrInfo[-rmInd,]
			methyData <- methyData[-rmInd,]
		}
		
		description <- with(chrInfo, paste("|@", CHROMOSOME, ":", START, "-", END, "|", sep=""))
		gct <- cbind(Name=chrInfo[,'PROBEID'], Description=description, methyData)
		
		## check version hgVersion is included in the filename
		 filename <- paste(savePrefix, "_", exportValue, "_", hgVersion, ".gct", sep='')
		
		cat("#1.2\n", file=filename)
		cat(paste(dim(gct), collapse='\t', sep=''), "\n", file=filename, append=TRUE)
		suppressWarnings(write.table(gct, sep="\t", file=filename, row.names=FALSE, append=TRUE, quote=FALSE))
		return(invisible(filename))

	} else if (file.format == 'bw') {
 		chrInfo <- getChromInfoFromUCSC(hgVersion)
		rownames(chrInfo) <- chrInfo[,'chrom']
 
		samplenames <- sampleNames(methyGenoSet)
		pdata <- pData(methyGenoSet)
		chr.info <- chrInfo(methyGenoSet)
		seq.lengths <- chr.info[,"stop"] - chr.info[,"offset"]
		for (i in 1:ncol(methyData)) {		
			score.i <- methyData[,i]
			cn.data.i <- GRanges(seqnames=space(locData(methyGenoSet)), ranges=IRanges(start=start(methyGenoSet),end=end(methyGenoSet)), strand='*', score=score.i)
			genome(cn.data.i) <- universe(methyGenoSet)
			cn.data.i <- cn.data.i[!is.na(score.i), ]
			if (savePrefix == '' || is.null(savePrefix)) {
				savePrefix.i <- samplenames[i]
			} else {
				savePrefix.i <- paste(savePrefix, samplenames[i], sep='_')
			}
			
			seqlengths(cn.data.i) <- chrInfo[as.character(seqlevels(cn.data.i)), 'length']
			filename.i <- paste(savePrefix.i, "_", exportValue, "_", hgVersion, ".bw", sep='')
			export.bw(cn.data.i, filename.i, dataFormat="auto")
		}
		return(invisible(filename.i))
	}
	
}


# This is the sliding test
## ratios.obj inlcudes c("PROBEID","CHROMOSOME","POSITION"), values 
## winSize is the half slide window size, by default it is 250
## 
detectDMR.slideWin <- function(methyGenoSet, sampleType, winSize=250, testMethod=c('ttest', 'wilcox'), p.adjust.method='fdr', ...) {
	
	testMethod <- match.arg(testMethod)
	if (is(methyGenoSet, 'MethyLumiM')) {
		methyGenoSet <- MethyLumiM2GenoSet(methyGenoSet, ...)
	}
	
	smooth <- TRUE
	if (!is.null(attr(methyGenoSet, 'windowIndex'))) {
		smooth <- FALSE
		if (!is.null(attr(methyGenoSet, 'windowSize'))) {
			if (attr(methyGenoSet, 'windowSize') != winSize) smooth <- TRUE
		}
	} 
	if (smooth) {
		methyGenoSet <- smoothMethyData(methyGenoSet, winSize=winSize, ...)
	}
	windowIndex <- attr(methyGenoSet, 'windowIndex')	
	windowRange <- attr(methyGenoSet, 'windowRange')
	
	chrInfo <- suppressWarnings(as(locData(methyGenoSet), 'GRanges'))
	smoothData <- exprs(methyGenoSet)
	
	## calculate class means
	uniClass <- unique(sampleType)
	selInd1 <- which(sampleType == uniClass[1])
	selInd2 <- which(sampleType == uniClass[2])
	classMean <- data.frame(rowMeans(smoothData[, selInd1, drop=FALSE]), rowMeans(smoothData[, selInd2, drop=FALSE]))
	names(classMean) <- paste('mean_', uniClass, sep='')

	if (!is.factor(sampleType)) sampleType <- as.factor(sampleType)

	## perform t-test
	if (is.function(testMethod)) {
		testInfo <- testMethod(smoothData, sampleType, ... )
	} else {
		if (testMethod == 'ttest') {
			tstats <- rowttests(smoothData, sampleType)
			p.value <- tstats$p.value
			difference <- tstats$dm
		} else if (testMethod == 'wilcox') {
			tmp <- split(as.data.frame(t(smoothData)), sampleType)
			xx <- t(tmp[[1]])
			yy <- t(tmp[[2]])
			p.value <- apply(cbind(xx, yy), 1, function(x) wilcox.test(x[1:ncol(xx)], x[(ncol(xx)+1):length(x)])$p.value, ...)
			difference <- classMean[,1]	- classMean[,2]	
		} else {
			stop('Unsupported "testMethod"!')
		}
		p.adjust <- p.adjust(p.value, method=p.adjust.method)
		testInfo <- data.frame(PROBEID=rownames(smoothData), difference=difference, p.value=p.value, p.adjust=p.adjust)
	}

	startInd <- unlist(lapply(windowIndex, function(x) sapply(x, min)))[rownames(smoothData)]
	endInd <- unlist(lapply(windowIndex, function(x) sapply(x, max)))[rownames(smoothData)]
	testInfo <- data.frame(testInfo, startWinIndex=startInd, endWinIndex=endInd, windowRange, classMean)
	testResult <- chrInfo
	
	## export as a GRanges object
	values(testResult) <- testInfo
	return(testResult)
}


identifySigDMR <- function(detectResult, p.adjust.method="fdr", pValueTh=0.01, fdrTh=pValueTh, diffTh=1, 
			minProbeNum=1, maxGap=2000, minGap=100, oppositeMethylation=FALSE, topNum=NULL) {
	
	if (is(detectResult, 'GRanges')) {
		detectResult.g <- detectResult
		detectResult <- data.frame(CHROMOSOME=as.character(seqnames(detectResult)), POSITION=start(detectResult), as(values(detectResult), 'data.frame'))
		rownames(detectResult) <- names(detectResult.g)
	} else if (is(detectResult, 'data.frame')) {
		if (!all(c("PROBEID","CHROMOSOME","POSITION") %in% names(detectResult)))
			stop(" PROBEID, CHROMOSOME and POSITION are required columns in the detectResult object!")
		# remove those probes lack of position information
		rmInd <- which(is.na(detectResult$POSITION))
		if (length(rmInd) > 0) {
			detectResult <- detectResult[-rmInd,]
		}
		## convert it as a GRanges object
		tmp <- with(detectResult, GRanges(seqnames=CHROMOSOME, ranges=IRanges(start=POSITION, end=POSITION), strand='*'))
		rmInd <- which(colnames(detectResult) %in% c('CHROMOSOME', 'POSITION', 'END', 'width'))
		values(tmp) <- detectResult[,-rmInd, drop=FALSE]
		detectResult.g <- tmp
		names(detectResult.g) <- rownames(detectResult)
	} else {
		stop('detectResult should be either a GRanges object or a data.frame with required columns!')
	}

	## No pValue constraint if we only interested in the top CpG-sites
	if (!is.null(topNum) && !is.na(topNum)) pValueTh <- 1
	
	## if there is no isSignificant column, then DMR should be identified first
	if (is.null(detectResult$isSignificant) || is.null(topNum) || is.na(topNum)) {
		detectResult <- .identifySigProbe(detectResult, p.adjust.method=p.adjust.method, pValueTh=pValueTh, fdrTh=fdrTh, diffTh=diffTh)
	}

	## Select the significant CpG-sites
	sigInd <- which(detectResult$isSignificant)

	## check the methylation direction of two conditions for each probe
	classMeanCol <- grep('^mean_', names(detectResult), value=TRUE)
	# if (length(classMeanCol) == 2) {
	# 	classMean <- detectResult[,classMeanCol]
	# 	if (oppositeMethylation) {
	# 		sigInd <- sigInd[classMean[sigInd,1] * classMean[sigInd,2] < 0]
	# 	}
	# } 

	if (!is.null(topNum) && !is.na(topNum)) {
		ord <- order(detectResult$p.value[sigInd], decreasing=FALSE)
		sigInd <- sigInd[ord]
		if (topNum < length(sigInd)) {
			rmInd <- sigInd[-(1:topNum)]
			detectResult$isSignificant[rmInd] <- FALSE
			sigInd <- sigInd[1:topNum]
		}
	} else {
		if (length(sigInd) == 0) {
			cat('No significant CpG-sites were identified based on current criteria!\n')
			return(NULL)
		}
	}
	sigProbe <- rownames(detectResult)[sigInd]
	detectResult$status <- rep(FALSE, nrow(detectResult))
	detectResult[sigProbe, 'status'] <- TRUE

	if (length(sigProbe) == 0) {
		return(list(sigDMRInfo=NULL, sigDataInfo=NULL))
	}

	## ---------------------------------------
	## identify the boundary of DMR (return a GRanges object)
	scoreFuns <- list(p.value=c(min=min), difference=c(max=function(x) x[which.max(abs(x))]), p.adjust=c(min=min))
	scoreColumns <- c('p.value', 'difference', 'p.adjust', classMeanCol)
	if (is.data.frame(detectResult)) {
		scoreColumns <- scoreColumns[scoreColumns %in% names(detectResult)]
	} else {
		scoreColumns <- scoreColumns[scoreColumns %in% names(values(detectResult))]
	}
	scoreFuns <- scoreFuns[names(scoreFuns) %in% scoreColumns]
	sigDMRInfo.g <- getContinuousRegion(detectResult, scoreColumns=scoreColumns, 
							 scoreFuns=scoreFuns, maxGap=maxGap, minGap=minGap)

	sigDMRInfo <- as(sigDMRInfo.g, 'data.frame')
	rownames(sigDMRInfo) <- names(sigDMRInfo.g)
	
	## filtering based on the minimum number of probes in each DMR
	if (minProbeNum > 1) {
		numProbe <- sigDMRInfo$NumOfProbe
		sigDMRInfo <- sigDMRInfo[numProbe >= minProbeNum, ,drop=FALSE]
		sigDMRInfo.g <- sigDMRInfo.g[numProbe >= minProbeNum]
	}
	
	## check oppositeMethylation direction
	if (length(classMeanCol) == 2 && oppositeMethylation) {
		dmrMean <- sigDMRInfo[, classMeanCol]	
		sigDMRInfo.g <- sigDMRInfo.g[dmrMean[,1] * dmrMean[,2] < 0]
	}

	## only keep the sigDataInfo which overlaps with sigDMRInfo
	sigDataInfo.g <- detectResult.g[sigProbe]
	sigDataInfo.g <- subsetByOverlaps(sigDataInfo.g, sigDMRInfo.g)

	return(list(sigDMRInfo=sigDMRInfo.g, sigDataInfo=sigDataInfo.g))
}


## Get continuous region from discrete probe measurements
# getContinuousRegion <- function(detectResult, scoreColumns=NULL, scoreFuns=mean, maxGap=2000, minGap=100, as.GRanges=TRUE) {
getContinuousRegion <- function(detectResult, scoreColumns=NULL, scoreFuns=c(mean=mean), maxGap=2000, minGap=100) {

	if (is.data.frame(detectResult)) {
		if (!all(c('CHROMOSOME', 'POSITION', 'status') %in% names(detectResult))) {
			stop('CHROMOSOME, POSITION and status columns are required!')
		}
		detectResult.g <- with(detectResult, GRanges(seqnames=CHROMOSOME, ranges=IRanges(start=POSITION, end=POSITION), strand='*'))
		values(detectResult.g) <- detectResult[, !(colnames(detectResult) %in% c('CHROMOSOME', 'POSITION'))]
	} else if (is(detectResult, 'GRanges')) {
		detectResult.g <- detectResult
		detectResult <- data.frame(CHROMOSOME=as.character(seqnames(detectResult)), POSITION=as.numeric(start(detectResult)), as(values(detectResult), 'data.frame'))
	}
	detectResult$POSITION <- as.numeric(detectResult$POSITION)
	## sort in chromsome order
	detectResult <- with(detectResult, detectResult[order(CHROMOSOME, POSITION, decreasing=FALSE),])
	
	if (is.null(detectResult$status)) {
		stop('"status" column is required!')
	} else {
		if (!is.logical(detectResult$status)) {
			warning('"status" column should be a logical vector! "status > 0" was used in calculation!')
			detectResult$status <- detectResult$status > 0
		}
	}
	## set NA as FALSE
	detectResult$status[is.na(detectResult$status)] <- FALSE
	if (!is.null(scoreColumns)) {
		if (!all(scoreColumns %in% colnames(detectResult))) {
			stop('Some "scoreColumns" does not match "detectResult"!')
		}
		if (is.null(scoreFuns)) scoreFuns <- c(mean=mean)
		# if (length(scoreFuns) < length(scoreColumns)) {
		# 	scoreFuns <- rep(scoreFuns, length(scoreColumns))
		# }
	}
	
	## split data by Chromosome 
	detectResult.chrList <- split(detectResult, detectResult$CHROMOSOME)
	
	## ---------------------------------------
	## identify the boundary of DMR
	dmr.list <- lapply(detectResult.chrList, function(detectResult.i) {
		len.significant.i <- length(which(detectResult.i$status))
		if (len.significant.i == 0) return(NULL)
		
		allPosition.i <- detectResult.i$POSITION
		status.i <- as.numeric(detectResult.i$status)
		diff.i <- diff(status.i)
		
		# identify the diff region
		startInd.i <- which(diff.i == 1) + 1
		if (status.i[1] == 1) 
			startInd.i <- c(1, startInd.i)
		endInd.i <- which(diff.i == -1)
		if (length(endInd.i) < length(startInd.i)) 
			endInd.i <- c(endInd.i, length(status.i))
		
		## check the gaps of sparse probe 
		gaps.i <- which(c(FALSE, diff(allPosition.i) > maxGap) & status.i == 1) 
		# the gap should be 1 in both sides
		gaps.i <- gaps.i[status.i[gaps.i + 1] == 1]
		if (length(gaps.i) > 0) {
			startInd.i <- sort(c(startInd.i, gaps.i))
			endInd.i <- sort(c(endInd.i, gaps.i))
		}
		
		if ("startWinIndex" %in% colnames(detectResult.i)) {
			startRegion.i <- detectResult.i[startInd.i, "startWinIndex"]
			endRegion.i <- detectResult.i[endInd.i, "endWinIndex"] 
			startLoc.i <- detectResult.i[startRegion.i, 'POSITION']
			endLoc.i <- detectResult.i[endRegion.i, 'POSITION']
			width.i <- endLoc.i - startLoc.i + 1
			regionSummary.i <- data.frame(CHROMOSOME=detectResult.i[startRegion.i, 'CHROMOSOME'],
				startLocation=startLoc.i, endLocation=endLoc.i, width=width.i, startInd=startRegion.i, endInd=endRegion.i)
		} else {
			startLoc.i <- detectResult.i[startInd.i, "POSITION"]
			endLoc.i <- detectResult.i[endInd.i, "POSITION"] 
			width.i <- endLoc.i - startLoc.i + 1
			regionSummary.i <- data.frame(CHROMOSOME=detectResult.i[startInd.i, 'CHROMOSOME'],
				startLocation=startLoc.i, endLocation=endLoc.i, width=width.i, startInd=startInd.i,
				endInd=endInd.i)
		}
		return(regionSummary.i)
	})
	
	sigDMRInfo <- do.call('rbind', dmr.list)
	if (length(sigDMRInfo) == 0) return(GRanges())

	## convert as GRanges
	tmp <- with(sigDMRInfo, GRanges(seqnames=CHROMOSOME, ranges=IRanges(start=startLocation, end=endLocation), strand='*'))
	rmInd <- which(colnames(sigDMRInfo) %in% c('CHROMOSOME', 'startLocation', 'endLocation', 'width'))
	values(tmp) <- sigDMRInfo[,-rmInd, drop=FALSE]
	sigDMRInfo <- tmp

	## combine overlapping regions
	sigDMRInfo.r <- reduce(sigDMRInfo, min.gapwidth=minGap) 
	
	matchInfo <- GenomicRanges::findOverlaps(detectResult.g, sigDMRInfo.r, select="first")
	ind <- which(!is.na(matchInfo))
	dmr2ind <- split(ind, matchInfo[ind])
	dmr2len <- sapply(dmr2ind, length)
	dmrInfo <- data.frame(NumOfProbe=dmr2len)

	## calculate the mean score of each region	
	if (!is.null(scoreColumns)) {
		if (!all(sapply(scoreFuns, is.function))) {
			if (!all(names(scoreFuns) %in% scoreColumns))
				stop('The names of scoreFuns list should be from scoreColumns!')
		}

		detectValue <- as.matrix(as(values(detectResult.g), 'data.frame'))[,scoreColumns, drop=FALSE]
		detectValue <- matrix(as.numeric(detectValue), nrow=nrow(detectValue), ncol=ncol(detectValue), dimnames=dimnames(detectValue))
		dmr2score <- lapply(dmr2ind, function(ind.i) {
			value.i <- NULL
			## If scoreFuns is a named list, then the named scoreColumns will be treated differently
			if (all(sapply(scoreFuns, is.function))) {
				for (fun.i in scoreFuns) {
					value.i <- cbind(value.i, apply(detectValue[ind.i,,drop=FALSE], 2, fun.i))
				}
			} else {
				for (scoreCol.j in scoreColumns) {
					scoreFuns.j <- scoreFuns[[scoreCol.j]]
					if (is.null(scoreFuns.j)) {
						value.ij <- mean(detectValue[ind.i, scoreCol.j])
						if (!grepl('^mean', scoreCol.j)) {
							names(value.ij) <- paste('mean', scoreCol.j, sep='_')
						} else {
							names(value.ij) <- scoreCol.j
						}
					} else {
						value.ij <- sapply(scoreFuns.j, function(x) x(detectValue[ind.i, scoreCol.j]))
						names(value.ij) <- paste(names(scoreFuns.j), scoreCol.j, sep='_')
					}
					value.i <- c(value.i, value.ij)
				}
			} 
			return(value.i)
		})

		dmr2score <- do.call('rbind', dmr2score)
		rownames(dmr2score) <- names(dmr2ind) 
		avg.score <- matrix(NA, nrow=length(sigDMRInfo.r), ncol=ncol(dmr2score))
		if (is.null(colnames(dmr2score))) {
			scoreFunName <- names(scoreFuns)
			if (is.null(scoreFunName)) scoreFunName <- paste('Fun', 1:length(scoreFuns), sep='')
			colnames(avg.score) <- unlist(lapply(scoreFunName, paste, scoreColumns, sep='_'))
		} else {
			colnames(avg.score) <- colnames(dmr2score)
		}
		avg.score[as.numeric(rownames(dmr2score)),] <- dmr2score
		dmrInfo <- data.frame(dmrInfo, as.data.frame(avg.score))
	}
	values(sigDMRInfo.r) <- dmrInfo
	
	return(sigDMRInfo.r)
}


annotateGRanges <- function(grange, annotationDatabase, CpGInfo=NULL, exons=FALSE, flankRange=0, promoterRange=2000, checkGeneBody=FALSE, EntrezDB='org.Hs.eg.db') {
	
	if (!is(grange, 'GRanges')) {
		stop('grange should be a GRanges object!')
	}
	grange <- checkChrName(grange, addChr=TRUE)	
	if (!is.null(names(grange))) {
		if (any(duplicated(names(grange)))) names(grange) <- NULL
	}
	
	## load human genome information and check overlap with known genes
	if (is.character(annotationDatabase)) {
		if (file.exists(annotationDatabase)) {
			annotationDatabase <- loadFeatures(annotationDatabase)		
		} else if (require(annotationDatabase, character.only=TRUE)) {
			annotationDatabase <- get(annotationDatabase)
		} else {
			stop('Provided annotationDatabase does not exist!')
		}
	} 

	if (is(annotationDatabase, 'TranscriptDb')) {
		tr <- transcripts(annotationDatabase, columns=c('gene_id', 'tx_id', 'tx_name'))		
	} else if (is(annotationDatabase, 'GRanges')) {
		tr <- annotationDatabase
	} else {
		stop('Wrong type of annotationDatabase! Please check help for more details.')
	}
	tr <- checkChrName(tr, addChr=TRUE)	

	if (is.logical(exons)) {
		if (exons && is(annotationDatabase, 'TranscriptDb')) {
			exons <- exons(annotationDatabase, columns=c('exon_id', 'exon_name', 'exon_rank', 'tx_name'))
		} else {
			exons <- NULL
		}
	} else if (!is(exons, 'GRanges')) {
		stop('"exons" should be a logical or GRanges object!')
	}
	
	## filter the transcripts without matching gene ids
	tr <- tr[elementLengths(values(tr)$gene_id) > 0]
	names(tr) <- values(tr)$tx_name
	tr <- checkChrName(tr, addChr=TRUE)
	tx2gene <- sub('GeneID:', '', unlist(values(tr)$gene_id))
	names(tx2gene) <- unlist(values(tr)$tx_name)

	if (is.character(CpGInfo)) {
		CpG.grange <- import.bed(CpGInfo[1], asRangedData=FALSE)
	} else if (is(CpGInfo, 'GRanges')) {
		CpG.grange <- CpGInfo
	} else if (is.null(CpGInfo) || is.na(CpGInfo)) {
		CpG.grange <- NULL
	} else {
		stop('CpGInfo should be a bed file or a GRanges object!')
	}
	if (!is.null(CpGInfo)) CpGInfo <- checkChrName(CpGInfo, addChr=TRUE)

	## expand the interested GRanges
	if (flankRange != 0) {
		grange.ext <- resize(grange, width=width(grange) + 2*flankRange, fix="center")
	} else {
		grange.ext <- grange
	}

	if (is.null(names(grange.ext))) {
		values(grange.ext)$id <- make.unique(paste(seqnames(grange), start(grange), end(grange), sep="_"))
	} else {
		values(grange.ext)$id <- names(grange.ext)
	}

	if (!is.null(values(grange)$Transcript)) {
		txName.known <- values(grange)$Transcript
		txName.known[!(txName.known %in% names(tr))] <- NA
	} else {
		txName.known <- rep(NA, length(grange))
	}
	
	## ----------------------------------------------------------
	## Find the nearest TSS and related Transcript
	tss <- flank(tr, width=-1, start=TRUE)
	if (any(is.na(txName.known))) {
		naInd <- which(is.na(txName.known))
		nearestInfo.tr <- IRanges::as.matrix(nearest(grange[naInd], tr, select="all"))
		
		nearestTx <- tapply(nearestInfo.tr[,2], nearestInfo.tr[,1], function(ind) names(tr)[ind])
		multiMapInd <- which(sapply(nearestTx, length) > 1)
		if (length(multiMapInd) > 0) {
			multiMapInd <- as.numeric(names(nearestTx)[multiMapInd])
			for (i in multiMapInd) {
				nearestTx.i <- nearestTx[[i]]
				nearestTx[[i]] <- nearestTx.i[nearest(grange[naInd][i], tss[nearestTx.i])]
			}
		}
		txName.known[naInd] <- unlist(nearestTx)
		values(grange)$Transcript <- txName.known
	}
	## add related EntrezID and GeneSymbol
	if (is.null(values(grange)$EntrezID)) {
		values(grange)$EntrezID <- tx2gene[values(grange)$Transcript]
	} else {
		values(grange)$EntrezID <- sub('^GeneID:', '', values(grange)$EntrezID)
	}
	if (is.null(values(grange)$GeneSymbol))
		values(grange)$GeneSymbol <- probe2gene(values(grange)$EntrezID, lib=EntrezDB, collapse=TRUE)

	## distance to TSS
	start2tss <- start(grange) - start(tss[txName.known]) 
	end2tss <- end(grange) - end(tss[txName.known]) 
	dist2tss <- pmin(abs(start2tss), abs(end2tss)) * sign(start2tss)
	dist2tss[sign(start2tss) * sign(end2tss) < 0] <- 0
	dist2tss[as.vector(strand(tss[txName.known])) == '-'] <- dist2tss[as.vector(strand(tss[txName.known])) == '-'] * -1
	values(grange)$distance2TSS <- dist2tss

	## nearest transcript based on TSS
	nearestInd.tss <- nearest(grange, tss)
	values(grange)$nearestTx <- values(tss[nearestInd.tss])$tx_name 
	dist2tss.nearestTx <- distanceToNearest(grange, tss[nearestInd.tss])$distance
	
	## mark promoter based on distance to TSS
	promoterStatus <- rep(FALSE, length(dist2tss.nearestTx))
	promoterStatus[dist2tss <= 0 & abs(dist2tss) <= promoterRange] <- TRUE
	values(grange)$PROMOTER <- promoterStatus

	## ----------------------------------------------------------
	## find the overlapping with gene body
	if (checkGeneBody) {
		OL.gene <- findOverlaps(grange.ext, tr)	
				
		dmr2gene <- cbind(values(grange.ext)$id[queryHits(OL.gene)], 
			sapply(values(tr)$gene_id[subjectHits(OL.gene)], paste, collapse=';'))
		dupInd <- duplicated(paste(dmr2gene[,1], dmr2gene[,2], sep='_'))
		if (any(dupInd)) dmr2gene <- dmr2gene[!dupInd,,drop=FALSE]
		dmr2gene.list <- split(dmr2gene[,2], as.factor(dmr2gene[,1]))
		dmr2gene <- sapply(dmr2gene.list, function(x) {
			x <- sub('GeneID:', '', unique(x), ignore.case=TRUE)
			paste(probe2gene(x, lib=EntrezDB, collapse=TRUE), collapse=';')
		})
		dmr2ll <- sapply(dmr2gene.list, function(x) paste(unique(x), collapse=';'))
		
		## map to Tx_name
		dmr2tx <- cbind(values(grange.ext)$id[queryHits(OL.gene)], 
			sapply(values(tr)$tx_name[subjectHits(OL.gene)], paste, collapse=';'))
		dupInd <- duplicated(paste(dmr2tx[,1], dmr2tx[,2], sep='_'))
		if (any(dupInd)) dmr2tx <- dmr2tx[!dupInd,,drop=FALSE]
		dmr2tx.list <- split(dmr2tx[,2], as.factor(dmr2tx[,1]))
		dmr2tx <- sapply(dmr2tx.list, paste, collapse=';')

		mapInfo <- matrix("", nrow=length(grange.ext), ncol=3)
		colnames(mapInfo) <- c("EntrezID_GeneBody", "GeneSymbol_GeneBody", "Tx_GeneBody")
		rownames(mapInfo) <- values(grange.ext)$id
		mapInfo[names(dmr2ll), 'EntrezID_GeneBody'] <- dmr2ll
		mapInfo[names(dmr2gene), 'GeneSymbol_GeneBody'] <- dmr2gene
		mapInfo[names(dmr2tx), 'Tx_GeneBody'] <- dmr2tx
		values(grange) <- data.frame(as(values(grange), 'data.frame'), mapInfo)		
	}
	
	## ----------------------------------------------------------
	## find the overlapping with exons
	if (is(exons, 'GRanges')) {
		exons <- checkChrName(exons, addChr=TRUE)	
		OL.exons <- findOverlaps(grange.ext, exons)	
		hitExons <- exons[subjectHits(OL.exons)]
		if (length(hitExons) > 0) {
			## combine tx_name and exon_rank to get a new exon name
			if (all(is.na(values(hitExons)$exon_name))) {
				ll <- sapply(as.list(values(hitExons)$tx_name), length)
				exon_id <- values(hitExons)$exon_id
				exon_id <- rep(exon_id, ll)
				tx_name <- unlist(values(hitExons)$tx_name)
				exon_rank <- unlist(values(hitExons)$exon_rank)
				exon_name <- paste(tx_name, exon_rank, sep='_')
				exon_name <- split(exon_name, exon_id)
				values(hitExons)$exon_name <- sapply(exon_name, paste, collapse=',')
			}
			dmr2exons <- cbind(values(grange.ext)$id[queryHits(OL.exons)], 
				values(hitExons)$exon_name)

			dupInd <- duplicated(paste(dmr2exons[,1], dmr2exons[,2], sep='_'))
			if (any(dupInd)) {
				dmr2exons <- dmr2exons[!dupInd,,drop=FALSE]
				hitExons <- hitExons[!dupInd]
			}
			dmr2exons.list <- split(dmr2exons[,2], as.factor(dmr2exons[,1]))
			dmr2exons.flat <- sapply(dmr2exons.list, function(x) paste(unique(x), collapse=';'))

			uniId <- unique(values(grange.ext)$id)
			mapInfo <- matrix("", nrow=length(uniId), ncol=2)
			colnames(mapInfo) <- c("Exons", "Exon1")
			rownames(mapInfo) <- uniId
			mapInfo[names(dmr2exons.flat), 'Exons'] <- dmr2exons.flat

			exon1_status <- sapply(as.list(values(hitExons)$exon_rank), function(x) any(x== 1))
			if (any(exon1_status)) {
				id_exon1 <- dmr2exons[exon1_status,1]
				tx_exon1 <- sapply(as.list(values(hitExons)$tx_name[exon1_status]), head, 1)
				mapInfo[id_exon1, 'Exon1'] <- tx_exon1
			}
			values(grange) <- suppressWarnings(data.frame(as(values(grange), 'data.frame'), mapInfo[values(grange.ext)$id,])) 		
		} 
	}	
	
	## ----------------------------------------------------------
	## check overlap with CpG islands
	if (!is.null(CpGInfo)) {
		OL.CpG <- findOverlaps(grange.ext, CpG.grange)	
		dmr2CpG <- cbind(values(grange.ext)$id[queryHits(OL.CpG)], 
			unlist(values(CpG.grange)$name)[subjectHits(OL.CpG)])		

		mapInfo <- rep('', length(grange.ext))
		names(mapInfo) <- values(grange.ext)$id
		mapInfo[dmr2CpG[,1]] <- dmr2CpG[,2]
		values(grange) <- data.frame(as(values(grange), 'data.frame'), CpG_island=mapInfo) 		
	}
	
	return(grange)
}


## input of annotateDMR is the output of identifySigDMR
annotateDMRInfo <- function(DMRInfo, annotationDatabase, CpGInfo=NULL, flankRange=500, promoterRange=2000, EntrezDB='org.Hs.eg.db', as.GRanges=TRUE) {
	
	if (all(c('sigDMRInfo', 'sigDataInfo') %in% names(DMRInfo))) {
		sigDMRInfo <- DMRInfo$sigDMRInfo
		sigDataInfo <- DMRInfo$sigDataInfo
	} else if (is(DMRInfo, 'GRanges')) {
		sigDMRInfo <- DMRInfo
		sigDataInfo <- NULL
	} else {
		stop('DMRInfo should be a GRanges object or a list of GRanges objects (sigDMRInfo and sigDataInfo)!')
	}
	
	## load human genome information and check overlap with known genes
	if (is.character(annotationDatabase)) {
		if (file.exists(annotationDatabase)) {
			annotationDatabase <- loadFeatures(annotationDatabase)		
		} else if (require(annotationDatabase, character.only=TRUE)) {
			annotationDatabase <- get(annotationDatabase)
		} else {
			stop('Provided annotationDatabase does not exist!')
		}
	} 
	
	if (is(annotationDatabase, 'TranscriptDb')) {
		## UCSC only includes reviewed genes!!!
		tr <- transcripts(annotationDatabase, columns=c('gene_id', 'tx_id', 'tx_name'))		
	} else if (is(annotationDatabase, 'GRanges')) {
		tr <- annotationDatabase
	} else {
		stop('Wrong type of annotationDatabase! Please check help for more details.')
	}

	## filter the transcripts without matching gene ids
	tr <- tr[elementLengths(values(tr)$gene_id) > 0]
	names(tr) <- values(tr)$tx_id

	if (is.character(CpGInfo)) {
		CpG.grange <- import.bed(CpGInfo[1], asRangedData=FALSE)
	} else if (is(CpGInfo, 'GRanges')) {
		CpG.grange <- CpGInfo
	} else if (is.null(CpGInfo) || is.na(CpGInfo)) {
		CpG.grange <- NULL
	} else {
		stop('CpGInfo should be a bed file or a GRanges object!')
	}
	
	## expand the transcript region
	## create a GRange object and check overlap
	if (!is.null(sigDMRInfo)) {
		sigDMRInfo <- annotateGRanges(sigDMRInfo, annotationDatabase=tr, CpGInfo=CpG.grange, flankRange=flankRange, promoterRange=promoterRange, EntrezDB=EntrezDB)
	}
	
	if (!is.null(sigDataInfo)) {
		sigDataInfo <- annotateGRanges(sigDataInfo, annotationDatabase=tr, CpGInfo=CpG.grange, flankRange=flankRange, promoterRange=promoterRange, EntrezDB=EntrezDB)
	}
	
	return(list(sigDMRInfo=sigDMRInfo, sigDataInfo=sigDataInfo))
}


export.DMRInfo	<- function(DMRInfo.ann, methyData=NULL, savePrefix='') {
	
	## output the DMR data
	sigDataInfo <- as.data.frame(DMRInfo.ann$sigDataInfo)
	sigProbe <- as.character(sigDataInfo$PROBEID)
	if (!is.null(methyData)) {
		sigDMRdataInfo <- cbind(sigDataInfo, exprs(methyData)[sigProbe, , drop=FALSE])
	} else {
		sigDMRdataInfo <- sigDataInfo
	}
	ord <- order(sigDMRdataInfo$seqnames, sigDMRdataInfo$start, sigDMRdataInfo$end, decreasing=FALSE)
	sigDMRdataInfo <- sigDMRdataInfo[ord,]
	
	if (all(c('startWinIndex', 'endWinIndex') %in% colnames(sigDMRdataInfo))) {
		numProbe <- sigDMRdataInfo[,'endWinIndex'] - sigDMRdataInfo[,'startWinIndex'] + 1
		ind.start <- which(colnames(sigDMRdataInfo) == 'startWinIndex')
		sigDMRdataInfo[,ind.start] <- numProbe
		colnames(sigDMRdataInfo)[ind.start] <- 'numberOfProbe'
		sigDMRdataInfo <- sigDMRdataInfo[,-which(colnames(sigDMRdataInfo) == 'endWinIndex')]
	}
	sigDMRdataInfo <- data.frame(CHROMOSOME=sigDMRdataInfo$seqnames, POSITION=sigDMRdataInfo$start, sigDMRdataInfo[,-c(1:5)])
	## remove columns
	rmCol <- c('isSignificant', 'p.value.adj', 'startLocation', 'endLocation')
	sigDMRdataInfo <- sigDMRdataInfo[, !(colnames(sigDMRdataInfo) %in% rmCol)]
	
	write.csv(sigDMRdataInfo, file=paste('DMRdata_', savePrefix, '_', Sys.Date(), '.csv', sep=''), row.names=FALSE)

	
	## output the sigDMRInfo as a bed file
	export(DMRInfo.ann$sigDMRInfo, paste('DMRInfo_', savePrefix, '_', Sys.Date(), '.bed', sep=''))
	
	## output as a csv file
	sigDMRInfo <- as.data.frame(DMRInfo.ann$sigDMRInfo)
	ord <- order(sigDMRInfo$seqnames, sigDMRInfo$start, sigDMRInfo$end, decreasing=FALSE)
	sigDMRInfo <- sigDMRInfo[ord,]
	if (all(c('startInd', 'endInd') %in% colnames(sigDMRInfo))) {
		numProbe <- sigDMRInfo[,'endInd'] - sigDMRInfo[,'startInd'] + 1
		ind.start <- which(colnames(sigDMRInfo) == 'startInd')
		sigDMRInfo[,ind.start] <- numProbe
		colnames(sigDMRInfo)[ind.start] <- 'numberOfProbe'
		sigDMRInfo <- sigDMRInfo[,-which(colnames(sigDMRInfo) == 'endInd')]
	}
	sigDMRInfo <- data.frame(CHROMOSOME=sigDMRInfo$seqnames, sigDMRInfo[,-1])
	## remove columns
	rmCol <- c('width', 'strand', 'isSignificant', 'p.value.adj')
	sigDMRInfo <- sigDMRInfo[, !(colnames(sigDMRInfo) %in% rmCol)]

	write.csv(sigDMRInfo, file=paste('DMRInfo_', savePrefix, '_', Sys.Date(), '.csv', sep=''), row.names=FALSE)
	

	return(invisible(TRUE))
}


## -----------------------------------------------------
## other internal functions
# detectResult is the return of detectDMR
.identifySigProbe <- function(detectResult, p.adjust.method="fdr", pValueTh=0.01, fdrTh=pValueTh, diffTh=log2(1.5)) {

	if (!is.numeric(pValueTh) && !is.na(pValueTh)) stop("pValueTh needs to be numeric.	Use 1 or NA if you want to turn it off.")
	if (!is.numeric(diffTh) && !is.na(pValueTh)) stop("diffTh needs to be numeric.	Use 0 or NA if you want to turn it off.")
	if (is.na(pValueTh)) pValueTh <- 1
	if (is.na(diffTh)) diffTh <- 0
	annotateColumn <- c("PROBEID","CHROMOSOME","POSITION")
	
	if (!all(annotateColumn %in% names(detectResult)))
		stop(" PROBEID, CHROMOSOME and POSITION are required columns!")
	
	annotateInd <- which(names(detectResult) %in% annotateColumn)
	# multiple testing correction
	p.val.cols <- grep("^p\\.",colnames(detectResult))
	p.val.adj.cols <- grep("p\\.adjust", colnames(detectResult))
	p.val.cols <- p.val.cols[!(p.val.cols %in% p.val.adj.cols)]
	if (length(p.val.cols) > 1) {
		p.value <- apply(detectResult[,p.val.cols], 1, min, na.rm=TRUE)
	} else {
		p.value <- detectResult[,p.val.cols]
	}
	p.value.adj <- p.adjust(p.value, method=p.adjust.method)
	
	if ('difference' %in% colnames(detectResult)) {
		max.diff <- detectResult[,'difference']
	} else {
		differences <- detectResult[, -c(annotateInd, p.val.cols, p.val.adj.cols, which(colnames(detectResult) %in% c("startWinIndex", "endWinIndex"))), drop=FALSE]
		max.diff <- rowMax(abs(as.matrix(differences)))
		max.diff <- max.diff * sign(differences[,1])
	}
	isSignificant <- rep(TRUE, length(max.diff))
	if (pValueTh < 1) isSignificant <- isSignificant & (p.value < pValueTh)
	if (p.adjust.method != 'none' && fdrTh < 1) isSignificant <- isSignificant & (p.value.adj < fdrTh)
	if (diffTh > 0) isSignificant <- isSignificant & (abs(max.diff) > diffTh)

	detectResult <- data.frame(detectResult, p.value.adj=p.value.adj, difference=max.diff, isSignificant=isSignificant)
	
	return(detectResult)
}


# This helper function sets up the positions for the sliding test based on the positions.
# Returns a list of index positions that fit within the sliding window with width winSize
.setupSlidingTests <- function(pos_data, winSize=250) {
	
	len <- length(pos_data)
	probesetList <- as.vector(seq(pos_data), mode='list')
	diff.pos.left <- c(Inf, diff(pos_data))
	growingHeight.left <- diff.pos.left
	growingStatus.left <- which(growingHeight.left <= winSize)
	diff.pos.right <- c(diff.pos.left[-1], Inf)
	growingHeight.right <- diff.pos.right
	growingStatus.right <- which(growingHeight.right <= winSize)
	iter <- 1
	## growing in both left and right directions until it reaches the half-window-size winSize
	while (length(c(growingStatus.left, growingStatus.right)) > 0) {

		growingHeight.left[-(1:iter)] <- growingHeight.left[-(1:iter)] + diff.pos.left[1:(len - iter)]		
		growingHeight.right[1:(len - iter)] <- growingHeight.right[1:(len - iter)] + diff.pos.right[-(1:iter)]
		
		## in left direction
		if (length(growingStatus.left) > 0) {
			probesetList[growingStatus.left] <- lapply(1:length(growingStatus.left), function(i)	 
				c(growingStatus.left[i] - iter, probesetList[[growingStatus.left[i]]]))
		}
		
		## in right direction
		if (length(growingStatus.right) > 0) {
			probesetList[growingStatus.right] <- lapply(1:length(growingStatus.right), function(i) c(probesetList[[growingStatus.right[i]]], growingStatus.right[i] + iter))
		}
		growingStatus.left <- which(growingHeight.left <= winSize)
		growingStatus.right <- which(growingHeight.right <= winSize)	
		
		iter <- iter + 1
	}

	return(probesetList)
}

probe2gene <- function(probe, lib="org.Hs.eg.db", simplify=TRUE, collapse=FALSE) {
	if (length(probe) == 0) return(NULL)
	if (!require(lib, character.only=TRUE)) stop(paste(lib, 'should be installed first!'))
	if (length(grep('.eg.db', lib)) > 0) {
		gene <- getEntrezAnnotation(probe, lib, from='eg', to='symbol')
		if (collapse) {
			gene <- sapply(gene,	paste, collapse=';')
		}				
	} else {
		llib <- paste(gsub('\\.db$', '', lib), 'SYMBOL', sep='')
		mapWithAllProbes <- do.call('toggleProbes', list(x=as.name(llib), value="all"))
		gene <- AnnotationDbi::mget(probe, mapWithAllProbes, ifnotfound=NA)
		if (collapse) {
			gene <- sapply(gene, paste, collapse=';')
		} else {
			if (simplify) gene <- sapply(gene, function(x) x[1])
		}
	}
	return(gene)
}


getEntrezAnnotation <- function(ll=NULL, lib=NULL, from='eg', to='symbol', species=c('human', 'rat', 'mouse', 'yeast', 'fly'), collapse=FALSE) {
	
	from <- toupper(from); to <- toupper(to)
	# from: 'eg', 'accnum', 'alias', 'chr', chrlengths', 'chrloc', 'enzyme', 'genename', 'go', 'map', 'mapcounts', 
	#	'path', 'pfam', 'pmid', 'prosite', 'refseq', 'symbol', 'unigene', 
	species <- match.arg(species)
	
	if (is.null(lib)) {
		lib = switch(species,
			'rat'='org.Rn.eg.db',
			'human'='org.Hs.eg.db',
			'mouse'='org.Mm.eg.db',
			'yeast'='org.Sc.eg.db',
			'fly'='org.Dm.eg.db'
			 )
	}
	if (!require(lib, character.only=TRUE)) {
		stop(paste('Package', lib, 'cannot be loaded!'))
	}

	ll.input <- ll
	na.ind <- which(is.na(ll) | ll == "")
	if (length(na.ind) > 0) {
		ll <- ll[-na.ind]
	}
		
	if (length(to) > 1) {
		allMapping <- NULL
		for (to.i in to) {
			mapping.i <- getEntrezAnnotation(ll, lib, from, to.i, species, collapse=TRUE)
			allMapping <- cbind(allMapping, mapping.i)
		}
		colnames(allMapping) <- to
		rownames(allMapping) <- ll
		return(allMapping)
	}
	libName <- sub('.db', '', lib)
	if (from == 'EG') {
		map <- paste(libName, to, sep='')		
	} else {
		map <- paste(libName, from, '2', to, sep='')		
		if (!exists(map)) {
			if (to == 'EG') {
				map <- paste(libName, from, sep='')
			}
			if (!exists(map)) stop(paste('Cannot find', map, 'in', lib, '!'))
		}
	}

	## 
	if (!is.null(ll)) {
		map <- AnnotationDbi::mget(ll, get(map), ifnotfound=NA)
	} else {
		map <- as.list(get(map))
	}
	if (to == 'EG' & !exists(paste(libName, from, '2', to, sep=''))) {
		len <- sapply(map, length)		
		map <- tapply(rep(names(map), len), unlist(map),	function(x) unique(x))
	}
	
	if (length(na.ind) > 0) {
		map.out <- rep(NA, length(ll.input))
		map.out[-na.ind] <- map
	} else {
		map.out <- map
	}

	if (collapse) {
		map.out <- sapply(map.out,	paste, collapse=';')
	}
	attr(map.out, 'lib') <- lib
	return(map.out)
}




