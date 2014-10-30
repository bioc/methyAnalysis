##' Identify the CpG-site locations from a genome library
##' 
##' @param bsgenome a BSgenome object or variant name in the genomeLib
##' @param seqnames chromosome names, if missing all chromosomes will be used.
##' @param genomeLib the BSgenome library in Bioconductor
##' @param pattern the sequence pattern to be matched.
##' @return a GRanges of CpG-site locations
##' 
##' @example
##' library(GenomicFeatures)
##' library(BSgenome)
##' seqnames <- paste('chr', c(1:22, 'X', 'Y', 'M'), sep='')
##'  
identifyCpG <- function(bsgenome='Hsapiens', seqnames, genomeLib="BSgenome.Hsapiens.UCSC.hg19", pattern='CG')
{
	if (class(bsgenome) != 'BSgenome') {
		require(genomeLib, character.only=TRUE) || stop(genomeLib, 'package is required!')
		bsgenome <- get(bsgenome)
	} 
	if (missing(seqnames)) seqnames <- GenomicRanges::seqnames(bsgenome)
	
	matchInfo_CG <- matchInfo_GC <- GRanges()
	for (seqname in seqnames) {
		subject <- bsgenome[[seqname]]
		cat("Finding CpG-sites in chromosome", seqname, "...\n")
    cg_matches <- matchPattern(pattern, subject)
 		mLen_cg <- length(cg_matches)
		if (mLen_cg > 0) {
		  match_cg <- GRanges(seqname=rep(seqname, mLen_cg), ranges=IRanges(start=start(cg_matches), end=end(cg_matches)),
		               strand=rep('+', mLen_cg), pattern=rep(pattern, mLen_cg))
			suppressWarnings(matchInfo_CG <- c(matchInfo_CG, match_cg))
		}
    # gc_matches <- matchPattern('GC', subject)
		# mLen_gc <- length(gc_matches)
		# if (mLen_gc > 0) {
		#   match_gc <- GRanges(seqname=rep(seqname, mLen_gc), ranges=IRanges(start=start(gc_matches), end=end(gc_matches)),
		#                strand=rep('+', mLen_gc), pattern=rep('GC', mLen_gc))
		# 	suppressWarnings(matchInfo_GC <- c(matchInfo_GC, match_gc))
		# }
		cat(">>> DONE\n")
	}
	#return(list(CG=matchInfo_CG, GC=matchInfo_GC))
	return(matchInfo_CG)
}

##' filtering the variant calls of Bisulfite converted sequencing data
##'
##' @param seqVariant a VRanges object output by NGS pipeline implemented in HTSeqGenie package
##' @param coverage the genome coverage (a RleList object) output by NGS pipeline
##' @param CpGInfo the precalculated CpG-site information (by identifyCpG function)
##' @param cleanVariant whether to filter those non-CpG with full CT and GA conversion, or non CT and GA variations
##' @param minCoverage minimum coverage for the variants
##' @param convertTh.nonCpG convert rate threshold used by filtering variants (cleanVariant is TRUE)
##'
##' @return a filtered VRanges object with two attributes: 'variantStats' and 'filterSettings'
##' @example
##' # seqVariants.filtered <- filterBisulfiteVariant(seqVariants, coverage, CpGInfo=cpgInfo)
filterBisulfiteVariant <- function(seqVariant, coverage, CpGInfo, cleanVariant=TRUE, minCoverage=1, convertTh.nonCpG=0.9) {
	
	if (!is(seqVariant, 'VRanges')) {
		stop('seqVariant should be a VRanges object!')
	}
	seqVariant <- checkChrName(seqVariant)
	coverage <- checkChrName(coverage)
		
	## For methylation, we only care about CT and GA conversion, so the "alt" should always be T or A
	filterStatus <- alt(seqVariant) %in% c('T', 'A')
	
	count <- altDepth(seqVariant) 
	totalCount <- totalDepth(seqVariant)
	filterStatus <- filterStatus & count > 0
	
	if (minCoverage > 1) {
		filterStatus <- filterStatus & (totalCount >= minCoverage)
		
		## subset the CpGInfo based on coverage
		view.cov <- slice(coverage, lower=minCoverage)
		CpGInfo <- suppressWarnings(CpGInfo[overlapsAny(CpGInfo, view.cov)])		
	} 
	
	## do filtering
	variantFreq <- signif(count / totalCount, 5)
	change <- paste0(ref(seqVariant), alt(seqVariant))
	location <- paste0(seqnames(seqVariant), start(seqVariant), end(seqVariant))
	
	## match with CpG-sites
	# strand(seqVariant) <- '*'
	cpg.match <- suppressWarnings(as.matrix(findOverlaps(seqVariant, CpGInfo, ignore.strand=TRUE)))
	cpg <- rep(0, length(seqVariant))
	cpg[cpg.match[,1]] <- 1
	seqVariant$CpG <- cpg
	## filter those non-informative items: non-CpG with full CT and GA conversion, or non CT and GA variations
	if (cleanVariant) {
		filterStatus <- filterStatus & !(variantFreq >= convertTh.nonCpG & cpg == 0 & change %in% c('CT', 'GA'))
	} 
	## final filtered results
	seqVariant <- seqVariant[filterStatus]

	## add some statistics of filtering
	nonCpG <- which(cpg == 0)
	CpG <- which(cpg == 1)
	converage.variant <- which(totalCount >= minCoverage) 
	## fully converted nonCpG
	converted.nonCpG <- which(variantFreq >= convertTh.nonCpG & cpg == 0 & change %in% c('CT', 'GA'))
	variantStats <- c(totalVariant=length(location), coverageTh=length(converage.variant), nonCpG=length(nonCpG), converted.nonCpG=length(converted.nonCpG), 
			finalVariant=length(seqVariant), CpG=length(CpG))
	variantStats.unique <- c(totalVariant=length(unique(location)), coverageTh=length(unique(location[converage.variant])),  
			nonCpG=length(unique(location[nonCpG])), converted.nonCpG=length(unique(location[converted.nonCpG])),  
			finalVariant=length(unique(location[filterStatus])), CpG=length(unique(location[CpG])))
	variantStats <- rbind(variantStats, variantStats.unique)
	rownames(variantStats) <- c('Variant number', 'Unique location')
	
	attr(seqVariant, 'variantStats') <- variantStats
	attr(seqVariant, 'filterSettings') <- list(cleanVariant=cleanVariant, minCoverage=minCoverage, convertTh.nonCpG=convertTh.nonCpG)
	return(seqVariant)
}


##' Estimate the methylation level (Beta-value) of Methyl-Seq data 
##'
##' @param seqVariant a VRanges object output by NGS pipeline implemented in HTSeqGenie package
##' @param coverage the genome coverage (a RleList object) output by NGS pipeline
##' @param CpGInfo the precalculated CpG-site information (by identifyCpG function)
##' @param mergeStrand whether to merge strand
##' @param cleanVariant whether to filter those non-CpG with full CT and GA conversion, or non CT and GA variations
##' @param minCoverage minimum coverage for the variants
##' @return a GRanges object with the Beta column shows the methylation levels
estimateMethySeq <- function(seqVariant, coverage, CpGInfo=NULL, mergeStrand=TRUE, cleanVariant=TRUE, minCoverage=10) {
	
	if (!is(seqVariant, 'VRanges')) {
		stop('seqVariant should be a VRanges object!')
	}
	seqVariant <- checkChrName(seqVariant)
	coverage <- checkChrName(coverage)
	
	## For methylation, we only care about CT and GA conversion, so the "alt" should always be T or A
	cat('Coverage and TA filtering ...\n')
	filterStatus <- alt(seqVariant) %in% c('T', 'A')
	
	count <- altDepth(seqVariant) 
	totalCount <- totalDepth(seqVariant)
	filterStatus <- filterStatus & count > 0
	
	if (minCoverage > 1) {
		filterStatus <- filterStatus & (totalCount >= minCoverage)
		
		## subset the CpGInfo based on coverage
		view.cov <- slice(coverage, lower=minCoverage)
		CpGInfo <- suppressWarnings(CpGInfo[GenomicRanges::overlapsAny(CpGInfo, view.cov)])		
	} 
	
	## do filtering
	seqVariant <- seqVariant[filterStatus]
	count <- count[filterStatus]
	totalCount <- totalCount[filterStatus]
	variantFreq <- signif(count / totalCount, 5)
	
	change <- paste(ref(seqVariant), alt(seqVariant), sep='')
	# seqInfo <- DataFrame(values(seqVariant), variantFreq=variantFreq, CpG=0, change=change)
	seqInfo <- DataFrame(ref=ref(seqVariant), alt=alt(seqVariant), count=count, totalCount=totalCount, variantFreq=variantFreq, CpG=0, change=change)
	## convert as GRanges because VRanges objects only allow "+" strand
	seqVariant <- as(seqVariant, 'GRanges')
	
	## match with CpG-sites
	cat('Matching CpG-sites ...\n')
	cpg.match <- suppressWarnings(as.matrix(findOverlaps(seqVariant, CpGInfo, ignore.strand=TRUE)))
	
	## check the C-->T conversion rate on the positive strand
	CT.ind <- which(change == 'CT')
	cpg.strand <- NULL
	if (length(CT.ind) > 0) {
		cpg.strand <- rep('+', length(CT.ind))
	}
	
	## check the G-->A conversion rate on the reverse strand
	GA.ind <- which(change == 'GA')
	if (length(GA.ind) > 0) {
		cpg.strand <- c(cpg.strand, rep('-', length(GA.ind)))
	}
	## mark CpG-sites
	convertInd <- c(CT.ind, GA.ind)
	if (!is.null(CpGInfo)) {
		cpg.ind <- intersect(convertInd, cpg.match[,1])
		if (length(cpg.ind) > 0)
			seqInfo[cpg.ind, 'CpG'] <- 1
	}
	if (length(convertInd) > 0)
		strand(seqVariant)[convertInd] <- cpg.strand
	
	## update values of seqVariant
	values(seqVariant) <- seqInfo

	## merge GA and CT if strand information is ignored
	cat('Merge strand ...\n')
	mergeInd <- NULL
	if (mergeStrand) {
		cpg.match <- cpg.match[as.vector(strand(seqVariant)[cpg.match[,1]]) != '*' & seqVariant$CpG[cpg.match[,1]] == 1, ]
		mergeInd <- split(cpg.match[,1], cpg.match[,2])
		mergeInd <- mergeInd[sapply(mergeInd, length) == 2]		
		if (length(mergeInd) > 0) {
			mergeInd <- do.call('rbind', mergeInd)
			## only keep those adjacent CT and GA variant pairs
			mergeInd <- mergeInd[mergeInd[,2] - mergeInd[,1] == 1, , drop=FALSE]
			mergeInd.flat <- as.vector(t(mergeInd))
			CT.ind <- mergeInd.flat[values(seqVariant)$change[mergeInd.flat] == 'CT']
			GA.ind <- mergeInd.flat[values(seqVariant)$change[mergeInd.flat] == 'GA']

			## Adding up the variant counts of AT and GC, with using the higher coverage of the two
			newCount <- count[CT.ind] + count[GA.ind]
			newCount.total <- pmax(totalCount[CT.ind], totalCount[GA.ind])
			## reset newCount as newCount.total if it is larger than the total
			resetInd <- which(newCount > newCount.total)
			if (length(resetInd) > 1)
				newCount[resetInd] <- newCount.total[resetInd]
			values(seqVariant)$count[CT.ind] <- newCount
			values(seqVariant)$totalCount[CT.ind] <- newCount.total
			values(seqVariant)$variantFreq[CT.ind] <- newCount / newCount.total
			strand(seqVariant[mergeInd.flat]) <- '*'
			values(seqVariant)$change[CT.ind] <- 'Both'
		}
	}
	
	## filter those non-informative items: non-CpG with full CT and GA conversion, or non CT and GA variations
	cat('Clean non-informative variants ...\n')
	if (cleanVariant) {
		rmInd <- which((values(seqVariant)$variantFreq > 0.99 & values(seqVariant)$CpG == 0) | strand(seqVariant) == '*')
		if (length(mergeInd) > 0) {
			rmInd <- rmInd[!(rmInd %in% CT.ind)]
		} 
		if (length(rmInd) > 0) seqVariant <- seqVariant[-rmInd]
	} else if (length(mergeInd) > 0) {
		seqVariant <- seqVariant[-GA.ind]
	}
	
	## convert GA as CT on the + strand
	cat('GA to CT conversion ...\n')
	GA.ind <- which(values(seqVariant)$change == 'GA' & values(seqVariant)$CpG == 1)
	if (length(GA.ind) > 0) {
		start(seqVariant)[GA.ind] <- start(seqVariant)[GA.ind] - 1
		end(seqVariant)[GA.ind] <- end(seqVariant)[GA.ind] - 1
		strand(seqVariant[GA.ind]) <- '*'
		values(seqVariant)$ref[GA.ind] <- 'C'
		values(seqVariant)$alt[GA.ind] <- 'T'
	}
	
	## check those CpG sites which has no varaint change (e.g. fully methylated CpG-sites)
	cat('Checking CpG-site without variant change ...\n')
	unchangedCpG <- CpGInfo[!suppressWarnings(overlapsAny(CpGInfo, seqVariant[values(seqVariant)$CpG == 1]))]
	if (length(unchangedCpG) > 0) {
		end(unchangedCpG) <- start(unchangedCpG)
		count.cpg <- getCoverage(unchangedCpG, coverage, startOnly=TRUE, as.GRanges=FALSE)
		# values(unchangedCpG) <- DataFrame(values(unchangedCpG), variantFreq=0, CpG=1, change='CC')
		values(unchangedCpG) <- DataFrame(ref='C', alt='T', count=count.cpg, totalCount=count.cpg, variantFreq=0, CpG=1, change='CC')
		if (minCoverage > 1) {
			unchangedCpG <- unchangedCpG[count.cpg >= minCoverage]
		}
		seqVariant <- suppressWarnings(c(seqVariant, unchangedCpG))
	}
	
	## calculate the Beta-value for CpG-sites
	beta <- rep(NA, length(seqVariant))
	cpg.ind <- (values(seqVariant)$CpG == 1) & (values(seqVariant)$change %in% c('Both', 'CT', 'GA'))
	beta[cpg.ind] <- 1 - values(seqVariant)$variantFreq[cpg.ind]
	values(seqVariant)$Beta <- beta
	
	return(seqVariant)
}


##' get the coverage based on given GRanges object
##'
##' @param grange GRanges object specify the genome locations to get coverage
##' @param coverage the genome coverage (a RleList object) output by NGS pipeline
##' @param startOnly whether to calculate the coverage based on the start location or the average of the entire GRanges element
##' @param as.GRanges whether return a GRanges object or a vector of coverage
##' @return
##'  a vector of coverage or a GRanges object
getCoverage <- function(grange, coverage, startOnly=FALSE, as.GRanges=FALSE) {
	if (!is(grange, 'GRanges')) 
		stop('grange should be a GRanges object!')
	
	allSeqname <- seqnames(grange)
	uniSeqname <- as.character(unique(seqnames(grange)))
	# grangeList <- split(grange, seqnames(grange))
	if (!all(uniSeqname %in% names(coverage))) 
		stop('names of coverage does not match grange object!')
	newGrange <- grange
	values(newGrange) <- DataFrame(coverage=rep(0, length(grange)))
	for (i in seq(uniSeqname)) {
		chr.i <- uniSeqname[i]
		selInd.i <- which(allSeqname == chr.i)
		if (startOnly) {
			cov.i <- as.vector(coverage[[chr.i]][start(grange)[selInd.i]])
		} else {
			cov.i <- sapply(selInd.i, function(ind) {
				mean(coverage[[chr.i]][start(grange[ind]):end(grange[ind])])
			})
		}
		values(newGrange)$coverage[selInd.i] <- cov.i
	}
	if (as.GRanges) {
		return(newGrange)
	} else {
		return(values(newGrange)$coverage)
	}
}


