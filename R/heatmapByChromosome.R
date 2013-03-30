
## Example:
## transTrack <- createTranscriptTrack(gene='55432', genomicFeature='TxDb.Hsapiens.UCSC.hg19.knownGene', lib='org.Hs.eg.db', extendRange=c(2000, 2000))
createTranscriptTrack <- function(gene, 
				genomicFeature='TxDb.Hsapiens.UCSC.hg19.knownGene', 
				lib='org.Hs.eg.db', 
				genome='hg19', 
				extendRange=c(2000, 2000), 
				includeOtherGene=FALSE,
				includeGeneBody=TRUE, 
				thinBox_utrOnly=FALSE,
				background.title='gray', 
				fill="#8282d2", ...) {
	
	if (length(extendRange) == 1) {
		extendRange <- rep(extendRange, 2)
	} else if (length(extendRange) == 0) {
		extendRange <- rep(0, 2)
	}
	if (missing(gene)) stop("Please provide 'gene' to plot the data!")
	if (length(gene) > 1) {
		warning('Current version can only support one gene per plot! \n Only the first gene was used!\n')
		gene <- gene[1]
	}

	## convert genomicFeature library as TranscriptDb 
	if (is.character(genomicFeature)) {
		if(!require(genomicFeature, character.only = TRUE)) {
			warning(paste(genomicFeature, 'library is not installed!'))
			genomicFeature <- NULL
		} else {
			genomicFeature <- get(genomicFeature)
		}
	}

	if (is(gene, 'GRanges')) {
		grange2show <- gene
		gene <- names(gene)
		chromosome <- as.character(seqnames(grange2show)[1])
	} else if (is.character(gene)) {
		grange2show <- NULL
		if (grepl('^chr', gene)) {
			chromosome <- gene
		} else if (require(lib, character.only=TRUE)) {
			chromosome <- lookUp(gene, lib, 'CHR')[[1]]
		}
	}
	if (is.na(chromosome)) {
		stop('Cannot find chromosome information for the gene!')
	} else {
		chromosome <- checkChrName(chromosome, addChr=TRUE)	
	}

	## convert TranscriptDb as "GeneRegionTrack" first
	if (is(genomicFeature, 'TranscriptDb')) {
		if (length(grep('^chr', seqlevels(genomicFeature), ignore.case=TRUE)) == 0) {
			options(ucscChromosomeNames=FALSE)
			genomicFeature <- GeneRegionTrack(genomicFeature, chromosome=checkChrName(chromosome, addChr=FALSE), showId=TRUE)
			chromosome(genomicFeature) <- checkChrName(chromosome, addChr=TRUE)
		} else {
			options(ucscChromosomeNames=TRUE)
			genomicFeature <- GeneRegionTrack(genomicFeature, chromosome=chromosome, showId=TRUE)
		}
		
		# ## convert a TranscriptDb as GeneRegionTrack for a selected Gene		
		# genomicFeature <- GeneRegionTrack(genomicFeature, gene=gene, showId=TRUE)
		genomicFeature <- checkChrName(genomicFeature, addChr=TRUE)		
		if (is(grange2show, 'GenomicRanges')) {
			attr(genomicFeature, 'grange2show') <- grange2show
			return(genomicFeature)
		} 
		if (grepl('^chr', gene)) {
			attr(genomicFeature, 'grange2show') <- NULL
			return(genomicFeature)
		}
	}

	## Create transcript track based on provided genomicFeatures
	transTrack <- NULL
	if (is(genomicFeature, 'GeneRegionTrack')) {
		## estimate grange2show and create a transTrack within selected grange2show
		genomicFeature <- checkChrName(genomicFeature, addChr=TRUE)		
		chromosome(genomicFeature) <- checkChrName(chromosome(genomicFeature), addChr=TRUE)	
		
		if (!includeOtherGene || is.null(grange2show)) {
			allGene <- gene(genomicFeature)
			selInd <- which(allGene == gene)
			if (length(selInd) == 0) {
				if (grepl('^GeneID:', gene, ignore.case=TRUE)) {
					gene <- sub('^GeneID:', '', gene)
				} else {
					gene <- paste('GeneID:', gene, sep='')
				}
				selInd <- which(allGene %in% gene)
				if (length(selInd) == 0) {
					warnings('No genes in the selected chromosome range!')
					transTrack <- genomicFeature
					names(transTrack) <- "Gene Model"
					displayPars(transTrack) <- list(background.title=background.title, fill=fill, ...)
					attr(transTrack, 'grange2show') <- NULL
					return(transTrack)
				}
			}
			transTrack <- genomicFeature[selInd]
			selTrans <- ranges(genomicFeature)[selInd]

			if (includeGeneBody) {
				rr <- cbind(start(selTrans), end(selTrans))
				ss <- min(rr) - extendRange[1]
				ee <- max(rr) + extendRange[2]
			} else {
				## only retrieve the region surrounding TSS
				rr <- ranges(genomicFeature[selInd])
				if (!is.null(values(rr)$transcript)) {
					tss <- unlist(sapply(split(rr, values(rr)$transcript), function(x) {
						tss.x <- ifelse(as.character(strand(x)[1]) == '-', max(end(x)), min(start(x)))
					}))
				} else {
					utr5.ind <- which(tolower(feature(genomicFeature)[selInd]) == 'utr5')
					if (length(utr5.ind) > 0) selTrans <- selTrans[utr5.ind]
					tss <- ifelse(as.character(strand(selTrans)) == '-', end(selTrans), start(selTrans))
				}
				ss <- min(tss) - extendRange[1]
				ee <- max(tss) + extendRange[2]
			}
			grange2show <- GRanges(seqnames=chromosome, strand='*', ranges=IRanges(start=ss, end=ee))
		} 
		if (includeOtherGene) {
			if (packageVersion('GenomicRanges') < '1.11.0') {
				transTrack <- genomicFeature[!is.na(match(ranges(genomicFeature), grange2show))]
			} else {
				transTrack <- genomicFeature[overlapsAny(ranges(genomicFeature), grange2show)]
			}
			transTrack <- genomicFeature[transcript(genomicFeature) %in% transcript(transTrack)]
		}
	}
		
	if (is.null(grange2show)) {
		## set the grange2show based on Entrez gene database, e.g., 'org.Hs.eg.db'
		gene <- sub('^GeneID:', '', gene)
		ss <- lookUp(gene, lib, 'CHRLOC')[[1]]
		ee <- lookUp(gene, lib, 'CHRLOCEND')[[1]]
		if (includeGeneBody) {
			rr <- sort(abs(c(ss, ee)), decreasing=FALSE)
			ss <- rr[1] - extendRange[1]
			ee <- rr[2] + extendRange[2]
		} else {
			tss <- ifelse(ss < 0, max(c(abs(ee), abs(ss))), min(c(abs(ee), abs(ss))))
			ss <- tss - extendRange[1]
			ee <- tss + extendRange[2]
		}
		grange2show <- GRanges(seqnames=chromosome, strand='*', ranges=IRanges(start=ss, end=ee))
	}

	## based on the biomart format
	if (is(genomicFeature, "Mart")) {
		## check the existing fields in the mart
		allAttr <- listAttributes(genomicFeature)
		## attributes used by BiomartGeneRegionTrack function 
		attributes <- c("ensembl_gene_id","ensembl_transcript_id","ensembl_exon_id","exon_chrom_start",
										"exon_chrom_end", "rank", "strand", "external_gene_id", "gene_biotype", "chromosome_name")
		if (!all(attributes %in% allAttr[,1])) {
			unmatchInd <- which(!(attributes %in% allAttr[,1]))
			attributes[unmatchInd] <- sub('^ensembl_', '', attributes[unmatchInd])
			unmatchInd <- which(!(attributes %in% allAttr[,1]))
			if (length(unmatchInd) > 0) {
				stop(paste('Attributes', attributes[unmatchInd], 'does not exist in the mart!', collapse=' '))
			} 
			filterNames <- c("chromosome_name", "start", "end")
			filterValues <- c(list(gsub("^chr", "", chromosome), start(grange2show), end(grange2show)))
			ens <- getBM(attributes, filters=filterNames, values=filterValues, mart=genomicFeature, uniqueRows=TRUE)
			colnames(ens) <- c("gene_id","transcript_id","exon_id","start", "end", "rank", "strand", "symbol", 
												 "biotype", "chromosome_name")
			# range <- GRanges(seqnames=.chrName(ens$chromosome_name), ranges=IRanges(start=ens$start, end=ens$end),
			range <- GRanges(seqnames=(ens$chromosome_name), ranges=IRanges(start=ens$start, end=ens$end),
											 strand=ens$strand, feature=as.character(ens$biotype), gene=as.character(ens$gene_id), 
											 exon=as.character(ens$exon_id), id=as.character(ens$exon_id),
											 transcript=as.character(ens$transcript_id), symbol=as.character(ens$symbol),
											 rank=as.numeric(ens$rank))
			# create transTrack
			transTrack <- GeneRegionTrack(ranges=range, start=start(grange2show)[1], end=end(grange2show)[1], genome=genome, chromosome=chromosome, 
						 name="Gene Model", showId=TRUE, background.title=background.title, ...)
			
		} else {
			## BiomartGeneRegionTrack requires several fields which may not exist for some marts
			transTrack <- BiomartGeneRegionTrack(start=start(grange2show)[1], end=end(grange2show)[1],
				chromosome=chromosome, genome=genome, biomart=genomicFeature, showId=TRUE, name='Gene Model', background.title=background.title, ...)
		}		
	}
	## if nothing provided, then create transTrack based on UCSC genome browser
	if (is.null(genomicFeature)) {
		transTrack <- UcscTrack(track="knownGene", trackType="GeneRegionTrack", from=start(grange2show)[1], to=end(grange2show)[1],
			chromosome=chromosome, genome=genome, rstarts="exonStarts", rends="exonEnds", gene="name",
			symbol="name", transcript="name", strand="strand", showId=TRUE, name="UCSC Genes")
	}	
	if (!is.null(transTrack)) {
		names(transTrack) <- "Gene Model"
		displayPars(transTrack) <- list(background.title=background.title, fill=fill, ...)
		attr(transTrack, 'grange2show') <- grange2show
	}
	## only show utr as thinBox
	if (thinBox_utrOnly) {
		thinBoxFeature <- displayPars(transTrack)$thinBoxFeature
		displayPars(transTrack)$thinBoxFeature <- thinBoxFeature[grepl('utr', thinBoxFeature, ignore.case=TRUE)]
	}
	return(transTrack)
}


transcriptDb2GeneRegionTrackByGene <- function(genomicFeature, selGene, extendRange=c(2000, 2000), includeGeneBody=TRUE, includeOtherGene=FALSE, ...) {
	
	geneRange <- NULL
	if (is(selGene, 'GRanges')) {
		geneRange <- selGene
	} else if (is.character(selGene)) {
		if (grepl('^chr', selGene)) {
			chromosome <- selGene
			if (length(grep('^chr', seqlevels(genomicFeature), ignore.case=TRUE)) == 0) {
				options(ucscChromosomeNames=FALSE)
				genomicFeature <- GeneRegionTrack(genomicFeature, chromosome=checkChrName(chromosome, addChr=FALSE), showId=TRUE)
			} else {
				options(ucscChromosomeNames=TRUE)
				genomicFeature <- GeneRegionTrack(genomicFeature, chromosome=chromosome, showId=TRUE)
			}
			attr(genomicFeature, 'grange2show') <- geneRange
			return(genomicFeature)
		} else {
			## get related transcripts
			tr <- transcripts(genomicFeature, vals=list(gene_id=selGene), columns=c('gene_id', 'tx_id', 'tx_name'))	
			if (length(tr) == 0) stop('No matched transcripts were found. Please check the format of Gene ID and genomicFeature!')
			geneRange <- GRanges(seqnames(tr)[1], ranges=IRanges(start=min(start(tr)), end=max(end(tr))), strand=strand(tr)[1])
		}
	}
	## expand the transcript region
	if (!is.null(extendRange)) {
		start(geneRange) <- start(geneRange) - extendRange[1]
		## expand the promoter region
		if (length(extendRange) == 2) {
			end(geneRange) <- end(geneRange) + extendRange[2]
		}
	} 
	genomicFeature <- GeneRegionTrack(genomicFeature, rstarts=start(geneRange), rends=end(geneRange), chromosome=seqnames(geneRange))
	
	if (includeGeneBody) {
		rr <- cbind(start(selTrans), end(selTrans))
		ss <- min(rr) - extendRange[1]
		ee <- max(rr) + extendRange[2]
	} else {
		## only retrieve the region surrounding TSS
		rr <- ranges(genomicFeature)
		if (!is.null(values(rr)$transcript)) {
			tss <- unlist(sapply(split(rr, values(rr)$transcript), function(x) {
				tss.x <- ifelse(as.character(strand(x)[1]) == '-', max(end(x)), min(start(x)))
			}))
		} else {
			utr5.ind <- which(tolower(feature(genomicFeature)[selInd]) == 'utr5')
			if (length(utr5.ind) > 0) selTrans <- selTrans[utr5.ind]
			tss <- ifelse(as.character(strand(selTrans)) == '-', end(selTrans), start(selTrans))
		}
		ss <- min(tss) - extendRange[1]
		ee <- max(tss) + extendRange[2]
	}
	grange2show <- GRanges(seqnames=chromosome, strand='*', ranges=IRanges(start=ss, end=ee))
	if (includeOtherGene) {
		if (packageVersion('GenomicRanges') < '1.11.0') {
			transTrack <- genomicFeature[!is.na(match(ranges(genomicFeature), grange2show))]
		} else {
			transTrack <- genomicFeature[overlapsAny(ranges(genomicFeature), grange2show)]
		}
	}
	
	attr(genomicFeature, 'grange2show') <- geneRange
	return(genomicFeature)
}


## .estimateStackLocation(annotationTracks[[3]], from=ranges[1], to=ranges[2])
.estimateStackLocation <- function(track, from, to, chromosome=NULL, ...) {
	track <- Gviz:::consolidateTrack(track, chromosome=chromosome, from=from, to=to, ...)

	## Now we can subset all the objects in the list to the current boundaries and compute the initial stacking
	track <- Gviz:::subset(track, from=from, to=to)
	track <- Gviz:::setStacks(track, from=from, to=to)
	return(track)
}

## build annotation tracks 
## return value: a list of Tracks (Gviz) with attribute "grange2show"
buildAnnotationTracks <- function(
					gene, 								# Entrez gene ids or a GRanges object with length equals one
					extendRange=c(2000, 2000), # extended range on each side of the	gene
					includeGeneBody=TRUE, # whether to include genebody of the provided gene
					cytobandInfo=NULL,		# cytoband information, 
					CpGInfo=NULL,					# CpG-island information, GRanges or bed file are supported
					genomeAxis=TRUE,			# whether to add genome axis or not
					lib='org.Hs.eg.db',		# gene annotation library
					genome='hg19',				# genome version
					genomicFeature='TxDb.Hsapiens.UCSC.hg19.knownGene',	# genomic features: "TranscriptDb" library or object, "Mart" object
					... ) {
	
	## check parameters
	if (length(gene) > 1) {
		warning('Current version can only support one gene per plot!')
		gene <- gene[1]
	}
	grange2show <- NULL
	if (is(gene, 'GRanges')) {
		chromosome <- as.character(seqnames(gene))
		grange2show <- gene
	} else if (require(lib, character.only=TRUE)) {
		chromosome <- lookUp(gene, lib, 'CHR')[[1]]
	}
	if (is.na(chromosome)) {
		return(NULL)
	}
	
	chromosome <- checkChrName(chromosome, addChr=TRUE)
		
	## create annotation tracks
	# Ideogram (cytoband) track
	allTracks <- iTrack <- NULL
	## check internet connection
	con <- url('http://genome.ucsc.edu')
	internetStatus <- length(readLines(con)) > 0
	close(con)
	if (is.null(cytobandInfo) && internetStatus) {
		iTrack <- IdeogramTrack(genome=genome, chromosome=chromosome)
	} else if (is(cytobandInfo, 'IdeogramTrack')) {
		iTrack <- cytobandInfo
	}
	if (!is.null(iTrack)) {
		allTracks <- c(allTracks, list(iTrack))
	}
	
	# chromosmoe location axis
	if (genomeAxis) {
		gTrack <- GenomeAxisTrack()
		allTracks <- c(allTracks, list(gTrack))
	}

	## CpG-island track
	if (!is.null(CpGInfo)) {
		if (is.character(CpGInfo)) {
			CpGInfo <- import.bed(CpGInfo, asRangedData=F)
		} else if (!is(CpGInfo, 'GRanges') && !is.na(CpGInfo)) {
			stop('CpGInfo should be either a bed file or a GRanges object!')
		}
	} else if (internetStatus) {
		## get CpG-island info from UCSC if CpGInfo is null
		CpGInfo <- getCpGIsland.ucsc(hgVersion=genome)
	} 

	if (is(CpGInfo, 'GRanges')) {
		values(CpGInfo) <- DataFrame(values(CpGInfo), feature='CpG_island', id=1:length(CpGInfo)) 
		cpgTrack <- AnnotationTrack(CpGInfo, chromosome=chromosome, genome=genome, name='CpG Island') # , cex.title=0.8)
		allTracks <- c(allTracks, list(cpgTrack))
	}

	## Create Transcript annotation track
	transTrack	<- createTranscriptTrack(gene, genomicFeature=genomicFeature, lib=lib, 
					extendRange=extendRange, includeGeneBody=includeGeneBody, genome=genome, ...)
	if (!is.null(transTrack)) {
		allTracks <- c(allTracks, list(transTrack))
		## update grange2show
		grange2show <- attr(transTrack, 'grange2show')
	}
	
	attr(allTracks, 'grange2show') <- grange2show
	return(allTracks)
}




heatmapByChromosome <- function(
					genoSet, 				# GenoSet object	
					gene, 								# Entrez gene ids or a GRanges object with length equals one
					annotationTracks=NULL,
					otherTrackList=NULL, 	# other Tracks objects supported by Gviz
					phenoData=NULL, 			# Expression profile
					phonoColorMap = NULL,	# a list of colormaps for every column-side rows
					extendRange=c(2000, 2000), # extended range on each side of the	gene
					includeGeneBody=TRUE, # wether to include genebody of the provided gene
					showFullModel=FALSE, # only valid when includeGeneBody is FALSE
					sortSample=TRUE,			# whether to sort samples based on intensity profiles
					cytobandInfo=NULL,		# cytoband information, 
					CpGInfo=NULL,					# CpG-island information, GRanges or bed file are supported
					genomeAxis=TRUE,			# whether to add genome axis or not
					dataTrackName='Methylation Profile',	# title of the data track
					lib='org.Hs.eg.db',		# gene annotation library
					genome='hg19',				# genome version
					genomicFeature='TxDb.Hsapiens.UCSC.hg19.knownGene',	# genomic features: "TranscriptDb" library or object, "Mart" object
					gradient=c("blue", "white", "red"), # gradient color map
					ncolor=16, 						# number of color levels
					ylim=NULL,
					th=0.99,
					main='',							# title
					... ) {
	
	## check parameters
	if (length(gene) > 1) {
		warning('Current version can only support one gene per plot!')
		gene <- gene[1]
	}
	
	## convert the genoSet  as a list of genoSet for more general purpose
	if (!is.list(genoSet)) {
		genoSetList <- list("Methylation Profile"=genoSet)
	} else {
		genoSetList <- genoSet
	}

	if (!all(sapply(genoSetList, function(x) is(x, 'GenoSet')))) 
		stop('"genoSet" must be a GenoSet object or a list of GenoSet objects.')

	if (is.null(ylim)) {
		ylim <- lapply(genoSetList, function(x) {
			if (th < 1) {
				up.lim <- quantile(abs(assayData(x)$exprs), th, na.rm=TRUE)
				ylim <- c(-up.lim, up.lim)
			} else {
				max.v <- max(abs(assayData(x)$exprs))
				ylim <- c(-max.v, max.v)
			}
			return(ylim)
		})
		ylim <- range(unlist(ylim))
	} 
	
	if (!is.null(phenoData)) {
		rn <- rownames(phenoData)
		## convert the phenoData as numeric matrix
		if (is.data.frame(phenoData) || is.list(phenoData)) {
			phenotypeLevels <- lapply(phenoData, function(x) {
				x <- levels(as.factor(x))
				return(x)
				})
			phenoData <- data.frame(lapply(phenoData, function(x) {
				x <- round(as.numeric(as.factor(x)))
				if (min(x, na.rm=TRUE) <= 0) x <- x - min(x, na.rm=TRUE) + 1
				return(x)
				}), check.names = FALSE)
		}
		rownames(phenoData) <- rn
		phenoData <- as.matrix(phenoData)
		if (is.null(colnames(phenoData)) && ncol(phenoData) == 1)
			colnames(phenoData) <- 'Expression Profile'

		# if (nrow(phenoData) != ncol(genoSetList[[1]])) {
		# 	if (ncol(phenoData) == ncol(genoSetList[[1]])) {
		# 		phenoData <- t(phenoData)
		# 	} else {
		# 		stop('The dimension of "phenoData" does not match with the genoSet!')
		# 	}
		# }
	}
	
	if (is.null(annotationTracks)) {
		annotationTracks <- buildAnnotationTracks(gene=gene, extendRange=extendRange, 
							includeGeneBody=includeGeneBody, cytobandInfo=cytobandInfo, CpGInfo=CpGInfo,
							genomeAxis=genomeAxis, lib=lib,	genome=genome, genomicFeature=genomicFeature, ...)
		#
		if (is.null(annotationTracks))
			stop('"annotationTracks" is missing!')
	}
	allTracks <- annotationTracks
	## get grange2show
	grange2show <- attr(annotationTracks, 'grange2show')
	if (is.null(grange2show)) {
		stop('"grange2show" of the annotationTracks should not be NULL!')
	} 
	chromosome <- as.character(seqnames(grange2show))
	
	## add otherTracks if provided
	if (!is.null(otherTrackList)) {
		if (is.list(otherTrackList)) {
			allTracks <- c(allTracks, otherTrackList)
		} else {
			allTracks <- c(allTracks, list(otherTrackList))
		}
	}
	
	for (i in 1:length(genoSetList)) {
		genoSet <- genoSetList[[i]]
		
		## select related methylation data	
		genoSet <- checkChrName(genoSet, addChr=TRUE)
		grange.data <- suppressWarnings(as(locData(genoSet), 'GRanges'))
		if (packageVersion('GenomicRanges') < '1.11.0') {
			selMethyData <- genoSet[!is.na(match(grange.data, grange2show)),]
		} else {
			selMethyData <- genoSet[overlapsAny(grange.data, grange2show),]
		}
		
		if (nrow(selMethyData) == 0) {
			warning("There is no methylation data exist in the selected grange2show!")
			return(NULL)
		}
		if ((sortSample || length(sortSample) == ncol(selMethyData)) && nrow(selMethyData) > 0 && ncol(selMethyData) > 1) {
			hcr <- hclust(dist(t(assayData(selMethyData)$exprs)))
			ddr <- as.dendrogram(hcr)
			if (length(sortSample) == ncol(selMethyData))
				ddr <- reorder(ddr, sortSample)
			ord <- order.dendrogram(ddr)
			selMethyData <- selMethyData[,ord]
		}

		## define data track
		dTrack <- DataTrack(range=suppressWarnings(as(locData(selMethyData), 'GRanges')), data=t(assayData(selMethyData)$exprs), 
			chromosome=chromosome, name=dataTrackName, type='heatmap',	gradient=gradient, ncolor=ncolor)		
		
		names(dTrack) <- names(genoSetList)[i]
		allTracks <- c(allTracks, list(dTrack))
	}

	if (!includeGeneBody && showFullModel) {
		## remove IdeogramTrack if exists
		allTracks <- allTracks[sapply(allTracks, class) != 'IdeogramTrack']
		
		## estimate the space of tracks
		trackHeights <- .estimateTrackHeight(allTracks, grange2show)
		geneModelTrackInd <- which(sapply(allTracks, class) == 'GeneRegionTrack')
		geneModelHeight <- trackHeights[geneModelTrackInd]

		if (geneModelHeight < 0.1) geneModelHeight <- 0.1
		if (main != '' && !is.null(main)) geneModelHeight <- geneModelHeight + 0.05
		layout.height <- c(geneModelHeight/(1.05+geneModelHeight), 0.05/(1.05+geneModelHeight), 1/(1.05+geneModelHeight))
		pushViewport(viewport(layout=grid.layout(3, 1, heights=layout.height)))

		## plot annotation or phenotype information
		pushViewport(viewport(layout.pos.col=1, layout.pos.row=3))	
		
		## change GenomeAxisTrack 
		axisTrackInd <- which(sapply(allTracks, class) == 'GenomeAxisTrack')
		if (length(axisTrackInd) > 0) {
			displayPars(allTracks[[axisTrackInd]])$littleTicks <- TRUE
		}
		names(allTracks[[geneModelTrackInd]]) <- 'TSS Gene Model'
		## plot Tracks with annotation of DataTrack
		plotInfo <- plotTracksWithDataTrackInfo(allTracks, grange2show=grange2show, dataInfo=phenoData, 
			dataColorMap=phonoColorMap, labelWidth=0.1, gradient=gradient, ncolor=ncolor, ylim=ylim, ...)
		popViewport(1)
		grange2show <- attr(plotInfo, 'grange2show')
		grange2show <- checkChrName(grange2show, addChr=TRUE)

		pushViewport(viewport(layout.pos.col=1, layout.pos.row=1))
		geneModelTrack <- allTracks[[geneModelTrackInd]]
		names(geneModelTrack) <- 'Full Gene Model'

		## make sure the transcript names are shown in right side of plot
		start.range <- min(start(geneModelTrack))
		end.range <- max(end(geneModelTrack))

		labelWidth <- max(nchar(transcript(geneModelTrack))) * getPar(geneModelTrack, 'fontsize') / 2 * getPar(geneModelTrack, 'cex')
		labelWidthRatio <- labelWidth/as.numeric(convertX(unit(0.95, 'npc'), 'points'))
		labelWidth.db <- round((end.range - start.range) * labelWidthRatio)
		# if (labelWidth.points < 0) labelWidth.points <- 0

		start.range <- start.range - labelWidth.db
		# end.range <- max(end(geneModelTrack))
		# start.range <- start.range - (end.range - start.range) / 5
		plotInfo.fullModel <- plotTracks(geneModelTrack, from=start.range, to=end.range, chromosome=chromosome, add=TRUE, main=main, ...)
			
		## plot rectangle of the grange2show
		X0 <- as.numeric(convertX(unit(0,'npc'), 'points'))
		X1 <- as.numeric(convertX(unit(1,'npc'), 'points'))
		Y0 <- as.numeric(convertY(unit(0,'npc'), 'points'))
		Y1 <- as.numeric(convertY(unit(1,'npc'), 'points'))
		totalWidth <- X1 - X0
		totalHeight <- Y1 - Y0
		
		# retrieve the plot coordinates		
		plotLoc.fullModel <- coords(plotInfo.fullModel$title)
		margin.x <- (plotLoc.fullModel[1, 'x1'] - X0)/totalWidth
		plotWidth <- (X1 - plotLoc.fullModel[1, 'x2'])/totalWidth - 2*margin.x
		x0 <- abs(plotLoc.fullModel[1, 'x2'] - X0)/totalWidth + margin.x
		y1 <- 0.01
		plotHeight <- 0.98
		# y1 <- abs(plotLoc.fullModel[1, 'y1'] - Y0)/totalHeight
		# plotHeight <- abs(plotLoc.fullModel[1, 'y2'] - plotLoc.fullModel[1, 'y1'])/totalHeight
		# if (y1 >= 0.01) y1 <- 0.05
		# if (plotHeight + y1 >= 0.95) plotHeight <- 0.95 - y1
		showRegion <- end(grange2show) - start(grange2show)
		x1 <- x0 + plotWidth * (start(grange2show) - start.range)/(end.range - start.range)
		if (x1 < x0) x1 <- x0
		rect.width <- plotWidth * (end(grange2show) - start(grange2show))/(end.range - start.range)
		if (x1 + rect.width > 1) rect.width <- 1 - x1 
		## set minimum width as 0.02
		if (rect.width < 0.02)  {
			x1 <- x1 - (0.02 - rect.width)/4
			rect.width <- 0.02
		}
		grid.rect(x1, y1, width=rect.width, height=plotHeight, gp=gpar(col=2, lwd=1.5, alpha=0.3, lty=2), just=c('left', 'bottom'))
		x1.points <- as.numeric(convertX(unit(x1, 'npc'), 'points'))
		x2.points <- as.numeric(convertX(unit(x1 + rect.width, 'npc'), 'points'))
		y1.points <- y2.points <- as.numeric(convertY(unit(1 - y1, 'npc'), 'points'))

		popViewport(2)
		x1.npc <- convertX(unit(x1.points, 'points'), 'npc')
		x2.npc <- convertX(unit(x2.points, 'points'), 'npc')
		y1.npc <- y2.npc <- convertY(unit(y1.points, 'points'), 'npc')
		## plot dash lines to the zoom-in region
		zoominCor <- coords(plotInfo$title)
		x1.zm <- convertX(unit(zoominCor[1, 'x2'], 'points'), 'npc')
		#y1.zm <- y2.zm <- as.numeric(convertY(unit(zoominCor[1, 'y1'], 'points'), 'npc')) # - 0.01
		y1.zm <- y2.zm <- sum(layout.height[1:2])
		x2.zm <- 1 - plotInfo$labelWidth
	
		grid.lines(x=c(x1.npc, x1.zm), y=1-c(y1.npc,y1.zm), default.units='npc', gp=gpar(col=2, lty=2))
		grid.lines(x=c(x2.npc, x2.zm), y=1-c(y2.npc,y2.zm), default.units='npc', gp=gpar(col=2, lty=2))

		newPlotInfo <- list(fullModel=plotInfo.fullModel, main=plotInfo, layout.height=layout.height)
		plotInfo <- newPlotInfo
	} else {
		## plot Tracks with annotation of DataTrack
		plotInfo <- plotTracksWithDataTrackInfo(allTracks,  grange2show=grange2show, dataInfo=phenoData, 
				dataColorMap=phonoColorMap, labelWidth=0.1, gradient=gradient, ncolor=ncolor, main=main, ylim=ylim, ...)
	}

	return(invisible(plotInfo))
}


## plot methylation heatmap by gene
##	 selGene: a vector of EntrezIDs or a list of gene2tx
##	 methyGenoSet: a GenoSet object for methylation data
##	 gene2tx: a gene to transcript mapping list, used for retrieving expression.tx data
##	 tx2exon: a transcript to exon mapping list, used for retrieving expression.exon data
##	 expression.tx: an ExpressionSet or data matrix for transcript expression
##	 expression.exon: an ExpressionSet or data matrix for exon expression
##	 phenoData: an ExpressionSet or data.frame for phenotype informaiton
##	 sortBy: sort the samples based on expression, methylation or NA (not sort)
##	 renameExon: whether to rename exons as "exon1_transcript"
##	 showAllTx: whether to show all transcript in gene2tx or just those provided in selGene
##	 includeGeneBody: if FALSE, then only shows the promoter region
##	 CpGInfo: a bed file or GRanges for CpG island information
##	 genomicFeature: used by buildAnnotationTracks function
##	 phenoColor: a vector of colors for pheno types.
plotMethylationHeatmapByGene <- function(selGene, methyGenoSet, gene2tx=NULL, tx2exon=NULL, expression.tx=NULL, expression.exon=NULL,  
	phenoData=NULL, sortBy=c('expression', 'methylation', NA), renameExon=FALSE, showAllTx=TRUE, useBetaValue=TRUE, includeGeneBody=FALSE,  
	CpGInfo=NULL, genomicFeature=NULL, phenoColor=list(gradient=c("green", "black", "red")), th=0.99, title.suffix=NULL, addLegend=TRUE, 
	methylationLegendTitle=NULL, expressionLegendTitle='Expression\n(log2-RPKM)', gradient=c("blue", "white", "red"), 
	ncolor=16, main=NULL, newPlot=TRUE, ...) {
		
	sortBy <- as.character(sortBy)
	sortBy <- match.arg(sortBy)
	if(is.na(sortBy)) sortBy <- 'NA'

	## convert gradient as a vector of colors if not
	if (!all(grepl("^#", gradient)) || length(gradient) < 5) {
		gradient <- colorRampPalette(gradient)(ncolor)
	} 
	ncolor <- length(gradient)

	if (!is(methyGenoSet, 'GenoSet')) 
		stop('"methyGenoSet" must be a GenoSet object.')
	
	if (useBetaValue) {
		methyGenoSet <- estimateBeta(methyGenoSet)
		if (th < 1) {
			up.lim <- max(quantile(abs(assayData(methyGenoSet)$exprs), th, na.rm=TRUE) - 0.5, 0.5 - quantile(abs(assayData(methyGenoSet)$exprs), 1-th, na.rm=TRUE))
			ylim <- c(0.5 - up.lim, 0.5 + up.lim)
		} else {
			ylim <- c(0, 1)
		}
	} else {
		if (th < 1) {
			up.lim <- quantile(abs(assayData(methyGenoSet)$exprs), th, na.rm=TRUE)
			ylim <- c(-up.lim, up.lim)
		} else {
			max.v <- max(abs(assayData(methyGenoSet)$exprs))
			ylim <- c(-max.v, max.v)
		}
	}

	if (!is.null(expression.tx)) {
		if (is(expression.tx, 'ExpressionSet')) {
			expMatrix <- exprs(expression.tx)
		} else {
			expMatrix <- expression.tx
		}
		if (max(expMatrix) > 100) {
			minV <- min(expMatrix[expMatrix > 0])
			expMatrix <- log2(expMatrix + minV/2)
		}
		rm(expression.tx)

		if (ncol(expMatrix) != ncol(methyGenoSet)) {
			warnings('Dimensions of expression.tx do not match methyGenoSet!')
			commSample <- intersect(colnames(expMatrix), sampleNames(methyGenoSet))
			expMatrix <- expMatrix[,commSample, drop=FALSE]
			methyGenoSet <- methyGenoSet[,commSample]
		}
		
	} else {
		expMatrix <- NULL
	}
	
	if (!is.null(expression.exon)) {
		if (is(expression.exon, 'ExpressionSet')) {
			expMatrix.exon <- exprs(expression.exon)
		} else {
			expMatrix.exon <- expression.exon
		}
		if (max(expMatrix.exon) > 100) {
			minV <- min(expMatrix.exon[expMatrix.exon > 0])
			expMatrix.exon <- log2(expMatrix.exon + minV/2)
		}
		rm(expression.exon); gc()
	
		if (ncol(expMatrix.exon) != ncol(expMatrix)) {
			warnings('Dimensions of exon data do not match transcript data!')
			commSample <- intersect(colnames(expMatrix.exon), colnames(expMatrix))
			if (ncol(methyGenoSet) == ncol(expMatrix)) {
				ind <- 1:ncol(expMatrix)
				names(ind) <- colnames(expMatrix)
			} else {
				ind <- 1:ncol(expMatrix.exon)
				names(ind) <- colnames(expMatrix.exon)
			}
			expMatrix.exon <- expMatrix.exon[,commSample]
			expMatrix <- expMatrix[,commSample]
			methyGenoSet <- methyGenoSet[,ind[commSample]]
		}
	} else {
		expMatrix.exon <- NULL
	}
	
	if (is.list(selGene)) {
		sigGene2tx <- selGene
		selGene <- names(selGene)
	} else {
		showAllTx <- FALSE
		if (!is.null(gene2tx)) {
			sigGene2tx <- gene2tx[selGene]
		} else {
			sigGene2tx <- as.list(selGene)
			showAllTx <- FALSE
		}
	}

	phenotypeLevels <- NULL
	if (!is.null(phenoData)) {
		if (is(phenoData, "ExpressionSet")) {
			phenoData <- exprs(phenoData)
		} 
		rn <- rownames(phenoData)
		## convert the phenoData as numeric matrix
		if (is.data.frame(phenoData) || is.list(phenoData)) {
			phenotypeLevels <- lapply(phenoData, function(x) {
				x <- levels(as.factor(x))
				return(x)
				})
			phenoData <- data.frame(lapply(phenoData, function(x) {
				x <- round(as.numeric(as.factor(x)))
				if (min(x, na.rm=TRUE) <= 0) x <- x - min(x, na.rm=TRUE) + 1
				return(x)
				}), check.names = FALSE)
		}
		rownames(phenoData) <- rn

		otherPhenoName <- colnames(phenoData)
		otherPhenoName <- otherPhenoName[!(otherPhenoName %in% names(phenoColor))]
		if (length(otherPhenoName) > 0) {
			allPhenoName <- names(phenoColor)
			for (phenoName.i in otherPhenoName) {
				maxLevel.i <- max(phenoData[[phenoName.i]], na.rm=TRUE)
				if (maxLevel.i <= 6) {
					phenoColor <- c(phenoColor, list(1:maxLevel.i))
					allPhenoName <- c(allPhenoName, phenoName.i)
				}
			}
			names(phenoColor) <- allPhenoName
		}		
	}

	plotResult <- lapply(1:length(sigGene2tx), function(i) {
	
		gene.i <- selGene[i] # '728591' # '63951'
		symbol <- unlist(lookUp(gene.i, 'org.Hs.eg.db', 'SYMBOL'))
		if (is.na(symbol) || (symbol == 'NA')) return(NULL)
		
		annotationTracks <- buildAnnotationTracks(gene=gene.i, includeGeneBody=includeGeneBody, CpGInfo=CpGInfo, genomicFeature=genomicFeature, ...)
	
		sigTx <- sigGene2tx[[i]]
		if (is.null(sigTx)) return(NULL)
		if (showAllTx) {
			selTx <- gene2tx[[gene.i]]
			selTx <- c(sigTx, selTx[!(selTx %in% sigTx)])
		} else {
			selTx <- sigTx
		}
	
		## sort the transcript based on annotation track
		geneRegionTrack <- annotationTracks[sapply(annotationTracks, class) == 'GeneRegionTrack'][[1]]
		## estimate the order of transcripts in the geneRegionTrack
		grange2show <- attr(annotationTracks, 'grange2show')
		if (is.null(grange2show)) {
			warnings('No region to plot!')
			return(NULL)
		}
		
		grange2show <- checkChrName(grange2show, addChr=TRUE)
		chromosome <- as.character(seqnames(grange2show))[1]
		geneRegionTrack <- .estimateStackLocation(geneRegionTrack, from=start(grange2show)[1], to=end(grange2show)[1], chromosome=chromosome)
		annTx <- split(values(geneRegionTrack)$transcript, stacks(geneRegionTrack))
		annTx <- unique(unlist(annTx))

		annTx.ind <- 1:length(annTx)
		names(annTx.ind) <- annTx
		orderedTx <- selTx[selTx %in% annTx]
		if (length(orderedTx) > 0) {
			orderedTx <- orderedTx[order(annTx.ind[orderedTx], decreasing=FALSE)]
			orderedTx <- c(orderedTx, selTx[!(selTx %in% annTx)])
			selTx <- orderedTx
		}
		if (!any(selTx %in% rownames(expMatrix))) {
			if (gene.i %in% rownames(expMatrix)) {
				expProfile <- t(expMatrix[gene.i,,drop=FALSE])
			} else {
				expProfile <- NULL
			}
		} else {
			selTx <- selTx[selTx %in% rownames(expMatrix)]
			expProfile <- t(expMatrix[selTx,,drop=FALSE])
		}
		if (!is.null(expProfile)) {
			expProfile.range <- range(expProfile, na.rm=TRUE)
		} else {
			expProfile.range <- NULL
		}
		
		## includeExon1
		ord <- 1:ncol(methyGenoSet)
		if (!is.null(expMatrix.exon) && !is.null(tx2exon)) {
			selExon1 <- tx2exon[selTx]
			selExon1 <- selExon1[selExon1 %in% rownames(expMatrix.exon)]
			expProfile.exon1 <- t(expMatrix.exon[selExon1,,drop=FALSE])
			if (renameExon) {
				# colnames(expProfile.exon1) <- paste(selTx, 'exon1', sep='_')
				# sigExon1 <- paste(sigTx, 'exon1', sep='_')
				colnames(expProfile.exon1) <- paste('exon1', selTx, sep='_')
				sigExon1 <- paste('exon1', sigTx, sep='_')
			} else {
				sigExon1 <- tx2exon[sigTx]
			}
			expProfile <- cbind(expProfile, expProfile.exon1)

			selCol <- c(sigExon1, sigTx)
			selCol <- selCol[selCol %in% colnames(expProfile)]
			## sort only based on significant exon1 and transcript
			ord <- order(rowMeans(expProfile[, selCol, drop=FALSE], na.rm=TRUE), decreasing=FALSE)
		} else if (!is.null(expProfile)) {
			## sort only based on significant Tx
			ord <- order(rowMeans(expProfile[,sigTx,drop=FALSE], na.rm=TRUE), decreasing=FALSE)
		}
		## combine phenoData with expProfile if both exist
		if (!is.null(phenoData)) {
			## combine the phenoData with the expProfile matrix
			if (!is.null(expProfile)) {
				expProfile <- cbind(expProfile, as.matrix(phenoData))
			} else {
				expProfile <- as.matrix(phenoData)
			}
		}

		## ------------------------------------
		## plot the heatmap of gene.i
		if (!is.null(title.suffix)) {
			title <- paste(symbol, ' (', title.suffix, ')', sep='')
		} else {
			title <- paste(symbol, ' (GeneID:', gene.i, ')', sep='')
		}
		cat("Ploting ", title, '\n')
		
		if (is.null(main)) main <- title
		
		## plotting legend
		if (newPlot) grid.newpage()
		if (addLegend) {
			legendWidth <- 0.12 
			sepWidth <- 0.08
			plotWidth <- 1 - legendWidth - sepWidth
			layout.width <- c(plotWidth, sepWidth, legendWidth)
			pushViewport(viewport(layout=grid.layout(1, 3, widths=layout.width)))
			pushViewport(viewport(layout.pos.col=1, layout.pos.row=1))
		} 

		if (sortBy == 'expression' && !is.null(expProfile)) {
			plotInfo <- heatmapByChromosome(methyGenoSet[, ord,drop=F],	gene.i,	annotationTracks=annotationTracks, phenoData=expProfile[ord,,drop=F], 
																 phonoColorMap=phenoColor, sortSample=F, dataTrackName='Methylation Profile', main=main, 
																 cex.main=1, ylim=ylim, newPlot=FALSE, gradient=gradient, ncolor=ncolor, includeGeneBody=includeGeneBody, ...)
		} else {
			sortSample <- ifelse(sortBy == 'methylation', TRUE, FALSE)
			plotInfo <- heatmapByChromosome(methyGenoSet,	gene.i,	annotationTracks=annotationTracks, phenoData=expProfile, 
																 phonoColorMap=phenoColor, sortSample=sortSample, dataTrackName='Methylation Profile', main=main, 
																 cex.main=1, ylim=ylim, newPlot=FALSE, gradient=gradient, ncolor=ncolor, includeGeneBody=includeGeneBody, ...)
		}
		
		## plot legendInfo
		if (addLegend) {
			popViewport(1)
			## plot legend information
			pushViewport(viewport(layout.pos.col=3, layout.pos.row=1))
			
			## determine the height of legends
			## the height of methylation and expression colorbars are 2*legendHeight
			numOfOtherLegend <- length(which(names(phenoColor) != 'gradient'))
			if (!is.null(expProfile)) {
				legendHeight <- (1 -	(3 + numOfOtherLegend) * 0.05) / (4 + numOfOtherLegend)
			} else {
				legendHeight <- (1 -	(2 + numOfOtherLegend) * 0.05) / (2 + numOfOtherLegend)
			}
			if (legendHeight > 1/8) legendHeight <- 1/8
			
			x0 <- 0.1 # 0.15 ## starting plotting point in x-axis
			colWidth <- 0.15
			
			## plot methylation legend
			stepHeight <- legendHeight * 2 * 0.9 / ncolor
			ystart <- 1 - (0.05 + legendHeight * 2 )
			methyColor <- gradient

			grid.rect(x0, stepHeight * (1:ncolor) + ystart, width=colWidth, height=stepHeight, 
								gp=gpar(col=methyColor, fill=methyColor), default.units="npc", just=c("left", "top"))
			# add tick information
			if (useBetaValue) {
				ytickLabel <- round(seq(0, 1, length=5), 2)
			} else {
				ytickLabel <- round(seq(ylim[1], ylim[2], length=5), 2)
			}
			ytickPos <- ystart + legendHeight * 2 * 0.9 * seq(0,1,length=5)
			grid.segments(x0 + colWidth, ytickPos, x0 + colWidth + 0.05, ytickPos, default.units = "npc")
			## add tick labels
			fontsize <- round(as.numeric(convertX(unit(0.8, 'npc'), 'points'))/6)
			grid.text(ytickLabel, x=x0 + colWidth + 0.08, y=ytickPos, just=c("left", "center"), default.units = "npc", gp=gpar(fontsize=fontsize))
			## add legend title
			if (is.null(methylationLegendTitle)) {
				methylationLegendTitle <- c('Methylation', ifelse(useBetaValue, '(Beta value)', '(M value)'))
			} else {
				methylationLegendTitle <- strsplit(methylationLegendTitle, '\n')[[1]]
			}
			grid.text(methylationLegendTitle[1], x=0.5, y=max(ytickPos) + 0.025 + as.numeric(convertY(unit(fontsize, 'points'), 'npc')), just=c('center', 'bottom'), 
						default.units = "npc", gp=gpar(fontsize=fontsize, fontface='bold'))
			if (length(methylationLegendTitle) > 1)
				grid.text(methylationLegendTitle[2], x=0.5, y=max(ytickPos) + 0.02, just=c('center', 'bottom'), default.units = "npc", gp=gpar(fontsize=fontsize, fontface='bold'))

			## plot expression legend
			if (!is.null(expProfile.range)) {
				if ('gradient' %in% names(phenoColor)) {
					gradient.exp <- phenoColor$gradient
				} else {
					gradient.exp <- gradient
				}
				ystart <- 1 - (0.1 + legendHeight * 4 )
				expColor <- colorRampPalette(gradient.exp)(ncolor)[1:ncolor]
				ytickLabel <- round(seq(expProfile.range[1], expProfile.range[2], length=5), 2)
				ytickPos <- ystart + legendHeight * 2 * 0.9 * seq(0,1,length=5)

				grid.rect(x0, stepHeight * (1:ncolor) + ystart, width=colWidth, height=stepHeight, 
									gp=gpar(col=expColor, fill=expColor), default.units="npc", just=c("left", "top"))
				# add tick information
				grid.segments(x0 + colWidth, ytickPos, x0 + colWidth + 0.05, ytickPos, default.units = "npc")
				## add tick labels
				grid.text(ytickLabel, x=x0 + colWidth + 0.08, y=ytickPos, just=c("left", "center"), default.units = "npc", gp=gpar(fontsize=fontsize))
				## add title
				if (is.null(expressionLegendTitle)) {
					expressionLegendTitle <- c('Expression', 'log2-RPKM')
				} else {
					expressionLegendTitle <- strsplit(expressionLegendTitle, '\n')[[1]]
				}
				grid.text(expressionLegendTitle[1], x=0.5, y=max(ytickPos) + 0.025 + as.numeric(convertY(unit(fontsize, 'points'), 'npc')), just=c('center', 'bottom'), 
							default.units = "npc", gp=gpar(fontsize=fontsize, fontface='bold'))
				if (length(expressionLegendTitle) > 1)			
					grid.text(expressionLegendTitle[2], x=0.5, y=max(ytickPos) + 0.02, just=c('center', 'bottom'), default.units = "npc", gp=gpar(fontsize=fontsize, fontface='bold'))
			}

			## Add other phenoColor legends
			if (numOfOtherLegend > 0) {
		
				phenoColor <- phenoColor[names(phenoColor) != 'gradient']
				## calculate the stepHeight based on font size,which is defined by the number of color levels
				totalLevel <- length(unlist(phenotypeLevels))
				if (!is.null(expProfile.range)) {
					ystart <- 1 - (0.1 + legendHeight * 4 )
				} else {
					ystart <- 1 - (0.1 + legendHeight * 2 )
				}
				stepHeight <- (ystart - 0.1 * length(phenoColor)) / totalLevel
				stepHeight <- min(stepHeight, as.numeric(convertY(unit(fontsize, 'points'), 'npc')) * 1.5)
				# stepWidth <- convertX(convertY(unit(stepHeight, 'npc'), 'points'), 'npc')
				rectWidth <- 0.08
				rectHeight <- min(stepHeight, as.numeric(convertY(convertX(unit(0.05, 'npc'), 'points'), 'npc')))
				fontsize.color <- floor(as.numeric(convertY(unit(stepHeight * 0.9, 'npc'), 'points')))
				fontsize.color <- min(fontsize.color, fontsize)
				for (i in 1:length(phenoColor)) {
					phenoColor.i <- phenoColor[[i]]
					phenoTypeName.i <- names(phenoColor)[i]
					phenotypeLevels.i <- phenotypeLevels[[phenoTypeName.i]]
					
					## add title
					grid.text(phenoTypeName.i, x=0.5, y=ystart - 0.05, just=c('center', 'top'), default.units = "npc", gp=gpar(fontsize=fontsize, fontface='bold'))

					for (j in 1:length(phenoColor.i)) {
						y0 <- ystart - 0.05 - stepHeight * (j + 3/4) 
						grid.rect(x0, y0, width=rectWidth, height=rectHeight, 
											gp=gpar(col=phenoColor.i[j], fill=phenoColor.i[j]), default.units="npc", just=c("left", "center"))
						## add tick labels
						grid.text(phenotypeLevels.i[j], x=x0 + rectWidth * 2, y=y0, just=c("left", "center"), default.units = "npc", gp=gpar(fontsize=fontsize.color))
					}
					## update ystart
					ystart <- y0
				}
			}
			
			## plot rectangle around legend
			# grid.rect(0, 1, width=1, height=min(1, abs(1 - ystart + 0.03)), gp=gpar(col=1, lty=1, lwd=1, fill=rgb(1,1,1, alpha=0)), default.units="npc", just=c("left", "top"))
			
			plotInfo <- c(plotInfo, layout.width=layout.width)
			popViewport(1)
		}			
		
	}) 
	
	return(invisible(plotResult))		
}




plotHeatmapByGene <- function(selGene, genoSet, phenoData=NULL, sortBy=c(NA, 'phenoData', 'data'), includeGeneBody=FALSE,  sortByTx=FALSE,
	CpGInfo=NULL, genomicFeature=NULL, phenoColor=list(gradient=c("green", "black", "red")), title.suffix=NULL, addLegend=TRUE, 
	genoSetLegendTitle=NULL, gradient=c("blue", "white", "red"), ncolor=16, main=NULL, newPlot=TRUE, ylim=NULL, ...) {
		
	sortBy <- as.character(sortBy)
	sortBy <- match.arg(sortBy)
	if(is.na(sortBy)) sortBy <- 'NA'
	if (length(selGene) > 1) {
		warnings('Only the first gene will be plotted!')
		selGene <- selGene[1]
	}
	
	## convert gradient as a vector of colors if not
	if (!all(grepl("^#", gradient)) || length(gradient) < 5) {
		gradient <- colorRampPalette(gradient)(ncolor)
	} 
	ncolor <- length(gradient)

	if (!is.list(genoSet)) {
		genoSetList <- list(genoSet)
	} else {
		genoSetList <- genoSet
	}
	if (!all(sapply(genoSetList, function(x) is(x, 'GenoSet')))) 
		stop('"genoSet" must be a GenoSet object or a list of GenoSet objects.')
	if (is.null(ylim)) {
		ylim <- range(unlist(lapply(genoSetList, function(x) range(assayData(x)$exprs))))
	}
	
	phenotypeLevels <- NULL
	if (!is.null(phenoData)) {
		if (is(phenoData, "ExpressionSet")) {
			phenoData <- exprs(phenoData)
		} 
		rn <- rownames(phenoData)
		if (is.data.frame(phenoData) || is.list(phenoData)) {
			phenotypeLevels <- lapply(phenoData, function(x) {
				x <- levels(as.factor(x))
				return(x)
				})
			phenoData <- data.frame(lapply(phenoData, function(x) {
				x <- round(as.numeric(as.factor(x)))
				if (min(x, na.rm=TRUE) <= 0) x <- x - min(x, na.rm=TRUE) + 1
				return(x)
				}), check.names = FALSE)
		}
		rownames(phenoData) <- rn

		otherPhenoName <- colnames(phenoData)
		otherPhenoName <- otherPhenoName[!(otherPhenoName %in% names(phenoColor))]
		if (length(otherPhenoName) > 0) {
			allPhenoName <- names(phenoColor)
			for (phenoName.i in otherPhenoName) {
				maxLevel.i <- max(phenoData[[phenoName.i]], na.rm=TRUE)
				if (maxLevel.i <= 6) {
					phenoColor <- c(phenoColor, list(1:maxLevel.i))
					allPhenoName <- c(allPhenoName, phenoName.i)
				}
			}
			names(phenoColor) <- allPhenoName
		}		
	}
	
	annotationTracks <- buildAnnotationTracks(gene=selGene, includeGeneBody=includeGeneBody, CpGInfo=CpGInfo, genomicFeature=genomicFeature, ...)

	## sort the transcript based on annotation track
	geneRegionTrack <- annotationTracks[sapply(annotationTracks, class) == 'GeneRegionTrack'][[1]]
	## estimate the order of transcripts in the geneRegionTrack
	grange2show <- attr(annotationTracks, 'grange2show')
	if (is.null(grange2show)) {
		warnings('No region to plot!')
		return(NULL)
	}

	grange2show <- checkChrName(grange2show, addChr=TRUE)
	chromosome <- as.character(seqnames(grange2show)[1])
	if (sortByTx) {
		geneRegionTrack <- .estimateStackLocation(geneRegionTrack, from=start(grange2show)[1], to=end(grange2show)[1], chromosome=chromosome)
		annTx <- split(values(geneRegionTrack)$transcript, stacks(geneRegionTrack))
		annTx <- rev(unique(as.character(unlist(annTx))))
		
		genoSetList <- lapply(genoSetList, function(x) {
			x.name <- sampleNames(x)
			ind <- 1:length(x.name)
			names(ind) <- x.name
			x.name.ord <- c(annTx[annTx %in% x.name], x.name[!(x.name %in% annTx)])
			x <- x[,ind[x.name.ord]]
		})
	}

	## ------------------------------------
	## plot the heatmap of selGene
	title <- NULL
	if (is.character(selGene)) {
		symbol <- unlist(lookUp(selGene, 'org.Hs.eg.db', 'SYMBOL'))
		if (!is.null(title.suffix)) {
			title <- paste(symbol, ' (', title.suffix, ')', sep='')
		} else {
			title <- paste(symbol, ' (GeneID:', selGene, ')', sep='')
		}
	} else if (is(selGene, 'GRanges')) {
		title <- paste(seqnames(selGene)[1], start(selGene)[1], end(selGene)[1], sep='_')
	}

	cat("Ploting ", title, '\n')
	
	if (is.null(main)) main <- title
	
	## plotting legend
	if (newPlot) grid.newpage()
	if (addLegend) {
		legendWidth <- 0.12 
		sepWidth <- 0.08
		plotWidth <- 1 - legendWidth - sepWidth
		layout.width <- c(plotWidth, sepWidth, legendWidth)
		pushViewport(viewport(layout=grid.layout(1, 3, widths=layout.width)))
		pushViewport(viewport(layout.pos.col=1, layout.pos.row=1))
	} 
	
	sortSample <- ifelse(sortBy == 'methylation', TRUE, FALSE)
	plotInfo <- heatmapByChromosome(genoSetList, selGene, annotationTracks=annotationTracks, phenoData=phenoData, 
				 phonoColorMap=phenoColor, sortSample=sortSample, dataTrackName='Methylation Profile', main=main, ylim=ylim,
				 cex.main=1, newPlot=FALSE, gradient=gradient, ncolor=ncolor, includeGeneBody=includeGeneBody, ...)
	
	## plot legendInfo
	if (addLegend) {
		
		popViewport(1)
		## plot legend information
		pushViewport(viewport(layout.pos.col=3, layout.pos.row=1))
		
		## determine the height of legends
		## the height of genoSet data 
		numOfOtherLegend <- length(which(names(phenoColor) != 'gradient'))
		legendHeight <- (1 - (1 + numOfOtherLegend) * 0.05) / (1 + numOfOtherLegend)
		if (legendHeight > 1/8) legendHeight <- 1/8
		x0 <- 0.1 # 0.15 ## starting plotting point in x-axis
		colWidth <- 0.15
		
		## plot genoset data legend
		stepHeight <- legendHeight * 2 * 0.9 / ncolor
		ystart <- 1 - (0.05 + legendHeight * 2 )
		methyColor <- gradient[1:ncolor]

		grid.rect(x0, stepHeight * (1:ncolor) + ystart, width=colWidth, height=stepHeight, 
							gp=gpar(col=methyColor, fill=methyColor), default.units="npc", just=c("left", "top"))
							
		# add tick information
		ytickLabel <- round(seq(ylim[1], ylim[2], length=5), 2)
		ytickPos <- ystart + legendHeight * 2 * 0.9 * seq(0,1,length=5)
		grid.segments(x0 + colWidth, ytickPos, x0 + colWidth + 0.05, ytickPos, default.units = "npc")
		## add tick labels
		fontsize <- round(as.numeric(convertX(unit(0.8, 'npc'), 'points'))/6)
		grid.text(ytickLabel, x=x0 + colWidth + 0.08, y=ytickPos, just=c("left", "center"), default.units = "npc", gp=gpar(fontsize=fontsize))
		## add legend title
		if (!is.null(genoSetLegendTitle)) {
			genoSetLegendTitle <- strsplit(genoSetLegendTitle, '\n')[[1]]
		}
		grid.text(genoSetLegendTitle[1], x=0.5, y=max(ytickPos) + 0.025 + as.numeric(convertY(unit(fontsize, 'points'), 'npc')), just=c('center', 'bottom'), 
					default.units = "npc", gp=gpar(fontsize=fontsize, fontface='bold'))
		if (length(genoSetLegendTitle) > 1)
			grid.text(genoSetLegendTitle[2], x=0.5, y=max(ytickPos) + 0.02, just=c('center', 'bottom'), default.units = "npc", gp=gpar(fontsize=fontsize, fontface='bold'))

		## Add other phenoColor legends
		if (numOfOtherLegend > 0) {
	
			phenoColor <- phenoColor[names(phenoColor) != 'gradient']
			## calculate the stepHeight based on font size,which is defined by the number of color levels
			totalLevel <- length(unlist(phenotypeLevels))
			ystart <- 1 - (0.1 + legendHeight * 2 )
			stepHeight <- (ystart - 0.1 * length(phenoColor)) / totalLevel
			stepHeight <- min(stepHeight, as.numeric(convertY(unit(fontsize, 'points'), 'npc')) * 1.5)
			# stepWidth <- convertX(convertY(unit(stepHeight, 'npc'), 'points'), 'npc')
			rectWidth <- 0.08
			rectHeight <- min(stepHeight, as.numeric(convertY(convertX(unit(0.05, 'npc'), 'points'), 'npc')))
			fontsize.color <- floor(as.numeric(convertY(unit(stepHeight * 0.9, 'npc'), 'points')))
			fontsize.color <- min(fontsize.color, fontsize)
			for (i in 1:length(phenoColor)) {
				phenoColor.i <- phenoColor[[i]]
				phenoTypeName.i <- names(phenoColor)[i]
				phenotypeLevels.i <- phenotypeLevels[[phenoTypeName.i]]
				
				## add title
				grid.text(phenoTypeName.i, x=0.5, y=ystart - 0.05, just=c('center', 'top'), default.units = "npc", gp=gpar(fontsize=fontsize, fontface='bold'))

				for (j in 1:length(phenoColor.i)) {
					y0 <- ystart - 0.05 - stepHeight * (j + 3/4) 
					grid.rect(x0, y0, width=rectWidth, height=rectHeight, 
										gp=gpar(col=phenoColor.i[j], fill=phenoColor.i[j]), default.units="npc", just=c("left", "center"))
					## add tick labels
					grid.text(phenotypeLevels.i[j], x=x0 + rectWidth * 2, y=y0, just=c("left", "center"), default.units = "npc", gp=gpar(fontsize=fontsize.color))
				}
				## update ystart
				ystart <- y0
			}
		}
		
		## plot rectangle around legend
		# grid.rect(0, 1, width=1, height=min(1, abs(1 - ystart + 0.03)), gp=gpar(col=1, lty=1, lwd=1, fill=rgb(1,1,1, alpha=0)), default.units="npc", just=c("left", "top"))

		plotInfo <- c(plotInfo, layout.width=layout.width)
		popViewport(1)
	}			

	return(invisible(plotInfo))		
}


## 
## plot the heatmap data tracks with sample (row) information
plotTracksWithDataTrackInfo <- function(trackList, labels=NULL, grange2show=NULL, dataTrackName=NULL, dataInfo=NULL, dataColorMap=NULL, 
			dataInfoRange=NULL, dataBackground=gray(0.9), minHeatmapColumnWidth=2, labelWidth=0.1, gradient=c("blue", "white", "red"), 
			ncolor=16, main='', newPlot=FALSE, sizes=NULL, ...) {
	
	if (missing(trackList)) {
		stop('Please provide "trackList"!')
	}
	
	## convert gradient as a vector of colors if not
	if (!all(grepl("^#", gradient)) || length(gradient) < 5) {
		gradient <- colorRampPalette(gradient)(ncolor)
	} 
	ncolor <- length(gradient)
	
	if (!is(trackList, 'list')) trackList <- list(trackList)
	dataTrackInd <- which(sapply(trackList, function(x) class(x) == 'DataTrack' && getPar(x, 'type') == 'heatmap'))
	
	if (length(dataTrackInd) == 0) {
		warning('No DataTrack was found!')
		dataTrackName <- NULL
	} else {
		if (is.null(dataTrackName)){
			dataTrackName <- sapply(trackList[dataTrackInd], names)
		}
		if (is.null(labels)) {
			labels <- lapply(trackList[dataTrackInd], function(x) rownames(values(x)))
			names(labels) <- dataTrackName
		}	else {
			if (is.list(labels)) {
				if (is.null(names(labels))) {
					stop('Please provide the list names of the labels! The names should be consistent with the dataTrack names.')
				} else {
					dataTrackName <- intersect(dataTrackName, names(labels))
					if (length(dataTrackName) < length(labels)) {
						warnings('Some label names are inconsistent with dataTrack names!')
					}
				}
			}	else {
				labels <- list(labels)
				names(labels) <- dataTrackName[1]
			}
		}
		
	} 

	## start plotting 
	if (newPlot)	grid.newpage()
	
	## determine the layout based on expression profile
	# labelWidth <- 0.1
	if (!is.null(dataInfo)) {
		## 
		if (is.data.frame(dataInfo)) {
			dataInfo <- as.matrix(dataInfo)
		} else if (!is.matrix(dataInfo)) {
			dataInfo <- matrix(dataInfo, ncol=1)
		}
		if (is.null(rownames(dataInfo))) {
			stop('Please provide labels to dataInfo, which should match the labels!') 
		}
		nn <- ncol(dataInfo)
		labelWidth <- labelWidth + 0.02 * nn
		if (labelWidth > 0.3) labelWidth <- 0.3
	} 

	legendWidth <- 0
	layout.col.num <- 2
	plotWidth <- 1 - labelWidth
	pushViewport(viewport(layout=grid.layout(1, layout.col.num, widths=c(plotWidth, labelWidth))))
	 
	chromosome <- as.character(seqnames(grange2show)[1])
	pushViewport(viewport(layout.pos.col=1, layout.pos.row=1))
	
	## estimate the space of tracks
	trackHeights <- .estimateTrackHeight(trackList, grange2show, sizes=sizes)
	
	## set the minimum width of the heatmap columns to have better visualization effects
	if (minHeatmapColumnWidth > 1) {
		newWidth <- ceiling(minHeatmapColumnWidth / as.numeric(convertX(unit(0.95, 'npc'), 'points')) * width(grange2show))
		trackList[dataTrackInd] <- lapply(trackList[dataTrackInd], function(x) {
			ww.points <- floor(width(x)/width(grange2show) * as.numeric(convertX(unit(0.95, 'npc'), 'points')))
			width(x)[ww.points < minHeatmapColumnWidth] <- newWidth
			return(x)
		})
	}
	
	## set background color if specified
	if (!is.na(dataBackground) && !is.null(dataBackground)) {
		for (dataTrackInd.i in dataTrackInd) {
			displayPars(trackList[[dataTrackInd.i]]) <- list(background.panel=dataBackground)
		}
	}
	
	## make sure the transcript names are shown in right side of plot
	geneRegionTrackInd <- which(sapply(trackList, class) == 'GeneRegionTrack') 
	if (length(geneRegionTrackInd) > 0) {
		selTrack <- trackList[[geneRegionTrackInd]]
		trackLabelWidth <- max(nchar(transcript(selTrack))) * getPar(selTrack, 'fontsize') * 3/5 * getPar(selTrack, 'cex')
		trackLabelWidthRatio <- trackLabelWidth/as.numeric(convertX(unit(0.9, 'npc'), 'points'))
		trackLabelWidth.points <- width(grange2show) * trackLabelWidthRatio - 2000
		if (trackLabelWidth.points < 0) trackLabelWidth.points <- 0
		start(grange2show) <- start(grange2show) - trackLabelWidth.points
	}
	ss <- start(grange2show)[1]
		
	plotInfo <- plotTracks(trackList, from=ss, to=end(grange2show)[1], sizes=trackHeights,
			 chromosome=chromosome,	add=TRUE, main=main, ...)
	## retrieve the plot coordinates		
	plotLoc <- coords(plotInfo$title)

	## plot a rectangle around the heatmap if necessary
	# if (!is.na(dataBackground) && !is.null(dataBackground)) {
	# 	totalHeight <- floor(as.numeric(convertY(unit(1, 'npc'), 'points')))
	# 	totalWidth <- floor(as.numeric(convertX(unit(1, 'npc'), 'points')))
	# 	margin <- 6
	# 	for (dataTrackName.i in dataTrackName) {
	# 		allHeight <- plotLoc[,'y2'] - plotLoc[,'y1']
	# 		y0 <- (plotLoc[dataTrackName.i, 'y1'] + margin) / totalHeight
	# 		height <- (allHeight[dataTrackName.i] - 2 * margin) / totalHeight
	# 		x0 <- plotLoc[dataTrackName.i, 'x2']/totalWidth 
	# 		width <- 1 - x0 - margin/totalWidth
	# 		grid.rect(x0, y=1-y0, width=width, height=height, just=c('left', 'top'), gp=gpar(lty=1, col='dark gray', fill=dataBackground, alpha=0.15))
	# 	}
	# }

	popViewport(1)
	## plot annotation or phenotype information
	pushViewport(viewport(layout.pos.col=2, layout.pos.row=1))

	defaultPar <- Gviz:::.parMappings$GdObject		

	allHeight <- plotLoc[,'y2'] - plotLoc[,'y1']
	margin <- plotLoc[1,'x1']
	y00 <- as.numeric(convertY(unit(1, 'npc'), 'points')) - sum(allHeight) - margin
	allHeight <- allHeight/as.numeric(convertY(unit(1, 'npc'), 'points'))
	## to handle more than one data tracks
	for (i in 1:length(dataTrackName)) {
		dataTrackName.i <- dataTrackName[i]
		labels.i <- labels[[dataTrackName.i]]

		## calculate the label positions 
		y0 <- (plotLoc[dataTrackName.i, 'y1'] - min(plotLoc[, 'y1']) + y00) / as.numeric(convertY(unit(1, 'npc'), 'points'))
		# 
		hh <- allHeight[dataTrackName.i]
		realHeight <- hh / 1.1
		# grid.rect(0, y=1 - y0 - hh/2, width=1, height=realHeight, just=c('left', 'center'), gp=gpar(lty='dashed', col=3))

		ystart <- 1 - y0 - hh/2 - realHeight / 2
		yend <- 1 - y0 - hh/2 + realHeight / 2
		stepHeight <- realHeight / length(labels.i)
		x0 <- 0.05
		colWidth <- 0.1
		if (!is.null(dataInfo)) {
			# plot the colnames first
			colWidth <- 1 / (5 + ncol(dataInfo))
			if (i == 1 && !is.null(colnames(dataInfo))) {
				fontsize <- floor(as.numeric(convertX(unit(colWidth * 0.95, 'npc'), 'points')))
				fontsize <- min(fontsize, floor(as.numeric(convertY(unit((1-yend-0.01)/6, 'npc'), 'points')))) # based on height
				grid.text(colnames(dataInfo), x=x0 + colWidth * (1:ncol(dataInfo)), y = yend + 0.01, rot=90, just=c('left', 'bottom'), 
					gp=gpar(fontfamily=defaultPar$fontfamily, fontsize=fontsize, fontface=defaultPar$fontface, col=1))
			}

			## check the consistency of dataInfo rownames and labels.i
			if (is.null(rownames(dataInfo))) {
				if (nrow(dataInfo) == length(labels.i)) {
					rownames(dataInfo) <- labels.i
				} else {
					stop('Please provide rownames to dataInfo, which should be consistent to the labels!')
				}
			} else {
				if (!all(labels.i %in% rownames(dataInfo))) {
					labels.i <- make.names(labels.i)
					rownames(dataInfo) <- make.names(rownames(dataInfo))
					if (!all(labels.i %in% rownames(dataInfo))) {
						stop("'labels' does not match rownames of dataInfo!\n")
					}
				}
			}

			dataColor <- vector(mode='list', length=ncol(dataInfo))
			names(dataColor) <- colnames(dataInfo)
			gradient.data <- gradient
			if (!is.null(dataColorMap)) {
				if (!is.list(dataColorMap)) stop('dataColorMap should be a named list!')

				if ('gradient' %in% names(dataColorMap)) {
					gradient.data <- dataColorMap$gradient
					if (!all(grepl("^#", gradient.data)) || length(gradient.data) < 5) 
						gradient.data <- colorRampPalette(gradient.data)(ncolor)
				}
				## convert colors to the format like "#FF0000"
				dataCols <- colnames(dataInfo)[(colnames(dataInfo) %in% names(dataColorMap))]
				otherCol <- colnames(dataInfo)[!(colnames(dataInfo) %in% names(dataColorMap))]
				for (dataCol.i in dataCols) {
					colorMap.i <- dataColorMap[[dataCol.i]]
					if (is.numeric(colorMap.i)) {
						paletteColor <- palette()
						colorMap.i <- rgb(t(col2rgb(paletteColor[round(colorMap.i) %% length(paletteColor)]))/255)
					} else if (length(grep("^#", colorMap.i)) < length(colorMap.i)) {
						colorMap.i <- rgb(t(col2rgb(colorMap.i))/255)
					}
					data.i <- dataInfo[labels.i,dataCol.i]
					## if there are larger than 10 color levels, the data will be scaled to the data range
					if (length(colorMap.i) > 10) {
						## dataInfoRange is to control the plot color range
						if (!is.null(dataInfoRange)) {
							if (is.list(dataInfoRange)) {
								dataInfoRange.i <- dataInfoRange[[dataCol.i]]
								if (is.null(dataInfoRange.i)) dataInfoRange.i <- dataInfoRange[[1]]
							} else {
								dataInfoRange.i <- dataInfoRange
							}
							data.i <- Gviz:::.z2icol(data.i, length(colorMap.i), xrange=dataInfoRange.i)
						} else {
							data.i <- Gviz:::.z2icol(data.i, length(colorMap.i), xrange=range(data.i))	
						}
					} 
					dataColor[[dataCol.i]] <- colorMap.i[round(data.i)]
				}

				if (length(otherCol) > 0) {
					for (i in 1:length(otherCol)) {
						valsScaled <- Gviz:::.z2icol(dataInfo[labels.i,otherCol, drop=FALSE], length(gradient.data), xrange=range(dataInfo[,otherCol,drop=FALSE]))
						dataColor[[otherCol[i]]] <- gradient.data[valsScaled[,i]]
					}
				}
			} else {
				## dataInfoRange is to control the plot color range
				if (!is.null(dataInfoRange)) {
					valsScaled <- Gviz:::.z2icol(dataInfo[labels.i,,drop=FALSE], ncolor, xrange=dataInfoRange)
				} else {
					valsScaled <- Gviz:::.z2icol(dataInfo[labels.i,,drop=FALSE], ncolor, xrange=range(dataInfo))	
				}
				for (i in 1:ncol(dataInfo)) {
					dataColor[[i]] <- gradient[valsScaled[,i]]
				}
			}

			## plot data matrix as a regular heatmap
			for (i in 1:ncol(dataInfo)) {
				dataCol.i <- colnames(dataInfo)[i]
				grid.rect(x0, stepHeight * (1:length(labels.i)) + ystart, width=colWidth, height=stepHeight, 
									gp=gpar(col=dataColor[[dataCol.i]], fill=dataColor[[dataCol.i]]),
									default.units="npc", just=c("left", "top"))
				x0 <- x0 + colWidth
			}
			x0 <- x0 + min(0.05, colWidth/2)
		}
		fontsize <- floor(as.numeric(convertY(unit(stepHeight * 0.9, 'npc'), 'points')))
		if (!is.null(dataInfo)) {
			realLabelWidth <- 1 - (colWidth) * ncol(dataInfo) - min(0.1, colWidth)
		} else {
			realLabelWidth <- 0.95
		}
		numchar <- min(max(nchar(labels.i)), 5)
		fontsizeW <- round(as.numeric(convertX(unit(realLabelWidth, 'npc'), 'points'))/ numchar)
		fontsize <- min(fontsize, fontsizeW)
		grid.text(labels.i, x=x0, y = stepHeight * (1:length(labels.i)) + ystart - stepHeight/2, just=c('left', 'center'), gp=gpar(fontfamily=defaultPar$fontfamily, fontsize=fontsize, fontface=defaultPar$fontface, col=1))
	}
	popViewport(2)
	
	plotInfo <- c(plotInfo, labelWidth=labelWidth)
	attr(plotInfo, 'grange2show') <- grange2show
	return(invisible(plotInfo))
}


## estimate the Gviz track heights
# 
# grange2show: a GRanges object specify the start and end of the plot range
.estimateTrackHeight <- function(trackList, grange2show, sizes=NULL, minPoints=50) {
	
	trackList <- lapply(trackList, Gviz::subset, from=start(grange2show), to=end(grange2show), chromosome=seqnames(grange2show))
  trackList <- lapply(trackList, Gviz:::setStacks, from=start(grange2show), to=end(grange2show))
	spaceSetup <- Gviz:::.setupTextSize(trackList, sizes=sizes)
	totalHeight <- as.numeric(convertY(unit(1, 'npc'), 'points'))
	space.points <- totalHeight * spaceSetup$spaceNeeded
	smallTrackInd <- which(space.points < minPoints)
	if (length(smallTrackInd) > 0) {
		space.points[smallTrackInd] <- minPoints
		remainPoints <- totalHeight - minPoints * length(smallTrackInd)
		space.points[-smallTrackInd] <- remainPoints * space.points[-smallTrackInd] / sum(space.points[-smallTrackInd])
		space.points[space.points < minPoints] <- minPoints
	}
	sizes.new <- space.points/sum(space.points)
	return(sizes.new)
}




## add chr to the chromosome names
checkChrName <- function(grange, addChr=TRUE) {
	if (is(grange, 'GRanges')) {
		chrName <- seqlevels(grange)
	} else if (is(grange, 'character')) {
		chrName <- grange
	} else if (is(grange, 'GenoSet')) {
		chrName <- chrNames(grange)
	} else if (is(grange, 'RangedData')) {
		chrName <- levels(space(grange))
	} else if (is(grange, 'RangeTrack')) {
		chrName <- seqlevels(ranges(grange))
	} else {
		chrName <- names(grange)
		if (is.null(chrName)) 
			stop('Un-supported data types!')
	}
	
	if (any(grepl('^chr[0-9XY][0-9]?', chrName))) {
		if (!addChr) {
			chrName <- sub('^chr', '', chrName)
		}
	} else if (addChr) {
		ind <- grep('^[0-9XY][0-9]?', chrName)
		ind <- c(ind, grep('^MT', chrName))
		chrName[ind] <- paste('chr', chrName[ind], sep='')
	}
	
	if (is(grange, 'GRanges')) {
		seqlevels(grange) <- chrName
	} else if (is(grange, 'character')) {
		grange <- chrName
	} else if (is(grange, 'GenoSet')) {
		names(locData(grange)) <- chrName
	} else if (is(grange, 'RangeTrack')) {
		seqlevels(ranges(grange)) <- chrName
	} else {
		names(grange) <- chrName
	} 
	
	return(grange)
}


getCytoBand.ucsc <- function(hgVersion='hg19') {
	
	session <- browserSession()
	query <- ucscTableQuery(session, "cytoBandIdeo", hgVersion)
	detailedInfo <- getTable(query)
	ranges <- GRanges(seqnames=detailedInfo$chrom, ranges=IRanges(start=detailedInfo$chromStart, end=detailedInfo$chromEnd),
								name=detailedInfo$name, type=detailedInfo$gieStain)
	return(ranges)
}


getCpGIsland.ucsc <- function(hgVersion='hg19') {
	
	session <- browserSession()
	query <- ucscTableQuery(session, "cpgIslandExt", hgVersion)
	# islands <- track(query)
	detailedInfo <- getTable(query)

	ranges <- GRanges(seqnames=detailedInfo$chrom, ranges=IRanges(start=detailedInfo$chromStart, end=detailedInfo$chromEnd),
									name=detailedInfo$name)
	return(ranges)
}

