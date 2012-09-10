
# transTrack <- createTranscriptTrack(gene='55432', genomicFeature='TxDb.Hsapiens.UCSC.hg19.knownGene', lib='org.Hs.eg.db', extendRange=c(2000, 2000))
createTranscriptTrack <- function(gene, 
				genomicFeature='TxDb.Hsapiens.UCSC.hg19.knownGene', 
				lib='org.Hs.eg.db', 
				genome='hg19', 
				extendRange=c(2000, 2000), 
				includeGeneBody=TRUE, 
				background.title='brown', ...) {
	
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
	if (is(gene, 'GRanges')) {
		grange2show <- gene
		gene <- names(gene)
		chromosome <- as.character(seqnames(grange2show)[1])
	} else if (require(lib, character.only=TRUE)) {
		chromosome <- lookUp(gene, lib, 'CHR')[[1]]
		grange2show <- NULL
	}
	chromosome <- checkChrName(chromosome, addChr=TRUE)	

	## identify related transcripts
	if (is.character(genomicFeature)) {
		if(!require(genomicFeature, character.only = TRUE)) {
			warning(paste(genomicFeature, 'library is not installed!'))
			genomicFeature <- NULL
		} else {
			genomicFeature <- get(genomicFeature)
		}
	}

	## Create transcript track based on provided genomicFeatures
	selTrans <- transTrack <- NULL
	if (is(genomicFeature, 'TranscriptDb')) {
		## set the grange2show based on the TranscriptDb object
		if (is.null(grange2show)) {
			# refseqs <- lookUp(gene, lib, 'REFSEQ')[[1]]
			# refseqs <- grep( "^NM", refseqs, value = TRUE )
			selTrans <- exons(genomicFeature, vals=list(gene_id=gene), columns=c('tx_name', 'exon_id'))	
			if (length(selTrans) == 0) {
				if (grepl('^GeneID:', gene, ignore.case=TRUE)) {
					gene <- sub('^GeneID:', '', gene)
				} else {
					gene <- paste('GeneID:', gene, sep='')
				}
				selTrans <- exons(genomicFeature, vals=list(gene_id=gene), columns=c('tx_name', 'exon_id'))
			}
			if (includeGeneBody) {
				rr <- cbind(start(selTrans), end(selTrans))
				ss <- min(rr) - extendRange[1]
				ee <- max(rr) + extendRange[2]
			} else {
				## only retrieve the region surrounding TSS
				tmp <- transcripts(genomicFeature, vals=list(gene_id=gene))
				tss <- ifelse(as.character(strand(tmp)) == '-', end(tmp), start(tmp))
				ss <- min(tss) - extendRange[1]
				ee <- max(tss) + extendRange[2]
			}
			grange2show <- GRanges(seqnames=chromosome, strand='*', ranges=IRanges(start=ss, end=ee))
			geneSymbol <- lookUp(gene, lib, 'SYMBOL')[[1]]
		} else {
			# if (is.null(grange2show)) stop('Please provide either gene id or grange2show parameter!')
			allTrans <- exons(genomicFeature, columns=c('tx_name', 'exon_id'))	
			selTrans <- allTrans[!is.na(GenomicRanges::match(allTrans, grange2show))]
			geneSymbol <- NULL
		}
		## flatten the tx_name (one exon to multiple transcripts)
		tx_name <- IRanges::as.list(values(selTrans)$tx_name)
		tx_len <- sapply(tx_name, length)
		selTrans <- rep(selTrans, tx_len)
		values(selTrans)$tx_name <- unlist(tx_name)
		selTrans <- checkChrName(selTrans, addChr=TRUE)
		id <- values(selTrans)$exon_id
		tx_name <- values(selTrans)$tx_name
		if (is.null(geneSymbol)) geneSymbol <- tx_name
		transTrack <- GeneRegionTrack(range=ranges(selTrans), strand=strand(selTrans), genome=genome, chromosome=chromosome, 
					 feature=tx_name, gene=geneSymbol, exon=id, transcript=tx_name, name="Gene Model", showId=TRUE, background.title=background.title, ...)
		
	} else if (is.null(grange2show)) {
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
			symbol="name", transcript="name", strand="strand", showId=TRUE, name="UCSC Genes", background.title=background.title, ...)
	}	
	attr(transTrack, 'grange2show') <- grange2show
	return(transTrack)
}


## .estimateStackLocation(annotationTracks[[3]], from=ranges[1], to=ranges[2])
.estimateStackLocation <- function(track, from, to, chromosome=NULL, ...) {
	track <- Gviz:::consolidateTrack(track, chromosome=chromosome, ...)

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
	if (is(gene, 'GRanges')) {
		chromosome <- as.character(seqnames(gene))
	} else if (require(lib, character.only=TRUE)) {
		chromosome <- lookUp(gene, lib, 'CHR')[[1]]
	}
	chromosome <- checkChrName(chromosome, addChr=TRUE)
		
	## create annotation tracks
	# Ideogram (cytoband) track
	allTracks <- iTrack <- NULL
	if (is.null(cytobandInfo) && length(readLines(url('http://genome.ucsc.edu'))) > 0) {
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
	} else if (length(readLines(url('http://genome.ucsc.edu'))) > 0) {
		## get CpG-island info from UCSC if CpGInfo is null
		CpGInfo <- getCpGIsland.ucsc(hgVersion=genome)
	} 

	if (is(CpGInfo, 'GRanges')) {
		cpgTrack <- AnnotationTrack(CpGInfo, chromosome=chromosome, genome=genome, name='CpG-island')
		allTracks <- c(allTracks, list(cpgTrack))
	}

	## Create Transcript annotation track
	transTrack	<- createTranscriptTrack(gene, genomicFeature=genomicFeature, lib=lib, 
					extendRange=extendRange, includeGeneBody=includeGeneBody, genome=genome, ...)
	if (!is.null(transTrack)) {
		allTracks <- c(allTracks, list(transTrack))
		## update grange2show
		attr(allTracks, 'grange2show') <- attr(transTrack, 'grange2show')
	}
	
	return(allTracks)
}



## require(Gviz)
## library(GenomicRanges)
## library(rtracklayer) # import.bed
## library(GenomicFeatures) # TranscriptDb
### importFrom(GenomicFeatures, 'exons', )
### importFrom(annoate, lookUp)
## library(genoset) # GenoSet, locData

heatmapByChromosome <- function(
					methyGenoSet, 				# GenoSet object	
					gene, 								# Entrez gene ids or a GRanges object with length equals one
					annotationTracks=NULL,
					otherTrackList=NULL, 	# other Tracks objects supported by Gviz
					expProfile=NULL, 			# Expression profile
					# expProfileColors=NULL, # a phenotype data.frame (required for legend) or color vector
					expColorMap = NULL,	# a list of colormaps for every column-side rows
					extendRange=c(2000, 2000), # extended range on each side of the	gene
					includeGeneBody=TRUE, # wether to include genebody of the provided gene
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
					main='',							# title
					... ) {
	
	## check parameters
	if (length(gene) > 1) {
		warning('Current version can only support one gene per plot!')
		gene <- gene[1]
	}
	if (!is.null(expProfile)) {
		if (!is.matrix(expProfile)) {
			expProfile <- as.matrix(expProfile)
			if (is.null(colnames(expProfile)) && ncol(expProfile) == 1)
				colnames(expProfile) <- 'Expression Profile'
		}		
		if (nrow(expProfile) != ncol(methyGenoSet)) {
			if (ncol(expProfile) == ncol(methyGenoSet)) {
				expProfile <- t(expProfile)
			} else {
				stop('The dimension of "expProfile" does not match with the methyGenoSet!')
			}
		}
	}
	
	if (is.null(annotationTracks)) {
		annotationTracks <- buildAnnotationTracks(gene=gene, extendRange=extendRange, 
							includeGeneBody=includeGeneBody, cytobandInfo=cytobandInfo, CpGInfo=CpGInfo,
							genomeAxis=genomeAxis, lib=lib,	genome=genome, genomicFeature=genomicFeature)
		#
		if (is.null(annotationTracks))
			stop('"annotationTracks" is missing!')
	}
	allTracks <- annotationTracks
	## get grange2show
	grange2show <- attr(annotationTracks, 'grange2show')
	grange2show <- checkChrName(grange2show, addChr=TRUE)
	chromosome <- as.character(seqnames(grange2show))[1]
	
	## add otherTracks if provided
	if (!is.null(otherTrackList)) {
		if (is.list(otherTrackList)) {
			allTracks <- c(allTracks, otherTrackList)
		} else {
			allTracks <- c(allTracks, list(otherTrackList))
		}
	}
	
	## select related methylation data	
	methyGenoSet <- checkChrName(methyGenoSet, addChr=TRUE)
	grange.data <- suppressWarnings(as(locData(methyGenoSet), 'GRanges'))
	selMethyData <- methyGenoSet[!is.na(GenomicRanges::match(grange.data, grange2show)),]
	if (nrow(selMethyData) == 0) {
		warning("There is no methylation data exist in the selected grange2show!")
		return(NULL)
	}
	if ((sortSample || length(sortSample) == ncol(selMethyData)) && nrow(selMethyData) > 0) {
		hcr <- hclust(dist(t(assayData(selMethyData)$exprs)))
		ddr <- as.dendrogram(hcr)
		if (length(sortSample) == ncol(selMethyData))
			ddr <- reorder(ddr, sortSample)
		ord <- order.dendrogram(ddr)
		selMethyData <- selMethyData[,ord]
		if (!is.null(expProfile)) {
			if (all(sampleNames(selMethyData) %in% rownames(expProfile))) {
				expProfile <- expProfile[colnames(selMethyData),,drop=FALSE]
			} else {
				expProfile <- expProfile[ord,,drop=FALSE]
			}
		}
	}
	
	## define data track
	dTrack <- DataTrack(range=suppressWarnings(as(locData(selMethyData), 'GRanges')), data=t(assayData(selMethyData)$exprs), 
		chromosome=chromosome, name=dataTrackName, type='heatmap',	gradient=gradient, ncolor=ncolor)		
	allTracks <- c(allTracks, list(dTrack))
	
	## plot Tracks with annotation of DataTrack
	plotInfo <- plotTracksWithDataTrackInfo(allTracks, sampleNames(selMethyData), grange2show=grange2show, dataInfo=expProfile, dataColorMap=expColorMap, labelWidth=0.1, gradient=gradient, ncolor=ncolor, main=main, ...)

	return(invisible(plotInfo))
}


## plot methylation heatmap by gene
##	 selGene: a vector of EntrezIDs or a list of gene2tx
##	 methyGenoSet: a GenoSet object for methylation data
##	 gene2tx: a gene to transcript mapping list, used for retrieving expression.tx data
##	 tx2exon: a transcript to exon mapping list, used for retrieving expression.exon data
##	 expression.tx: an ExpressionSet or data matrix for transcript expression
##	 expression.exon: an ExpressionSet or data matrix for exon expression
##	 phenotype: an ExpressionSet or data.frame for phenotype informaiton
##	 sortBy: sort the samples based on expression, methylation or NA (not sort)
##	 renameExon: whether to rename exons as "transcript_exon1"
##	 showAllTx: whether to show all transcript in gene2tx or just those provided in selGene
##	 includeGeneBody: if FALSE, then only shows the promoter region
##	 CpGInfo: a bed file or GRanges for CpG island information
##	 genomicFeature: used by buildAnnotationTracks function
##	 phenoColor: a vector of colors for pheno types.
plotMethylationHeatmapByGene <- function(selGene, methyGenoSet, gene2tx=NULL, tx2exon=NULL, expression.tx=NULL, expression.exon=NULL, 
	phenotype=NULL, sortBy=c('expression', 'methylation', 'NA'), renameExon=FALSE, showAllTx=TRUE, useBetaValue=TRUE, includeGeneBody=F, CpGInfo=NULL, genomicFeature=NULL, 
	phenoColor=NULL, th=0.99, title.suffix=NULL, addLegend=TRUE, gradient=c("blue", "white", "red"), ncolor=16, main=NULL, newPlot=TRUE, ...) {
	sortBy <- as.character(sortBy)
	sortBy <- match.arg(sortBy)
	
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
	if (!is.null(phenotype)) {
		if (is(phenotype, "ExpressionSet")) {
			phenotype <- exprs(phenotype)
		} 
		if (is.data.frame(phenotype) || is.list(phenotype)) {
			phenotypeLevels <- lapply(phenotype, function(x) {
				x <- levels(as.factor(x))
				return(x)
				})
			phenotype <- as.data.frame(lapply(phenotype, function(x) {
				x <- round(as.numeric(as.factor(x)))
				if (min(x, na.rm=TRUE) <= 0) x <- x - min(x, na.rm=TRUE) + 1
				return(x)
				}))
		}

		otherPhenoName <- colnames(phenotype)
		otherPhenoName <- otherPhenoName[!(otherPhenoName %in% names(phenoColor))]
		if (length(otherPhenoName) > 0) {
			allPhenoName <- names(phenoColor)
			for (phenoName.i in otherPhenoName) {
				maxLevel.i <- max(phenotype[[phenoName.i]], na.rm=TRUE)
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
			expProfile.exon1 <- t(expMatrix.exon[selExon1,,drop=FALSE])
			if (renameExon) {
				colnames(expProfile.exon1) <- paste(selTx, 'exon1', sep='_')
				sigExon1 <- paste(sigTx, 'exon1', sep='_')
			} else {
				sigExon1 <- tx2exon[sigTx]
			}
			expProfile <- cbind(expProfile, expProfile.exon1)
			## sort only based on significant exon1 and transcript
			ord <- order(rowMeans(expProfile[, c(sigExon1, sigTx), drop=FALSE], na.rm=TRUE), decreasing=FALSE)
		} else if (!is.null(expProfile)) {
			## sort only based on significant Tx
			ord <- order(rowMeans(expProfile[,sigTx,drop=FALSE], na.rm=TRUE), decreasing=FALSE)
		}
		
		if (!is.null(phenotype)) {
			## combine the phenotype with the expProfile matrix
			expProfile <- cbind(expProfile, as.matrix(phenotype))
		}

		phenotypeLevel.i <- NULL
		if (!is.null(phenotypeLevels)) {	
			if (selTx[i] %in% names(phenotypeLevels)) {
				phenotypeLevel.i <- phenotypeLevels[[selTx[i]]]
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
			plotWidth <- 1 - legendWidth
			pushViewport(viewport(layout=grid.layout(1, 2, widths=c(plotWidth, legendWidth))))
			pushViewport(viewport(layout.pos.col=1, layout.pos.row=1))
		} 
		
		if (sortBy == 'expression' && !is.null(expProfile)) {
			plotInfo <- heatmapByChromosome(methyGenoSet[, ord,drop=F],	gene.i,	annotationTracks=annotationTracks, expProfile=expProfile[ord,,drop=F], 
																 expColorMap=phenoColor, sortSample=F, dataTrackName='Methylation Profile', main=main, 
																 cex.main=1, ylim=ylim, newPlot=FALSE, gradient=gradient, ncolor=ncolor, ...)
		} else {
			sortSample <- ifelse(sortBy == 'methylation', TRUE, FALSE)
			plotInfo <- heatmapByChromosome(methyGenoSet,	gene.i,	annotationTracks=annotationTracks, expProfile=expProfile, 
																 expColorMap=phenoColor, sortSample=sortSample, dataTrackName='Methylation Profile', main=main, 
																 cex.main=1, ylim=ylim, newPlot=FALSE, gradient=gradient, ncolor=ncolor, ...)
		}
		
		## plot legendInfo
		if (addLegend) {
			popViewport(1)
			## plot legend information
			pushViewport(viewport(layout.pos.col=2, layout.pos.row=1))
			
			## determine the height of legends
			## the height of methylation and expression colorbars are 2*legendHeight
			numOfOtherLegend <- length(which(names(phenoColor) != 'gradient'))
			if (!is.null(expProfile)) {
				legendHeight <- (1 -	(3 + numOfOtherLegend) * 0.05) / (4 + numOfOtherLegend)
			} else {
				legendHeight <- (1 -	(2 + numOfOtherLegend) * 0.05) / (2 + numOfOtherLegend)
			}
			x0 <- 0.15
			colWidth <- 0.15
			
			## plot methylation legend
			stepHeight <- legendHeight * 2 * 0.9 / ncolor
			ystart <- 1 - (0.05 + legendHeight * 2 )
			methyColor <- colorRampPalette(gradient)(ncolor)[1:ncolor]

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
			## add title
			methySubtitle <- ifelse(useBetaValue, '(Beta-value)', '(M-value)')
			grid.text('Methylation', x=0.5, y=max(ytickPos) + 0.025 + as.numeric(convertY(unit(fontsize, 'points'), 'npc')), just=c('center', 'bottom'), 
						default.units = "npc", gp=gpar(fontsize=fontsize, fontface='bold'))
			grid.text(methySubtitle, x=0.5, y=max(ytickPos) + 0.02, just=c('center', 'bottom'), default.units = "npc", gp=gpar(fontsize=fontsize, fontface='bold'))

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
				grid.text('Expression', x=0.5, y=max(ytickPos) + 0.025 + as.numeric(convertY(unit(fontsize, 'points'), 'npc')), just=c('center', 'bottom'), 
							default.units = "npc", gp=gpar(fontsize=fontsize, fontface='bold'))
				grid.text('(log2-RPKM)', x=0.5, y=max(ytickPos) + 0.02, just=c('center', 'bottom'), default.units = "npc", gp=gpar(fontsize=fontsize, fontface='bold'))
			}

			## Add other phenoColor legends
			if (length(phenoColor) > 0) {
				## calculate the stepHeight based on font size,which is defined by the number of color levels
				totalLevel <- length(unlist(phenotypeLevels))
				if (!is.null(expProfile)) {
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
			
			plotInfo <- c(plotInfo, legendWidth=legendWidth)
			popViewport(1)
		}			
		
	}) 
	
	return(invisible(plotResult))		
}

# grid.legend.general <- function (labels, fill, line, pch, frame = TRUE, hgap = unit(0.5, "lines"), 
#		 vgap = unit(0.5, "lines"), default.units = "lines", gp = gpar(), 
#		 draw = TRUE, vp = NULL) 
# {
#		 labels <- as.character(labels)
#		 nkeys <- length(labels)
#		 if (length(pch) != nkeys) 
#				 stop("'pch' and 'labels' not the same length")
#		 if (!is.unit(hgap)) 
#				 hgap <- unit(hgap, default.units)
#		 if (length(hgap) != 1) 
#				 stop("'hgap' must be single unit")
#		 if (!is.unit(vgap)) 
#				 vgap <- unit(vgap, default.units)
#		 if (length(vgap) != 1) 
#				 stop("'vgap' must be single unit")
#		 legend.layout <- grid.layout(nkeys, 3, widths = unit.c(unit(2, 
#				 "lines"), max(unit(rep(1, nkeys), "strwidth", as.list(labels))), 
#				 hgap), heights = unit.pmax(unit(2, "lines"), vgap + unit(rep(1, 
#				 nkeys), "strheight", as.list(labels))))
#		 fg <- frameGrob(layout = legend.layout, vp = vp, gp = gp)
#		 for (i in 1L:nkeys) {
#				 fg <- placeGrob(fg, pointsGrob(0.5, 0.5, pch = pch[i]), 
#						 col = 1, row = i)
#				 fg <- placeGrob(fg, textGrob(labels[i], x = 0, y = 0.5, 
#						 just = c("left", "centre")), col = 2, row = i)
#		 }
#		 if (draw) 
#				 grid.draw(fg)
#		 fg
# }

plotTracksWithDataTrackInfo <- function(trackList, labels, grange2show=NULL, dataTrackName=NULL, dataInfo=NULL, 
			dataColorMap=NULL, dataInfoRange=NULL, labelWidth=0.1, gradient=c("blue", "white", "red"), ncolor=16, main='', newPlot=FALSE, ...) {
	
	if (missing(trackList) || missing(labels)) {
		stop('Please provide "trackList" and "labels"!')
	}
	if (!is(trackList, 'list')) trackList <- list(trackList)
	trackClass <- sapply(trackList, class)
	dataTrackInd <- which(trackClass == 'DataTrack')
	if (length(dataTrackInd) == 0) {
		warning('No DataTrack was found!')
		dataTrackName <- NULL
	} else if (is.null(dataTrackName)){
		dataTrackName <- names(trackList[[dataTrackInd[1]]])
	}
	
	## start plotting 
	if (newPlot)	grid.newpage()
	
	## determine the layout based on expression profile
	# labelWidth <- 0.1
	if (!is.null(dataInfo)) {
		## 
		if (!is.matrix(dataInfo)) {
			dataInfo <- matrix(dataInfo, ncol=1)
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
	plotInfo <- plotTracks(trackList, from=start(grange2show)[1], 
			to=end(grange2show)[1], chromosome=chromosome,	add=TRUE, main=main, ...)
	## retrieve the plot coordinates		
	plotLoc <- coords(plotInfo$title)
	popViewport(1)

	## plot annotation or phenotype information
	pushViewport(viewport(layout.pos.col=2, layout.pos.row=1))

	defaultPar <- Gviz:::.parMappings$GdObject		
	## calculate the label positions 
	hh <- (plotLoc[dataTrackName, 'y2'] - plotLoc[dataTrackName, 'y1'])/(plotLoc[nrow(plotLoc), 'y2'] + plotLoc[1, 'x1'])
	y0 <- (plotLoc[dataTrackName, 'y1'])/(plotLoc[nrow(plotLoc), 'y2'] + plotLoc[1, 'x1'])
	realHeight <- hh / 1.1
	# grid.rect(0, y=1 - y0 - hh/2, width=1, height=realHeight, just=c('left', 'center'), gp=gpar(lty='dashed', col=3))
	ystart <- 1 - y0 - hh/2 - realHeight / 2
	yend <- 1 - y0 - hh/2 + realHeight / 2
	stepHeight <- realHeight / length(labels)
	x0 <- 0.05
	colWidth <- 0.1
	if (!is.null(dataInfo)) {
		# plot the colnames first
		colWidth <- 1 / (5 + ncol(dataInfo))
		if (!is.null(colnames(dataInfo))) {
			fontsize <- floor(as.numeric(convertX(unit(colWidth * 0.99, 'npc'), 'points')))
			grid.text(colnames(dataInfo), x=x0 + colWidth * (1:ncol(dataInfo)), y = yend + 0.01, rot=90, just=c('left', 'bottom'), 
				gp=gpar(fontfamily=defaultPar$fontfamily, fontsize=fontsize, fontface=defaultPar$fontface, col=1))
		}
		
		dataColor <- vector(mode='list', length=ncol(dataInfo))
		names(dataColor) <- colnames(dataInfo)
		gradient.data <- gradient
		if (!is.null(dataColorMap)) {
			if (!is.list(dataColorMap)) stop('dataColorMap should be a named list!')
			
			if ('gradient' %in% names(dataColorMap)) {
				gradient.data <- dataColorMap$gradient
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
				data.i <- dataInfo[,dataCol.i]
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
				valsScaled <- Gviz:::.z2icol(dataInfo[,otherCol, drop=FALSE], ncolor, xrange=range(dataInfo[,otherCol]))
				for (i in 1:length(otherCol)) {
					dataColor[[otherCol[i]]] <- colorRampPalette(gradient.data)(ncolor)[valsScaled[,i]]
				}
			}
		} else {
			## dataInfoRange is to control the plot color range
			if (!is.null(dataInfoRange)) {
				valsScaled <- Gviz:::.z2icol(dataInfo, ncolor, xrange=dataInfoRange)
			} else {
				valsScaled <- Gviz:::.z2icol(dataInfo, ncolor, xrange=range(dataInfo))	
			}
			dataColorMap <- colorRampPalette(gradient)(ncolor)
			for (i in 1:ncol(dataInfo)) {
				dataColor[[i]] <- colorRampPalette(gradient)(ncolor)[valsScaled[,i]]
			}
		}

		## plot data matrix as a regular heatmap
		for (i in 1:ncol(dataInfo)) {
			dataCol.i <- colnames(dataInfo)[i]
			grid.rect(x0, stepHeight * (1:length(labels)) + ystart, width=colWidth, height=stepHeight, 
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
		realLabelWidth <- 0.9
	}
	fontsizeW <- round(as.numeric(convertX(unit(realLabelWidth / 5, 'npc'), 'points')))
	fontsize <- min(fontsize, fontsizeW)
	grid.text(labels, x=x0, y = stepHeight * (1:length(labels)) + ystart - stepHeight/2, just=c('left', 'center'), gp=gpar(fontfamily=defaultPar$fontfamily, fontsize=fontsize, fontface=defaultPar$fontface, col=3))
	popViewport(2)
	
	plotInfo <- c(plotInfo, labelWidth=labelWidth)

	return(invisible(plotInfo))
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
		chrName[ind] <- paste('chr', chrName[ind], sep='')
	}
	
	if (is(grange, 'GRanges')) {
		seqlevels(grange) <- chrName
	} else if (is(grange, 'character')) {
		grange <- chrName
	} else if (is(grange, 'GenoSet')) {
		names(locData(grange)) <- chrName
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

