#===================================================================
#===================================================================
#===================================================================
#'  Create a table linking species names with their index in CBS matrices; can optionally be used as a template for a	`speciesAttributeTable`
#' @notes Columns: species name, integer index
#' @export

speciesIndexTable=function(allSpeciesMaps,sumDirs){
	# old way with lists
	#sp.ind=data.frame(species=unique(unlist(inputsFromSDMWorkflow$sp.names, recursive=T)))
	sp.ind=data.frame(species=allSpeciesMaps$sp.names)
	sp.ind$index=1:nrow(sp.ind)
	saveRDS(sp.ind,file=paste0(sumDirs$sumBaseDir,'/speciesIndexTable.rds'))
	sp.ind
}

#===================================================================
#===================================================================
#===================================================================
#' @title Calculate taxon richness
#' @description
#' @param
#' @param
#' @param
#' @param
#' @param
#' @param
#' @param
#' @param
# @examples
#
#' @return
#' @author Cory Merow <cory.merow@@gmail.com>
#' @note optionally uses an attribute table to subset by each column not called `species` or `index` and creates a richness file based on the name of the attribute table column. the attribute table must include the species and index columns, as generated by `speciesIndexTable`
# @seealso
# @references
# @aliases - a list of additional topic names that will be mapped to
# this documentation when the user looks them up from the command
# line.
# @family - a family name. All functions that have the same family tag will be linked in the docum
#' @export
richnessFromCBS=function(cbsDir,
												 scenario,
												 envGrid,
												 mc.cores=1,
												 attrTable=NULL,
												 attrName=NULL,
											   outDir=NULL,
											   verbose=F){

	t1=proc.time()
	message(paste0('starting ',scenario))
	cbs.f=changeRangeR:::.getCBS(cbsDir,scenario)
	if(Sys.info()["sysname"]== "Windows") mclapply=parallelsugar::mclapply
	if(Sys.info()["sysname"]!= "Windows") mclapply=parallel::mclapply

	if(is.null(attrTable)){
		richByCell=mclapply(seq_along(cbs.f), function(x){
			if(verbose) message(x)
			cbs=readRDS(cbs.f[x])
			a=data.frame(cellID=as.numeric(rownames(cbs)),
								 	 rich=textTinyR::sparse_Sums(cbs, rowSums = T))
		},mc.cores=mc.cores)
		rich.vec=do.call('rbind',richByCell)
		rich.r=raster(envGrid[[1]])
		values(rich.r)[rich.vec$cellID]= rich.vec$rich
    names(rich.r)=scenario
		if(!is.null(outDir))
		  writeRaster(rich.r,file=paste0(outDir,'/richness_', scenario,'.tif'), overwrite=T)
		t2=proc.time()-t1
		message( paste0(round(t2[3],2),' s') )
		return(rich.r)
	} else {
		# CM: not needed now, since manually specified
	  # attrNames=names(attrTable)
		# # if species and index were included, toss them
		# toss=unlist(mapply(function(x){grep(x,attrNames)}, c('species','index')))
		# if(length(toss) > 0 ) attrNames=attrNames[-toss]

	  # make dummy cols for each factor
	  dum=attrTable %>% select(!!attrName) %>% fastDummies::dummy_cols() %>% select(-!!attrName)

		out=lapply(seq_along(dum),function(y){
			if(verbose) message(names(dum)[y])
			keep=attrTable$index[dum[y]==1]
			richByCell=mclapply(seq_along(cbs.f), function(x){
				if(verbose) message(x)
				cbs.tmp=readRDS(cbs.f[x])
				cbs=cbs.tmp[,keep]
				data.frame(cellID=as.numeric(rownames(cbs)),rich=textTinyR::sparse_Sums(cbs, rowSums = T))
			},mc.cores=mc.cores)
			rich.vec=do.call('rbind',richByCell)
			rich.r=raster(envGrid[[1]])
			values(rich.r)[rich.vec$cellID]= rich.vec$rich

			if(!is.null(outDir))
			  writeRaster(rich.r,file=paste0(outDir,'/', names(dum[y]),'.tif'), overwrite=T)
			rich.r
		})
		out1=stack(out)
		names(out1)=names(dum)
		t2=proc.time()-t1
		message( paste0(round(t2[3],2),' s') )
		return(out1)
	}
}

#===================================================================
#===================================================================
#===================================================================
# #' @notes optionally uses an attribute table to subset by each column not called `species` or `index` and creates a richness file based on the name of the attribute table column. the attribute table must include the species and index columns, as generated by `speciesIndexTable`
# #' @export
# richnessFromCBSAttr=function(cbsDir,
# 														 scenario,
# 														 env,
# 														 mc.cores=1,
# 														 attrTable=NULL,
# 														 outDir=NULL,
# 														 verbose=F){
# # 	cbs.f=list.files(paste0(cbsDir,'/',scenario),full.names=T)
# # 	toss=grep('temp_long',cbs.f)
# # 	cbs.f=cbs.f[-toss]
# 	cbs.f=list.files(paste0(cbsDir,'/',scenario),full.names=T,pattern='chunk')
# 	attrNames=names(attrTable)
# 	# if species and index were included, toss them
# 	toss=unlist(mapply(function(x){grep(x,attrNames)}, c('species','index')))
# 	if(length(toss) > 0 ) attrNames=attrNames[-toss]
#
# 	if(Sys.info()["sysname"]== "Windows") mclapply=parallelsugar::mclapply
# 	if(Sys.info()["sysname"]!= "Windows") mclapply=parallel::mclapply
#
# 	out=lapply(seq_along(attrNames),function(y){
# 		if(verbose) message(attrNames[y])
# 		keep=attrTable$index[attrTable[attrNames[y]]==1]
# 		richByCell=mclapply(seq_along(cbs.f), function(x){
# 		  if(verbose) message(x)
# 			cbs.tmp=readRDS(cbs.f[x])
# 			cbs=cbs.tmp[,keep]
# 			#fuck=data.frame(spID=as.numeric(colnames(cbs)),rich=textTinyR:: sparse_Sums(cbs, rowSums = F))
# 			#print(fuck[fuck[,2]>0,])
# 			data.frame(cellID=as.numeric(rownames(cbs)),rich=textTinyR:: sparse_Sums(cbs, rowSums = T))
# 		},mc.cores=mc.cores)
# 		rich.vec=do.call('rbind',richByCell)
# 		rich.r=env[[1]]
# 		values(rich.r)=NA
# 		values(rich.r)[rich.vec$cellID]= rich.vec$rich
# 		if(!is.null(outDir))		writeRaster(rich.r,file=paste0(outDir,'/', attrNames[y],'.tif'), overwrite=T)
# 		rich.r
# 	})
# 	out1=stack(out)
# 	names(out1)=attrNames
#   out1
# }

#===================================================================
#===================================================================
#===================================================================
#' @title Calculate taxon range area
#'
#' @description
#' @param
#' @param
#' @param
#' @param
#' @param
#' @param
#' @param
#' @param
# @examples
#
#' @return
#' @author Cory Merow <cory.merow@@gmail.com>
#' @note Units of range area are number of pixels. So conversion to area requires an equal area projection. If you're not using an equal area projection, consider using the area of the cell as a cell attribute with `raster::area`
# @seealso
# @references
# @aliases - a list of additional topic names that will be mapped to
# this documentation when the user looks them up from the command
# line.
# @family - a family name. All functions that have the same family tag will be linked in the docum
#' @export

rangeArea=function(cbsDir,outDir=NULL,scenario,sp.ind,mc.cores=1,verbose=F){
	t1=proc.time()
	message(paste0('starting ',scenario))
	cbs.f=changeRangeR:::.getCBS(cbsDir,scenario)
	if(Sys.info()["sysname"]== "Windows") mclapply=parallelsugar::mclapply
	if(Sys.info()["sysname"]!= "Windows") mclapply=parallel::mclapply

	rangeSize.tmp=mclapply(seq_along(cbs.f),function(x){
	  if(verbose) message(x)
		cbs=readRDS(cbs.f[x])
		data.frame(rich=textTinyR::sparse_Sums(cbs, rowSums = F))
	},mc.cores=mc.cores)
		# assumes the columns line up perfectly
	rangeSize=data.frame(sp.ind,rangeArea=apply(do.call('cbind', rangeSize.tmp), 1,sum))
	if(!is.null(outDir)) saveRDS(rangeSize,file=paste0(outDir,'/RangeSize_',scenario, '.rds'))
	# check against base range maps (leaving this here in case we find a need)
	# checkVsRangeMaps=F
	# 	if(scenario=='Present' & checkVsRangeMaps){
	# 		range.f=list.files( '/Users/ctg/Documents/SDMs/BIEN41/NWPlants_BinaryOnly/BIEN41_outputs/PPM/BinaryMaps',full.names=T,recursive=T,pattern='TP05')
	# 		keep=sample(seq_along(range.f),50)
	# 		out=mclapply(seq_along(keep),function(x){
	# 			r=raster(range.f[keep[x]])
	# 			true=cellStats(r,sum,na.rm=T)
	# 			sp=strsplit(tools::file_path_sans_ext(basename(range.f[keep[x]])), '__')[[1]][2]
	# 			keep1=grep(sp, rangeSize$species)
	# 			print(x)
	# 			data.frame(sp=sp,true=true,cbsRangesize=rangeSize[keep1,3])
	# 		},mc.cores=5)
	# 		(out1=do.call(rbind,out))
	# 	}
	t2=proc.time()-t1
	message(paste0(round(t2[3],2),' s'))
	rangeSize
}


#===================================================================
#===================================================================
#===================================================================
#' @export
rangeAreaInBinaryStack=function(someStack,cbsDir,scenario,sp.ind,cell.ind,mc.cores=1){
    #  for testing
  	#  someStack=someRaster2

	t1=proc.time()
	message(paste0('starting ',scenario))

	cbs.f=changeRangeR:::.getCBS(cbsDir,scenario)
	bigMat=values(someStack)[cell.ind$cellID,]
	#bigMat=matrix(values(someStack)[cell.ind$cellID,],ncol=nlayers(someStack))
	# bigMat= someStack %>% values %>% as.matrix # seems as slow
	rangeSize.tmp=lapply(seq_along(cbs.f),function(x){
		message(x)
		cbs=readRDS(cbs.f[x])
		# row numbers of cell.ind for cells in this cbs matrix
		keep=	which(cell.ind$cellID %in% as.numeric(dimnames(cbs)[[1]]) )
		# cellIDs we need from the stack
		keep1=cell.ind$cellID[keep]
		someStack.mat=bigMat[keep1,]
		#matrix(values(someStack),ncol=1)[keep1,] # for 1 layer
		someStack.mat[is.na(someStack.mat)]=0
		out=t(cbs) %*% someStack.mat
    out
	})#,mc.cores=mc.cores)
	out=Reduce('+',rangeSize.tmp) %>% as.matrix %>% as.data.frame
	names(out)=names(someStack)
	out=data.frame(sp.ind,out)

			# assumes the columns line up perfectly
	#aa=data.frame(rangeSize=apply(do.call('cbind', rangeSize.tmp), 1,sum))  # for 1 layer
# 	all.rs.burned.by.year=data.frame(do.call('cbind', lapply(1:ncol(rangeSize.tmp[[1]]),function(x){
# 			thisScen=do.call('cbind',lapply(rangeSize.tmp,function(y) y[[x]]))
# 			apply(thisScen,1,sum)
# 		})))
# 	names(all.rs.burned.by.year)=names(someStack)#basename(ba.f)
	t2=proc.time()-t1
	message(paste0(round(t2[3],2),' s'))
	out
}

# different version i can probably toss
# # #for(i in 1:length(inputsFromSDMWorkflow$scenarios)){
# # #	scenario=inputsFromSDMWorkflow$scenarios[i]
# # #	print(scenario)
# # 	scenario="Present"
# # 	cbs.f=list.files(paste0(sumDirs$cbsDir,'/',scenario),full.names=T)
# # 	toss=grep('temp_',cbs.f)
# # 	cbs.f=cbs.f[-toss]
# # 	rangeSizeBurned.tmp=mclapply(seq_along(cbs.f),function(x){
# # 		message(x)
# # 		cbs=readRDS(cbs.f[x])
# # 		# subset to just the cells in this chunk
# # 		# this should work, but sometimes a cell is missing
# # 		# burn.sub=burn.val[cell.ind$cellID[cell.ind$chunkID==x],]
# # 		burn.sub=burn.val[as.numeric(dimnames(cbs)[[1]]),]
# # 		t(cbs) %*% burn.sub
# # 	},mc.cores=mc.cores)
# # 	out=Reduce('+',rangeSizeBurned.tmp) %>% as.matrix %>% as.data.frame
# # 		# assumes the columns line up perfectly
# # 	rangeSizeBurned=data.frame(sp.ind,out)
# #
# # 	rangeSizeOutDir=paste0(summaryBaseDir,'/RangeSize')
# # 	if(!file.exists(rangeSizeOutDir)) dir.create(rangeSizeOutDir)
# # 	saveRDS(rangeSizeBurned, file=paste0(rangeSizeOutDir,'/RangeSizeBurnt_',scenario,version,'.rds'))
# # #}

#===================================================================
#===================================================================
#===================================================================
#' @title Calculate range area intersecting a raster
#'
#' @description For use with a single categorical raster
#' @param
#' @param
#' @param
#' @param
#' @param
#' @param
#' @param
#' @param
# @examples
#
#' @return
#' @author Cory Merow <cory.merow@@gmail.com>
#' @note
# @seealso
# @references
# @aliases - a list of additional topic names that will be mapped to
# this documentation when the user looks them up from the command
# line.
# @family - a family name. All functions that have the same family tag will be linked in the docum
#' @export

rangeAreaCategoricalRaster=function(someRaster,catNames=NULL,cbsDir,
                                    scenario,sp.ind,cell.ind, mc.cores=1){
	#  for testing
  #  someRaster=ecoAus2; cbsDir=sumDirs$cbsDir; catNames=NULL

  t1=proc.time()
	cats=sort(na.omit(unique(values(someRaster))))
  n.layers=length(unique(values(ecoAus2)))
  if(n.layers > 50) message(paste0("Warning - you're asking to layerize ",
                                   n.layers," which can be pretty slow."))
	someRaster2=layerize(someRaster)
	if(!is.null(catNames)) names(someRaster2)=catNames
	out=rangeAreaInBinaryStack(someRaster2,cbsDir=cbsDir,scenario=scenario,
	                           sp.ind=sp.ind, cell.ind=cell.ind,mc.cores=mc.cores)

	t2=proc.time()-t1
	message(paste0(round(t2[3],2),' s'))
	out
}

#===================================================================
#===================================================================
#===================================================================
# CM: an initial try to make cellAttributes work with dataframes, but probably not needed because rangeAreaCategoricalRaster already does this
#' @param	cellAttr a cell attribute table (`data.frame`)
#' @param	toToss a vector of column names to omit from the cellAttribute table. The default value omits the automatically generated indexing columns used internally by `changeRanger`
# @example
# rbc=rangeByCellAttr(cellAttr=cellEcoAttr, toToss=c('x','y','cellID','chunkID','cellChunkID'),sp.ind)
#
# rangeByCellAttr=function(cellAttr,
# 												 toToss=c('x','y','cellID','chunkID','cellChunkID'),
# 												 sp.ind){
# 	#  for testing
# 	#  cellAttr=cellEcoAttr; toToss=c('x','y','cellID','chunkID','cellChunkID')
#
# 	t1=proc.time()
# 	cbs.f=.getCBS(cbsDir,scenario)
#
# 	if(Sys.info()["sysname"]== "Windows") mclapply=parallelsugar::mclapply
# 	if(Sys.info()["sysname"]!= "Windows") mclapply=parallel::mclapply
#
# 	tmp=mclapply(seq_along(cbs.f),function(x){
# 		message(paste0('chunk ',x))
# 		cbs=readRDS(cbs.f[x])
# 		keep=which(cellAttr$cellID %in% as.numeric(dimnames(cbs)[[1]]) )
# 		# not ideal that i manually remove columns by name.
# 		keepCol=which(!names(cellAttr) %in% toToss)
# 		some.mat=as.matrix(cellAttr[keep,keepCol] )
# 		some.mat[is.na(some.mat)]=0
# 		out=t(cbs) %*% some.mat # species X #cells in ecoregion
# 		out
# 	},mc.cores=mc.cores)
# 	out=Reduce('+',tmp) %>% as.matrix %>% as.data.frame
#
# 	names(out)=names(cellAttr)[which(!names(cellAttr) %in% toToss)]
# 	out=data.frame(sp.ind,out)
# 	t2=proc.time()-t1
# 	message(paste0(round(t2[3],2),' s'))
# 	out
# }

#===================================================================
#===================================================================
#===================================================================

# DONT USE THIS; IT IS REPLACED BY rangeAreaInBinaryStack
# @export
# subsetRangeByCellProportion=function(someRaster,cbs.f,cell.ind,mc.cores=1){
# 	rangeSize.tmp=mclapply(seq_along(cbs.f),function(x){
# 		cbs=readRDS(cbs.f[x])
# 		keep=	which(cell.ind$cellID %in% cbs@Dimnames[[1]] )
# 		keep1=cell.ind$cellID[keep]
# 		someRaster.mat=matrix(values(someRaster),ncol=nlayers(someRaster))[keep1,] #matrix(values(someRaster),ncol=1)[keep1,] # for 1 layer
# 		someRaster.mat[is.na(someRaster.mat)]=0
# 		out=t(cbs) %*% someRaster.mat
# 		message(x)
# 		tmp=data.frame(as.matrix(out))
# 		names(tmp)=names(someRaster)
# 		tmp
# 	},mc.cores=mc.cores)
# 			# assumes the columns line up perfectly
# 	#aa=data.frame(rangeSize=apply(do.call('cbind', rangeSize.tmp), 1,sum))  # for 1 layer
# 	all.rs.burned.by.year=data.frame(do.call('cbind', lapply(1:ncol(rangeSize.tmp[[1]]),function(x){
# 			thisScen=do.call('cbind',lapply(rangeSize.tmp,function(y) y[[x]]))
# 			apply(thisScen,1,sum)
# 		})))
# 	names(all.rs.burned.by.year)=names(someRaster)#basename(ba.f)
# 	all.rs.burned.by.year
# }


