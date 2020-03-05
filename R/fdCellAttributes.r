

chunkFinder=function(someRaster,cell.ind,chunkFiles){
	notNAs=which(!is.na(values(someRaster)))
	keep=which(cell.ind$cellID %in% notNAs )
	chunks=unique(cell.ind[keep,'chunkID'])
	keep=mapply(function(x){grep(paste0('chunk_',x), basename(chunkFiles))},chunks)
	chunkFiles[keep]
}

subsetRangeByCellProportion=function(someRaster,cbs.f,cell.ind,mc.cores=1){
	rangeSize.tmp=mclapply(seq_along(cbs.f),function(x){
		cbs=readRDS(cbs.f[x])
		keep=	which(cell.ind$cellID %in% cbs@Dimnames[[1]] )
		keep1=cell.ind$cellID[keep]
		someRaster.mat=matrix(values(someRaster),ncol=nlayers(someRaster))[keep1,] #matrix(values(someRaster),ncol=1)[keep1,] # for 1 layer
		someRaster.mat[is.na(someRaster.mat)]=0
		out=t(cbs) %*% someRaster.mat
		message(x)
		tmp=data.frame(as.matrix(out))
		names(tmp)=names(someRaster) 
		tmp
	},mc.cores=mc.cores)
			# assumes the columns line up perfectly
	#aa=data.frame(rangeSize=apply(do.call('cbind', rangeSize.tmp), 1,sum))  # for 1 layer
	all.rs.burned.by.year=data.frame(do.call('cbind', lapply(1:ncol(rangeSize.tmp[[1]]),function(x){
			thisScen=do.call('cbind',lapply(rangeSize.tmp,function(y) y[[x]]))
			apply(thisScen,1,sum)
		})))
	names(all.rs.burned.by.year)=names(someRaster)#basename(ba.f)
	all.rs.burned.by.year
}