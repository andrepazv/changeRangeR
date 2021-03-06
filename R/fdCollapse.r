#===================================================================
#===================================================================
#===================================================================
#' @title Collapse columns of a sparse matrix
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

collapseCols=function(inDir,
                      outDir,
                      scenario,
                      spAttrTable,
                      attrName,
                      type='binary',
                      keepChunks=F,
                      mc.cores=mc.cores,
                      overwrite=F,
                      verbose=T){
  ## for testing
  ## attrTable=genusSpAttr; attrName='genus'; scenario='present'; inDir=sumDirs$cbsDir; outDir=sumDirs$cbgDir; keepChunks=F

  t1=proc.time()
  message(paste0('starting ',scenario))
  cbs.f=changeRangeR:::.getCBS(inDir,scenario)
  if(Sys.info()["sysname"]== "Windows") mclapply=parallelsugar::mclapply
  if(Sys.info()["sysname"]!= "Windows") mclapply=parallel::mclapply

  anames=unique(spAttrTable[,attrName]) %>% na.omit
  message(paste(length(anames), 'groups to aggregate over'))
  out=mclapply(seq_along(cbs.f), function(x){
    if(verbose) message(x)
    out.f=paste0(outDir,'/chunk_',x,'.rds')
    if(file.exists(out.f) & !overwrite) {
    	message(paste('chunk',x,'done, skipping'))
    	return()
    }
    cbs.tmp=readRDS(cbs.f[x])
    colByName=lapply(seq_along(anames),function(y){
      keep=which(spAttrTable[,attrName]==anames[y])
      if(length(keep)==1) return(cbs.tmp[,keep])
      o=cbs.tmp[,keep] %>% textTinyR::sparse_Sums(rowSums = T)
      if(type=='binary') return(o %>% pmin(1))
      if(type=='count') return(o)
    }) # returns the cols of the cbs matrix where each is collapsed to a unique attribute name
    sp.mat=do.call(cbind,colByName) %>% Matrix::Matrix(sparse=T)
    rownames(sp.mat)=rownames(cbs.tmp)
    colnames(sp.mat)=anames
    if(keepChunks) saveRDS(sp.mat,file=paste0(outDir,'/chunk_',x,'.rds'))
    sp.mat
  },mc.cores=mc.cores)
  if(!keepChunks){
    out1=do.call(rbind,out)
    saveRDS(out1,file=paste0(outDir,'/chunk_all.rds'))
  }

  t2=proc.time()-t1
  message( paste0(round(t2[3],2),' s') )
  return()
}

#===================================================================
#===================================================================
#===================================================================
#' @title Collapse rows of a sparse matrix
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

collapseRows=function(inDir,
                      outDir,
                      scenario,
                      cellAttrTable,
                      attrName,
                      type='binary',
                      keepChunks=F,
                      mc.cores=mc.cores,
                      verbose=F){
  ## for testing
  ## cellAttrTable=cellAttr; attrName='ecoregion'; scenario='present'; inDir=sumDirs$cbsDir; outDir=sumDirs$cbgDir; keepChunks=F; type='binary'

  t1=proc.time()
  message(paste0('starting ',scenario))
  cbs.f=changeRangeR:::.getCBS(inDir,scenario)
  if(Sys.info()["sysname"]== "Windows") mclapply=parallelsugar::mclapply
  if(Sys.info()["sysname"]!= "Windows") mclapply=parallel::mclapply

  anames=unique(cellAttrTable[,attrName]) %>% na.omit
  out=mclapply(seq_along(cbs.f), function(x){
    if(verbose) message(x)
    cbs.tmp=readRDS(cbs.f[x])
    rowByName=lapply(seq_along(anames),function(y){
      # make sure to get the cellID (not just the row number from the attribute table)
      keep=cellAttrTable[which(cellAttrTable[,attrName]==anames[y]),'cellID']
      keep.in.this.chunk=which(keep %in% as.numeric(rownames(cbs.tmp)))
      # when you collapse spatially reference them  based on the cell attribute table, where there's a field linking the (e.g.) ecoregion ID to the row name of the ecoregion by species matrix
      if(length(keep.in.this.chunk)==0) { return(rep(0,ncol(cbs.tmp)))}
      o=cbs.tmp[keep.in.this.chunk,] %>% textTinyR::sparse_Sums(rowSums = F)
      if(type=='binary') return(o %>% pmin(1))
      if(type=='count') return(o)
    }) # returns the cols of the cbs matrix where each is collapsed to a unique attribute name
    sp.mat=do.call(rbind,rowByName) %>% Matrix::Matrix(sparse=T)
    not0=which(apply(sp.mat,1,sum)>0)
    #if(all(null.attrs)) return(NULL)
    rownames(sp.mat)=anames
    colnames(sp.mat)=colnames(cbs.tmp)
    # probably a bad option because the same spatial unit can be distributed across a few chunks
    if(keepChunks) saveRDS(sp.mat,file=paste0(outDir,'/chunk_',x,'.rds'))
    sp.mat
  },mc.cores=mc.cores)
  # need to combine spatial units which are split across chunks. just use reduce, because all the output matrices should be the same size
  out1=Reduce('+',out)

  if(!keepChunks){
    saveRDS(out1,file=paste0(outDir,'/chunk_all.rds'))
  }

  t2=proc.time()-t1
  message( paste0(round(t2[3],2),' s') )
  return()
}


#unfinished first try to be deleted when the real thing works
# cbsCollapseCells=function(cbsDir,scenario,cellAttr,attrName,verbose=F){
#
#   t1=proc.time()
#   message(paste0('starting ',scenario))
#   cbs.f=changeRangeR:::.getCBS(cbsDir,scenario)
#   if(Sys.info()["sysname"]== "Windows") mclapply=parallelsugar::mclapply
#   if(Sys.info()["sysname"]!= "Windows") mclapply=parallel::mclapply
#   attrs=unique(cellAttr[,attrName])
#
#   richByFact=mclapply(seq_along(cbs.f), function(fu){
#     if(verbose) message(fu)
#     cbs.tmp=readRDS(cbs.f[fu])
#     # 1. split CBS by factor
#     a=cellAttr %>% filter(chunkID==fu)
#     ua=unique(a[,attrName]) %>% na.omit
#     if(length(ua)==0) return(NULL)
#     # 2. colSums within factor
#     cbs %>% mutate(cellID=row.names(cbs)) %>% left_join(cellAttr,by='cellID') %>%
#       group_by(ecoregion) %>% summarize(somethingAboutColSums) %>% somethingToMakeSparseAgain
#     # make a sparse matrix where the rows are the factor values and cols are species
#     out=lapply(ua,function(xx){
#       cid=cellAttr %>% filter(.data[[attrName]]==xx) %>% select(cellID)
#       # match cellIDs to the rowNames of cbs
#       keep=which(as.integer(row.names(cbs.tmp)) %in% cid[,1])
#       out=textTinyR::sparse_Sums(cbs.tmp[keep,], rowSums = F)
#       # 3. convert to binary with vals > 0
#       out[out>0]=1 #there's probably a faster way
#       out
#     }) %>% data.frame
#     names(out)=ua
#     out
#   },mc.cores=mc.cores)
#   # 4. combine across cbs
#   str(richByFact)
#   richByFact %>% reduce(bind_rows)
#   richByFact %>% bind_rows
#
#   out=Reduce('+',richByFact) %>% as.matrix %>% as.data.frame
#   names(out)=names(someStack)
#   out=data.frame(sp.ind,out)
#   # 5. row sums for richness (or any other operation on these reduced factor by species (FBS) matrices)
#   data.frame(cellID=as.numeric(rownames(cbs)),rich=textTinyR::sparse_Sums(cbs, rowSums = T))
#
# }
