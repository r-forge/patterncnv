#===  Print session/data information

.patCNV.ShowSession <- function(obj)
{
  cat('An object of "PatCNVSession" class \n')
  cat('\t Total number of exons:', obj@N.exon,'\n')
  cat('\t Total number of samples:', obj@N.sample,'\n')
  cat('\t Sample type distribution:\n')
  .patCNV.print.table(obj@sample.info$sample.type)
  
  if(file.exists(obj@pattern.file$avg_wig) & 
       file.exists(obj@pattern.file$avg_wig) & 
        file.exists(obj@Misc$median_RPKM_file) &
          file.exists(obj@Misc$mean_RPKM_file)) {
    cat('Pattern files have been generated \n')
  } else {
    cat('Pattern files have not been generated yet \n')
  }
    
  cat('\n')
}


.patCNV.ShowData <- function(obj)
{
  cat('An object of "PatCNVData" class \n')
  cat('\t Total number of exons:', obj@N.exon,'\n')
  cat('\t Total number of samples:', obj@N.sample,'\n')
  cat('\t Sample type distribution:\n')
  #.patCNV.print.table(obj@sample.info$sample.type)
  .patCNV.print.table(colData(obj@eSet)$sample.type)
  cat('\n')
  if(!is.element("covg",obj@datatype.vec))
  {cat('Coverage information is empty \n')} else{
    cat('Coverage information has been generated \n')  
  }
  if(!is.element("CNV",obj@datatype.vec))
  {cat('CNV information is empty \n')} else{
    cat('CNV information has been generated \n')  
  }
}

setMethod('summary','PatCNVSession', function(object)     .patCNV.ShowSession(object))
setMethod('show', 'PatCNVSession', function(object)     .patCNV.ShowSession(object))

setMethod('summary','PatCNVData', function(object)     .patCNV.ShowData(object))
setMethod('show', 'PatCNVData', function(object)     .patCNV.ShowData(object))


#===  exonNum 

setMethod("exonNum", c("obj" = "PatCNVData"),
          function(obj)
          {           return(obj@N.exon)       } ) 
setMethod("exonNum", c("obj" = "PatCNVSession"),
          function(obj)
          {           return(obj@N.exon)       } ) 


#===  sampleNum

setMethod("sampleNum", c("obj" = "PatCNVData"),
          function(obj)
          {           return(obj@N.sample)       } ) 
setMethod("sampleNum", c("obj" = "PatCNVSession"),
          function(obj)
          {           return(obj@N.sample)       } )

#===  exonInfo

.getExonInfo <- function(obj,gene.name,capture.only,attri)
{
  if(is.null(gene.name))
  { sel.exon.idx <- seq(1,exonNum(obj)) } else {
    sel.exon.idx <- which(is.element(obj@exon.info$gene.name,gene.name))
  }
  tmp.list <- list()
  tmp.list$gene.name <- obj@exon.info$gene.name[sel.exon.idx]
  tmp.list$chr <- as.character(obj@exon.info@seqnames[sel.exon.idx])
  tmp.list$start <- obj@exon.info@ranges@start[sel.exon.idx]
  tmp.list$end <- obj@exon.info@ranges@start[sel.exon.idx] + 
                      obj@exon.info@ranges@width[sel.exon.idx]
  tmp.list$is.captured <- as.numeric(obj@exon.info$is.captured[sel.exon.idx])
  tmp.list$N.exon.bin <- as.numeric(obj@exon.info$N.exon.bin[sel.exon.idx])
  #sim.session@exon.info$N.exon.bin
  
  if(is.null(attri))
  {
    # default attributes to return 
    attri <- c("gene.name","chr","start","end","is.captured") 
    tmp.list <- data.frame(tmp.list)[,attri]  
  } else {
    tmp.list <- data.frame(tmp.list)[,attri]
  }
  
  if(capture.only==TRUE)  
    { captured.idx <- which(tmp.list$is.captured==1)
      return(tmp.list[captured.idx,])
    } else {
      return(tmp.list)
    }
}


.getExonInfoData <- function(obj,gene.name,capture.only,attri)
{
  if(is.null(gene.name))
  { sel.exon.idx <- seq(1,exonNum(obj)) } else {
    sel.exon.idx <- which(is.element(rowData(obj@eSet)$gene.name,gene.name))
  }
  tmp.list <- list()
  tmp.list$gene.name <- rowData(obj@eSet)$gene.name[sel.exon.idx]
  tmp.list$chr <- as.character(rowData(obj@eSet)@seqnames[sel.exon.idx])
  tmp.list$start <- rowData(obj@eSet)@ranges@start[sel.exon.idx]
  tmp.list$end <- rowData(obj@eSet)@ranges@start[sel.exon.idx] + 
    rowData(obj@eSet)@ranges@width[sel.exon.idx]
  tmp.list$is.captured <- as.numeric(rowData(obj@eSet)$is.captured[sel.exon.idx])
  tmp.list$N.exon.bin <- as.numeric(rowData(obj@eSet)$N.exon.bin[sel.exon.idx])
    
  if(is.null(attri))
  {
    # default attributes to return 
    attri <- c("gene.name","chr","start","end","is.captured") 
    tmp.list <- data.frame(tmp.list)[,attri]  
  } else {
    tmp.list <- data.frame(tmp.list)[,attri]
  }
  
  tmp.list <- data.frame(tmp.list)
  if(capture.only==TRUE)  
  { captured.idx <- which(tmp.list$is.captured==1)
    return(tmp.list[captured.idx,])
  } else {
    return(tmp.list)
  }
}

setMethod("exonInfo", c("obj" = "PatCNVSession"),
          function(obj, gene.name=NULL, capture.only=FALSE, attri=NULL)
          {    return(.getExonInfo(obj,gene.name,capture.only,attri))  
                   } ) 
setMethod("exonInfo", c("obj" = "PatCNVData"),
          function(obj, gene.name=NULL, capture.only=FALSE, attri=NULL)
          {    return(.getExonInfoData(obj,gene.name,capture.only,attri))
                   } ) 


#========= sampleInfo(.)

.getSampleInfo <- function(obj,sample.name,attrib)
{
  if(is.null(sample.name))
  { sel.sample.idx <- seq(1,sampleNum(obj)) } else {
    sel.sample.idx <- which(is.element(obj@sample.info$sample.name,sample.name))
  }
  return(obj@sample.info[sel.sample.idx,attrib])
}

.getSampleInfoData <- function(obj,sample.name,attrib)
{
  if(is.null(sample.name))
  { sel.sample.idx <- seq(1,sampleNum(obj)) } else {
    sel.sample.idx <- which(is.element(colData(obj@eSet)$sample.name,sample.name))
  }
  return(colData(obj@eSet)[sel.sample.idx,attrib])
}

setMethod("sampleInfo", c("obj" = "PatCNVSession"),
          function(obj,sample.name=NULL,attrib=c('sample.name','sample.type'))
          { return(.getSampleInfo(obj,sample.name,attrib)) } ) 
setMethod("sampleInfo", c("obj" = "PatCNVData"),
          function(obj,sample.name=NULL,attrib=c('sample.name','sample.type'))
          { return(.getSampleInfoData(obj,sample.name,attrib)) } ) 

#========= sampleCoverage(.)

setMethod("sampleCoverage", c("obj" = "PatCNVData"),
          function(obj,sample.type=NULL)
          { 
            if(is.null(sample.type))
            { sel.sample.idx <- seq(1,sampleNum(obj)) } else {
              sel.sample.idx <- which(obj@sample.info$sample.type==sample.type)
            }
            return(obj@sample.covg.vec[sel.sample.idx]) 
          } ) 

#========= .formatExonID()

.formatExonID <- function(exon.info,type='gene+pos',gene.sep='@')
  # type = c('gene+pos','gene+start',NULL)
  # 'gene+pos':     GeneX@chrA: startB-endC
  # 'gene+start':   GeneX@chrA: startB
  # 'pos':          chrA: startB-endC    
  # NULL:           chrA: startB
{
  
  if(is.null(type)) # chr: start position only
  {
    pos.vec <- paste(exon.info$chr,":",exon.info$start,sep="")  
  } else {
    
            if(type=="pos")
            {
              pos.vec <- paste(exon.info$chr,":",
                               exon.info$start,"-",exon.info$end,sep="")  
            }
            
            if(type=='gene+pos')
            { 
              pos.vec <- paste(exon.info$chr,":",
                               exon.info$start,"-",exon.info$end,sep="")  
              pos.vec <- paste(exon.info$gene.name,pos.vec,sep=gene.sep)
            } 
            
            if(type=='gene+start')
            { 
              pos.vec <- paste(exon.info$chr,":",exon.info$start,sep="")
              pos.vec <- paste(exon.info$gene.name,pos.vec,sep=gene.sep)
            } 
            
  } # not NULL
  
  return(pos.vec)
}

#========= avgExonCoverage(.)

setMethod("avgExonCoverage", c("obj" = "PatCNVData"),
          function(obj,gene.name=NULL,type=c('RPKM'), method='median')
          { 
            if(is.null(gene.name))
            { sel.exon.idx <- seq(1,exonNum(obj)) } else {
              sel.exon.idx <- which(is.element(rowData(obj@eSet)$gene.name,gene.name))
            }
            
            if(type=='raw')
              {  tmp.mtx <- assay(obj@eSet,"covg")[sel.exon.idx,] }
            if(type=='RPKM')
              {  tmp.mtx <- assay(obj@eSet,"RPKM")[sel.exon.idx,] }
            
                  
            if(method=='median')
            {
              exon.covg.vec <- apply(tmp.mtx,1,median,na.rm=TRUE)
            }
            if(method=='mean')
            {
              exon.covg.vec <- apply(tmp.mtx,1,mean,na.rm=TRUE)
            }
            
            exon.ID.vec <- .formatExonID(exonInfo(obj))[sel.exon.idx]
            names(exon.covg.vec) <- exon.ID.vec
            return(exon.covg.vec)
          } ) 


#=============================   coverageMatrix, cnvMatrix, fdrMatrix

.getDataMatrix <- function(obj,sample.name=NULL,
               gene.name=NULL,chr=NULL,pos.range=NULL,capture.only=FALSE,
               exon.score.vec=NULL,min.exon.score=NULL,
               data.type=c('raw','RPKM','CNV'), exon.ID.format='gene+pos')
{
  if(is.null(sample.name))
  { sel.sample.idx <- seq(1,sampleNum(obj)) } else {
    sel.sample.idx <- which(is.element(colData(obj@eSet)$sample.name,sample.name))
    if(is.null(sel.sample.idx)) { stop('input sample.name cannot be located.')}
  }
  
  if(is.null(gene.name))
  { sel.exon.idx1 <- seq(1,exonNum(obj)) } else {
    sel.exon.idx1 <- which(is.element(rowData(obj@eSet)$gene.name,gene.name))
  }
  
  if(!is.null(chr))
  {
    chr.vec <- exonInfo(obj)$chr
    start.vec <- exonInfo(obj)$start
    end.vec <- exonInfo(obj)$end
    if(!is.null(pos.range))
    {
      tmp.start <- pos.range[1]
      tmp.end <- pos.range[2]
      sel.exon.idx2 <- which( chr.vec==chr &
                                start.vec>=tmp.start & end.vec<=tmp.end )
    } else {
      sel.exon.idx2 <- which( chr.vec==chr )
    }
    sel.exon.idx1 <- intersect(sel.exon.idx1,sel.exon.idx2) 
    # in list of gene.name and in chr of given position
  }
  
  if(!is.null(exon.score.vec))
  {
    sel.exon.idx2 <- which(exon.score.vec>=min.exon.score)
    sel.exon.idx1 <- intersect(sel.exon.idx1,sel.exon.idx2)
  }
  
  if(capture.only)
  {
    sel.exon.idx <- intersect(sel.exon.idx1,
                              which(exonInfo(obj)$is.captured==1))
  } else {
    sel.exon.idx <- sel.exon.idx1
  }
  
  #print(sel.exon.idx)
  #print(sel.sample.idx)
  
  tmp.mtx <- 
    data.matrix(assay(obj@eSet,data.type)[sel.exon.idx,sel.sample.idx])
  
  exon.ID.vec <- .formatExonID(exon.info=exonInfo(obj),type=exon.ID.format)[sel.exon.idx]
  #print(exon.ID.vec)
  
  rownames(tmp.mtx) <- exon.ID.vec
  colnames(tmp.mtx) <- as.character(sampleInfo(obj,attrib='sample.name'))[sel.sample.idx]
    
  return(tmp.mtx)
}


setMethod("coverageMatrix", c("obj" = "PatCNVData"),
          function(obj,sample.name=NULL,
                   gene.name=NULL,chr=NULL,pos.range=NULL,
                   exon.score.vec=NULL,min.exon.score=NULL, capture.only=FALSE,
                   type=c('RPKM'), exon.ID.format='gene+pos') 
          { 
            if(type=='RPKM')
            {
              tmp.datatype <- 'RPKM'
            }
            if(type=='raw')
            {
              tmp.datatype <- 'covg'
            }
            return(.getDataMatrix(obj,sample.name=sample.name,
                          gene.name=gene.name,chr=chr,pos.range=pos.range,
                          exon.score.vec=exon.score.vec,min.exon.score=min.exon.score,capture.only=capture.only,
                          data.type=tmp.datatype, exon.ID.format=exon.ID.format))
          } ) 



setMethod("cnvMatrix", c("obj" = "PatCNVData"),
          function(obj,sample.name=NULL,
                   gene.name=NULL,chr=NULL,pos.range=NULL,
                   exon.score.vec=NULL,min.exon.score=NULL,capture.only=FALSE,
                   exon.ID.format='gene+pos') 
          { 
            if(!is.element("CNV",obj@datatype.vec))
            {
              stop('\n CNV matrix has not been computed yet. Please run computeMultiCNV(.) ')
            } else {
              return(.getDataMatrix(obj,sample.name=sample.name,
                                    gene.name=gene.name,chr=chr,pos.range=pos.range,
                                    exon.score.vec=exon.score.vec,min.exon.score=min.exon.score,capture.only=capture.only,
                                    data.type='CNV', exon.ID.format=exon.ID.format))
            }
          } ) 



setMethod("fdrMatrix", c("obj" = "PatCNVData"),
          function(obj,sample.name=NULL,
                   gene.name=NULL,chr=NULL,pos.range=NULL,
                   exon.score.vec=NULL,min.exon.score=NULL,capture.only=FALSE,
                   type='FDR',exon.ID.format='gene+pos') 
          { 
            if(!is.element("FDR",obj@datatype.vec))
            {
              stop('\n FDR matrix has not been computed yet. Please run estimateFDR(.) ')
            } else {
              
              if(type=='FDR')
                {
                    return(.getDataMatrix(obj,sample.name=sample.name,
                                    gene.name=gene.name,chr=chr,pos.range=pos.range,
                                    exon.score.vec=exon.score.vec,min.exon.score=min.exon.score,capture.only=capture.only,
                                    data.type='FDR', exon.ID.format=exon.ID.format))
                } 
              
              if(type=='pval')
              {
                return(.getDataMatrix(obj,sample.name=sample.name,
                                      gene.name=gene.name,chr=chr,pos.range=pos.range,
                                      exon.score.vec=exon.score.vec,min.exon.score=min.exon.score,capture.only=capture.only,
                                      data.type='pval', exon.ID.format=exon.ID.format))
              } 
              
            }
          } ) 

