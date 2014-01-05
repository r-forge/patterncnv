fitNullModel <- function(refdata.name,
                            type='autosome',is.plot=FALSE,capture.only=TRUE)

{

  if(class(refdata.name)!='PatCNVData')
  {
    stop('input refdata.name should be class of "PatCNVData" \n')
  }  
  
  if(!is.element("CNV",refdata.name@datatype.vec))
  {
    stop('CNV results have not been computed yet \n')
  }
  
  cnv.mtx <- cnvMatrix(refdata.name,capture.only=capture.only,exon.ID.format='pos')	  

  if(type=='all')
  {
    res <- list()
    res$type='all'
    v1 <- as.vector(cnv.mtx)
    
    if(is.plot) {plot(density(v1), xlim=c(-2,2),main='all',lwd=2,xlab='',ylab='')}
    res$location <- median(v1)
    res$scale <- mad(v1)
    res$SD <- sd(v1)
    return(res)
  }
  
  if(type=='autosome')
  {
  res <- list()
  res$type='autosome'
  res$Chr <- paste('chr',seq(1,22),sep='')
  res$location <- mat.or.vec(22,1)
  res$scale <- mat.or.vec(22,1) + 1
  res$SD <- mat.or.vec(22,1) + 1
  
  tmp.chr.pos <- unlist(strsplit(rownames(cnv.mtx),':'))
  tmp.chr.vec <- tmp.chr.pos[seq(1,2*exonNum(refdata.name),2)]
  tmp.pos.vec <- as.numeric(tmp.chr.pos[seq(2,2*exonNum(refdata.name),2)])
  
  #====== begin plotting
  if(is.plot)
  {
    par(mfrow = c(4, 6))
    
    for ( k in 1:22)
    {
      tmp_chr <- res$Chr[k]
      plot(density(as.vector(cnv.mtx[which(tmp.chr.vec==tmp_chr),])),
           xlim=c(-1,1),ylim=c(0,5),main=tmp_chr,lwd=2,xlab='',ylab='')
    }  
    par(mfrow = c(1, 1))
  } #===== end plotting
  
  
  for ( k in 1:22)
  {
    tmp_chr <- res$Chr[k]
    v1 <- as.vector(cnv.mtx[which(tmp.chr.vec==tmp_chr),])
    res$location[k] <- median(v1)
    res$scale[k] <- mad(v1)
    res$SD[k] <- sd(v1)
  }

  return(res)
  
  }
  
}


