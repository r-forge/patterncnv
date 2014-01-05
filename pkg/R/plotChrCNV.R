plotChrCNV <- function( data.name, sample.name, 
				sel.chr='chr1', pos.range=NULL,
				capture.only=TRUE, min.ref.avgRPKM=2,ref.avg.method='median',
				cex=0.65,col='steelblue3',ylim=c(-3,3),
				xlab='Mb',ylab='log2-ratio',
        ideogram=FALSE,ideogram.cex=0.75,...)
#===== 
{
	
  if(class(data.name)!='PatCNVData')
  {
    stop('input data.name should be class of "PatCNVData" \n')
  }  
  
  
	sel.sample.idx <- which(colData(data.name@eSet)$sample.name==sample.name)
	N_exon <- exonNum(data.name)
	if(!length(sel.sample.idx))
	 { stop(paste(sample.name,'cannot be located in input','\n')) }	
  
    baseline.median.file <- data.name@method.info$median.covg.file
    baseline.mean.file <- data.name@method.info$mean.covg.file
  
    	if(!file.exists(baseline.median.file) & 
           !(file.exists(baseline.mean.file))) {
    	  # no coverage files, ignoring this information
    	  ref.avg.RPKM <- mat.or.vec(N_exon,1) + 1e10
    	  
    	} else {
    	  if(ref.avg.method=='median')  
    	  {ref.avg.RPKM <- unlist(read.delim(baseline.median.file,header=FALSE,sep='\t')$V2)}
    	  
    	  if(ref.avg.method=='mean')	
    	  {ref.avg.RPKM <- unlist(read.delim(baseline.mean.file,header=FALSE,sep='\t')$V2)} 
    	}
    
      tmp.CNV <- cnvMatrix(obj=data.name,sample.name=sample.name,exon.ID.format=NULL,
                           chr=sel.chr,pos.range=pos.range,
                           exon.score.vec=ref.avg.RPKM,min.exon.score=min.ref.avgRPKM,
                           capture.only=capture.only)    
      gnm_pos <- 
          as.numeric(gsub(pattern=paste(sel.chr,':',sep=''),replacement='',x=rownames(tmp.CNV)))/1e6
      plot(gnm_pos, tmp.CNV, type='p',
		          main=paste(sel.chr,' log2-ratio CNV @ ',sample.name), 
		          cex=cex,col=col,xlab=xlab,ylab=ylab,ylim=ylim,...)  
	  axis(side=4)  
  
   if(ideogram==TRUE)
   {
     par(cex=ideogram.cex)
     #require(SNPchip)
     chr_j <- as.numeric(unlist(strsplit(tolower(sel.chr),'chr'))[2])
     plotIdiogram(chr_j,build='hg19',ylim=ylim,cytoband.ycoords=c(ylim[1]+0.1, ylim[1]+0.3),
                  new=FALSE,label.y=ylim[1]+0.6,unit='Mb',cex.axis=0.45)
     par(cex=1)   
   }
     
  
  
}
