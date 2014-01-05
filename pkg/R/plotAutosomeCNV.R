plotAutosomeCNV <- function( data.name, sample.name, chr.Len.vec,
                                    capture.only=TRUE, min.ref.avgRPKM=2,ref.avg.method='median',
                                    ylim=c(-3,3),cex=0.6,	xlab='Mb',ylab='log2-ratio',
                                    color.vec=c('red','blue','green','orange','brown','purple','black'),... )
{
    if(class(data.name)!='PatCNVData')
    {
      stop('input data.name should be class of "PatCNVData" \n')
    }  
    
    if(!is.element("CNV",data.name@datatype.vec))
    {
      stop('CNV results have not been computed yet \n')
    }
    
    chr_length_vec <- chr.Len.vec  
 
    sel.sample.idx <- which(sampleInfo(data.name,attrib="sample.name")==sample.name)
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
    

    chr_color.vec <- rep(color.vec,ceiling(22/length(color.vec))) # repeat color patterns
  
    chr_start_vec <- chr_length_vec*0
    chr_start_vec[1] <- 0
    for(k in 2:22)
    {
      chr_start_vec[k] <- chr_start_vec[k-1] + chr_length_vec[k-1]
    }
  
 
    plot(-1e3,-1e3,xlim=c(0,(chr_start_vec[22]+chr_length_vec[22])/1e6),
         ylim=ylim,xlab='Mb',ylab='log2-ratio',
         main=paste('Whole exome CNV of',sample.name))
    
    tmp.CNV <- cnvMatrix(obj=data.name,sample.name=sample.name,exon.ID.format=NULL,
                         exon.score.vec=ref.avg.RPKM,min.exon.score=min.ref.avgRPKM,
                         capture.only=capture.only)  
    tmp.chr.pos <- unlist(strsplit(rownames(tmp.CNV),':'))
    tmp.chr.vec <- tmp.chr.pos[seq(1,2*exonNum(data.name),2)]
    tmp.pos.vec <- as.numeric(tmp.chr.pos[seq(2,2*exonNum(data.name),2)])
    
    for(k in 1:22)
    {
      sel_chr <- paste('chr',k,sep='')
      tmp.sel.idx <- which(tmp.chr.vec==sel_chr)
      
      if(length(tmp.sel.idx)!=0) # as long as this chr is not empty
      {
        gnm_pos <- chr_start_vec[k]/1e6 + tmp.pos.vec[tmp.sel.idx]/1e6
        
        lines(gnm_pos, tmp.CNV[tmp.sel.idx], type='p',
              cex=cex,col=chr_color.vec[k] )    
      }
      
    }

  
}
