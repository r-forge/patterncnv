estimateFDR <- function(data.name, null.model, FDR.type='q.value',is.plot=FALSE )
# FDR.type= 'q.value' or 'localFDR'
{
  if(class(data.name)!='PatCNVData')
  {
    stop('input data.name should be class of "PatCNVData" \n')
  }  
  
  if(!is.element("CNV",data.name@datatype.vec))
  {
    stop('CNV results have not been computed yet \n')
  }
  
  cnv.mtx <- cnvMatrix(data.name,exon.ID.format='pos')
  N_sample <- sampleNum(data.name)
  N_exon <- exonNum(data.name)
  
  pval.mtx <- mat.or.vec(N_exon,N_sample) + 1
  FDR.mtx <- mat.or.vec(N_exon,N_sample) + 1
  
  if(null.model$type=='all')
  {
    for (k in 1:N_sample)
    {
        sel_vec <- cnv.mtx[,k]
        p_vec <- .patCNV.lap.pval(q=sel_vec,location=null.model$location,scale=null.model$scale)
        fdr_res <- fdrtool::fdrtool(x=p_vec,statistic='pvalue',plot=FALSE,verbose=FALSE)
        
        pval.mtx[,k] <- p_vec
        if(FDR.type=='q.value') {FDR.mtx[,k] <- fdr_res$qval}
        if(FDR.type=='localFDR') {FDR.mtx[,k] <- fdr_res$lfdr}
    }  
      if(is.plot)       { hist(p_vec,50,xlab='p-value',ylab='#') }
  } # end of if(null.model$type=='all
    
  if(null.model$type=='autosome')
  {
    tmp.chr.pos <- unlist(strsplit(rownames(cnv.mtx),':'))
    tmp.chr.vec <- tmp.chr.pos[seq(1,2*exonNum(data.name),2)]
    tmp.pos.vec <- as.numeric(tmp.chr.pos[seq(2,2*exonNum(data.name),2)])
    for (k in 1:N_sample)
    {
      for(j in 1:22)
      {
        sel_idx <- which(tmp.chr.vec==null.model$Chr[j])
        sel_vec <- cnv.mtx[sel_idx,k]
        p_vec <- .patCNV.lap.pval(q=sel_vec,location=null.model$location[j],scale=null.model$scale[j])
        fdr_res <- fdrtool::fdrtool(x=p_vec,statistic='pvalue',plot=FALSE,verbose=FALSE)
        
        pval.mtx[sel_idx,k] <- p_vec
        if(FDR.type=='q.value') {FDR.mtx[sel_idx,k] <- fdr_res$qval}
        if(FDR.type=='localFDR') {FDR.mtx[sel_idx,k] <- fdr_res$lfdr}
              }
      
      if(is.plot) { hist(p_vec,50,xlab='p-value',ylab='#') }
    }

  } # end of if(null.model$type=='autosome
  
  #data.name@pval.mtx <- pval.mtx
  #data.name@FDR.mtx <- FDR.mtx
  
  assay(data.name@eSet,"pval") <- pval.mtx
  assay(data.name@eSet,"FDR") <- FDR.mtx
  
  data.name@datatype.vec <- union(data.name@datatype.vec,c("pval","FDR"))
  data.name@method.info$FDR.type <- FDR.type
  
  return(data.name)
}


