exportCNVtable <- function(session.name, data.name,
                                  output_suffix='_CNV_table.txt',
				                          ref.avg.method='median', min.ref.avgRPKM=2,
                                  capture.only=TRUE,FDR.threshold=0.05,logR.cut=0.1,
                                  is.verbose=TRUE)
# ref.avg.method = 'mean' or 'median'
{
  if(class(session.name)!='PatCNVSession')
  {
    stop('input session.name should be class of "PatCNVSession" \n')
  }
  
  if(class(data.name)!='PatCNVData')
  {
    stop('input data.name should be class of "PatCNVData" \n')
  }  
  
  if(!is.element("CNV",data.name@datatype.vec))
  {
    stop('CNV results have not been computed yet \n')
  }
  
  sample_ID_vec <- as.character(sampleInfo(data.name,attrib='sample.name'))
  N_sample <- sampleNum(data.name)
  N_exon <- exonNum(data.name)
  
  tmp_DIR <- session.name@output.DIR$txt_output_DIR
  .patCNV.create.DIR(DIR_name=tmp_DIR)
    
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
  
  FDR.type <- NULL
  if(!is.null(data.name@method.info$FDR.type))  { FDR.type <- data.name@method.info$FDR.type }
  tmp.exoninfo <- exonInfo(data.name)
  
  for(smpl_idx in 1:N_sample)
  {
    
    sel_sample_ID <- sample_ID_vec[smpl_idx]
    tmp_output_file <- paste(tmp_DIR,sel_sample_ID,output_suffix,sep='')
    
    if(is.verbose)
    { cat('Exporting',smpl_idx,'-th sample',sel_sample_ID,' CNV table:',tmp_output_file,'\n') }
    
    tmp.CNVvec <- as.vector(cnvMatrix(data.name,sample.name=sel_sample_ID))
    tmp.rawcovg <- as.vector(coverageMatrix(data.name,sample.name=sel_sample_ID,type='raw'))
    tmp.RPKMcovg <- as.vector(coverageMatrix(data.name,sample.name=sel_sample_ID,type='RPKM'))
    if(is.null(FDR.type)) 
    { tmp.FDRvec <- rep(0,N_exon) } else {
      tmp.FDRvec <- as.vector(fdrMatrix(data.name,sample.name=sel_sample_ID,type='FDR'))
      tmp.pvalvec <- as.vector(fdrMatrix(data.name,sample.name=sel_sample_ID,type='pval'))
    }
    sel.exon.idx <- which(tmp.FDRvec <= FDR.threshold &
                            abs(tmp.CNVvec) >= logR.cut)
    #cat(1,length(sel.exon.idx),'\n')
    
    if(capture.only)
    { non_cap_idx <- which(tmp.exoninfo$is.captured==0)
      sel.exon.idx <- setdiff(sel.exon.idx,non_cap_idx) }
    #cat(2,length(sel.exon.idx),'\n')
    
    efficient_idx <- which(ref.avg.RPKM >= min.ref.avgRPKM)
    sel.exon.idx <- intersect(sel.exon.idx,efficient_idx)
    #cat(3,length(sel.exon.idx),'\n')
    
    outmtx1 <- tmp.exoninfo[sel.exon.idx,]
    
    sel.CNV <- tmp.CNVvec[sel.exon.idx]
    sel.raw <- tmp.rawcovg[sel.exon.idx]
    sel.RPKM<- tmp.RPKMcovg[sel.exon.idx]
    sel.baseline <- ref.avg.RPKM[sel.exon.idx]
    
    outmtx2 <- cbind(sel.CNV,sel.raw,sel.RPKM,sel.baseline)
    colnames(outmtx2) <- c('CNV.log2R','BP.covg',
                           'RPKM.covg',paste('Ref',ref.avg.method,'RPKM',sep='.'))
    
    if(is.null(FDR.type))
      {
        out.mtx <- cbind(outmtx1,outmtx2)
      } else {
        sel.FDR <- tmp.FDRvec[sel.exon.idx]
        sel.pval<- tmp.pvalvec[sel.exon.idx]
        FDRoutmtx <- cbind(sel.pval,sel.FDR) 
        colnames(FDRoutmtx) <- c('p.value',FDR.type)
        out.mtx <- cbind(outmtx1,outmtx2,FDRoutmtx)
      }
    
    write.table(x=out.mtx,file=tmp_output_file,
                quote=FALSE,sep='\t',row.names=FALSE)
  }
  
}

