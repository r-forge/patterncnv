segmentCNV <- function( session.name, data.name,
                               output_suffix='CNV_seg',
                               ref.avg.method='median', min.ref.avgRPKM=2,
                               capture.only=TRUE,
                               is.plot=TRUE, ylim=c(-3,3),cex=0.6,...)
# ... for wrapping segment() in DNAcopy package  
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
  
  #require(DNAcopy)
  
  sample_ID_vec <- as.character(sampleInfo(data.name,attrib='sample.name'))
  N_sample <- sampleNum(data.name)
  N_exon <- exonNum(data.name)
  
  tmp_DIR <- session.name@output.DIR$txt_output_DIR
  tmp_plot_DIR <- session.name@output.DIR$plot_output_DIR
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

  tmp.exoninfo <- exonInfo(data.name)
  
  sel.exon.idx <- which(ref.avg.RPKM >= min.ref.avgRPKM)
  if(capture.only==TRUE)
    { 
      non_cap_idx <- which(tmp.exoninfo$is.captured==0)
      sel.exon.idx <- setdiff(sel.exon.idx,non_cap_idx)  
      sel.exon.idx <- setdiff(sel.exon.idx,
                        which(tmp.exoninfo$chr=='chrX' | 
                              tmp.exoninfo$chr=='chrY')) # exclude sex Chr
  
      chrom <- tmp.exoninfo$chr[sel.exon.idx]
      maploc <- tmp.exoninfo$start[sel.exon.idx]
      maploc.end <-  tmp.exoninfo$end[sel.exon.idx]
  
      
  for(smpl_idx in 1:N_sample)
  {
    
    sel_sample_ID <- sample_ID_vec[smpl_idx]
    tmp_output_bed_file <- paste(tmp_DIR,sel_sample_ID,'.bed',sep='')
    tmp_output_txt_file <- paste(tmp_DIR,sel_sample_ID,'.txt',sep='')
    
    tmp.CNVvec <- as.vector(cnvMatrix(data.name,sample.name=sel_sample_ID))
    
    cna.res_My <- DNAcopy::CNA(tmp.CNVvec[sel.exon.idx], chrom, maploc,
                      data.type=c("logratio"),sampleid=sel_sample_ID)
    seg.res_My <- DNAcopy::segment(cna.res_My, verbose=0,...)
    
    seg_mtx <- seg.res_My$output
    #bed_mtx <- cbind(as.character(seg_mtx$chrom),seg_mtx$loc.start,seg_mtx$loc.end,seg_mtx$seg.mean)
    bed_mtx <- cbind(as.character(seg_mtx$chrom),seg_mtx$loc.start,
                     maploc.end[seg.res_My$segRows$endRow],seg_mtx$seg.mean) # end of segment is modified according to end of exon
    write.table(x=bed_mtx,file=tmp_output_bed_file,quote=FALSE,
                row.names=FALSE,col.names=FALSE,sep='\t')
    
    txt_mtx <- cbind(as.character(seg_mtx$chrom),seg_mtx$loc.start,
                     maploc.end[seg.res_My$segRows$endRow],seg_mtx$seg.mean,seg_mtx$num.mark) # end of segment is modified according to end of exon
    write.table(x=txt_mtx,file=tmp_output_txt_file,quote=FALSE,
                row.names=FALSE,col.names=FALSE,sep='\t')
    #the only difference between .bed and .txt files is "seg_mtx$num.mark"
    
    if (is.plot)
    {
      .patCNV.create.DIR(DIR_name=tmp_plot_DIR)
      CNV_seg_pdf_filename <- paste(tmp_plot_DIR,sel_sample_ID,output_suffix,'.pdf',sep='')
      pdf(CNV_seg_pdf_filename)
      suppressWarnings(
            DNAcopy::plotSample(seg.res_My,ylim=ylim,cex=cex,plot.type="chrombysample",xmaploc=FALSE,...)
            #xmaploc=TRUE x-axis should be according to chr position rather than bin number
        )
      dev.off()  
    } # if (is.plot)
    
  } # end of for(smpl_idx
 } # end of CNV_segment function
}  