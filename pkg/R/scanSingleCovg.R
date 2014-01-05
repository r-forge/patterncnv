scanSingleCovg <- function(wig.filename,sample.ID,exon.bin.vec,is.capture.vec,
                               is.plot=FALSE,log.for.plot=TRUE,plot.output.DIR,
                               bin.size=10, ylab='# of exons',xlab='log2(coverage per base)')
{
  #========= It does not expect users to directly use this function
  
  r_int <- suppressWarnings(as.integer(readLines(wig.filename)))
  # suppressWarnings(.) suppress the warning of "NAs introduced by coercion", which is casued integer conversion of wig text header	  

  N_exons <- length(exon.bin.vec)
  exon_header_vec <- which(is.na(r_int))
  
  if(length(exon_header_vec)!=N_exons)  
    # checking if provided WIG file contains expected number of exons
  {
    cat('Number of exons differs between exon key file and wig file',wig.filename,'\n')
    cat('Expected # of exons from key file:',N_exons,'\n')
    cat('Observed # of exons in wig file:',length(exon_header_vec),'\n')
    cat('Please check and re-run BAM2WIG conversion if needed. \n')
    stop('Number of exons mismatch \n')
  }
  
  exon_start_vec <- exon_header_vec + 1
  exon_end_vec <- c(exon_header_vec[2:N_exons]-1,length(r_int)) # orginal
  
    
  exon_vec <- mat.or.vec(N_exons,1)
  for (k in 1:N_exons)
  {
     exon_vec[k] <- sum(r_int[exon_start_vec[k]:exon_end_vec[k]])
  }
  
    
  #============ plot exon coverage distributions
  if(is.plot)
  {
   
    .patCNV.create.DIR(plot.output.DIR)
    captured_idx <- which(is.capture.vec==1)
    non_captured_idx <- which(is.capture.vec==0)
    
    
    pdf_filename <- paste(plot.output.DIR,sample.ID,'_exon_covg.pdf',sep='')
    pdf(pdf_filename)
    
    if(log.for.plot)
    {
      hist_res <- hist(log2(exon_vec/(bin.size*exon.bin.vec)), 50,
                       main= paste('All exons ',sample.ID,sep='@'), ylab=ylab,xlab=xlab )
      hist_res <- hist(log2(exon_vec[captured_idx]/(bin.size*exon.bin.vec[captured_idx])), 50,
                       main= paste('Captured exons ',sample.ID,sep='@'), ylab=ylab, xlab=xlab )
      if(length(non_captured_idx)>0)
      {
        hist_res <- hist(log2(exon_vec[non_captured_idx]/(bin.size*exon.bin.vec[non_captured_idx])), 50,
                         main= paste('Non-captured exons ',sample.ID,sep='@'), ylab=ylab, xlab=xlab )  
      }
      
    } else {
      
      hist_res <- hist((exon_vec/(bin.size*exon.bin.vec)), 50,
                       main= paste('All exons ',sample.ID,sep='@'), ylab=ylab,xlab=xlab )
      hist_res <- hist((exon_vec[captured_idx]/(bin.size*exon.bin.vec[captured_idx])), 50,
                       main= paste('Captured exons ',sample.ID,sep='@'), ylab=ylab, xlab=xlab )
      if(length(non_captured_idx)>0)
      {
        hist_res <- hist((exon_vec[non_captured_idx]/(bin.size*exon.bin.vec[non_captured_idx])), 50,
                       main= paste('Non-captured exons ',sample.ID,sep='@'), ylab=ylab, xlab=xlab )
      }
    }
    dev.off()
  } # end of plot
  
  #============= return exon counts
  return(exon_vec)
  
}

