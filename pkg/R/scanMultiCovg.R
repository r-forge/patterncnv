scanMultiCovg <- function(session.name, sample.type=NULL, bin.size=10, is.verbose=TRUE, is.plot=FALSE)
  
{
  if(class(session.name)!='PatCNVSession')
  {
    stop('input session.name should be class of "PatCNVSession" \n')
  }
  

  #=========== loading configuration information
  plot_output_DIR <- session.name@output.DIR$plot_output_DIR
  txt_output_DIR <- session.name@output.DIR$txt_output_DIR
    
  exon_bin_vec <- as.numeric(values(session.name@exon.info)[,'N.exon.bin'])
  is_capture_vec <- as.numeric(values(session.name@exon.info)[,'is.captured'])
    
  
  if(is.null(sample.type)) # using all the samples
  {
    wig_filename_vec <- session.name@sample.info$file.name
    sample_ID_vec <- session.name@sample.info$sample.name  
    sel_sample_idx <- seq(1,length(sample_ID_vec))
  } else {                 
    # using selected samples with given sample type
    sel_sample_idx <- which(session.name@sample.info$sample.type==sample.type)
    wig_filename_vec <- session.name@sample.info$file.name[sel_sample_idx]
    sample_ID_vec <- session.name@sample.info$sample.name[sel_sample_idx]
  }
  
  
  N_exons <- length(exon_bin_vec)
  N_samples <- length(wig_filename_vec)
  
  tmp_obj <- new("PatCNVData")
  
  tmp_obj@N.exon <- N_exons
  tmp_obj@N.sample <- N_samples
  
  #tmp_obj@exon.info <- session.name@exon.info
  
 
  
  tmp_obj@eSet@rowData <- session.name@exon.info
  tmp_obj@eSet@colData <- DataFrame(list(sample.name=session.name@sample.info$sample.name[sel_sample_idx],
                              sample.type=session.name@sample.info$sample.type[sel_sample_idx],
                              subject.ID=session.name@sample.info$subject.ID[sel_sample_idx]))
  
  
  
  tmp.zeromtx <- mat.or.vec(N_exons,N_samples)
  assays(tmp_obj@eSet) <- SimpleList(covg=tmp.zeromtx,RPKM=tmp.zeromtx,
             CNV=tmp.zeromtx,pval=tmp.zeromtx,FDR=tmp.zeromtx)
  
  #=========== end of loading configuration information
  
  total_count_vec <- mat.or.vec(N_samples,1)
  names(total_count_vec) <- sample_ID_vec
  exon_count_mtx <- mat.or.vec(N_exons,N_samples)
  colnames(exon_count_mtx) <- sample_ID_vec
  exon_RPKM_mtx <- mat.or.vec(N_exons,N_samples)
  colnames(exon_RPKM_mtx) <- sample_ID_vec

  
  for(k in 1:N_samples) 
  {
    tmp_wig_filename <- wig_filename_vec[k]
    tmp_sample_ID <- sample_ID_vec[k]
    if (is.verbose)
    {
      cat(paste('scanning ',k,'-th sample:',tmp_sample_ID,'\n'))  
    }
    count_vec <- 
      scanSingleCovg(wig.filename=tmp_wig_filename,sample.ID=tmp_sample_ID,
                              exon.bin.vec=exon_bin_vec,is.capture.vec=is_capture_vec,
                              bin.size=bin.size,is.plot=is.plot,plot.output.DIR=plot_output_DIR)
                                            
    exon_count_mtx[,k] <- count_vec
    total_count_vec[k] <- sum(count_vec,na.rm=TRUE) # total bp counts
    exon_RPKM_mtx[,k] <- (count_vec*1e9)/(bin.size*exon_bin_vec*total_count_vec[k])
      # RPKM = count_in_exon/ ( (kb) * (total_counts/1e6) )
      #      = count_in_exon/ ( (N_bins * bin.size/1e3) * total_counts/1e6)
      #      = count_in_exon / (N_bins * bin.size * total_counts/1e9)
      #      = count_in_exon*1e9 / (N_bins * bin.size * total_counts)
  }
  
  #============= return exon counts
  
  tmp_obj@sample.covg.vec <- total_count_vec
  assay(tmp_obj@eSet,"covg") <- exon_count_mtx
  assay(tmp_obj@eSet,"RPKM") <- exon_RPKM_mtx
  
  tmp_obj@datatype.vec <- c("covg","RPKM")
  
  return(tmp_obj)
}



