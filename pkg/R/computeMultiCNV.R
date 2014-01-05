computeMultiCNV <- function( session.name, data.name=NULL, sample.type=NULL,
				   ref.type='average.pattern',	
           episl=1, small.delta=1e-5, bin.size=10,
				   zero.median.adjust=TRUE,is.verbose=FALSE)

{
 
  if(class(session.name)!='PatCNVSession')
  {
    stop('input session.name should be class of "PatCNVSession" \n')
  }
  
  if(class(data.name)!='PatCNVData' & !is.null(data.name))
  {
    stop('input data.name should be either class of "PatCNVData" or NULL \n')
  }  
  
  if(is.null(data.name)) # create an object of "PatCNVData" class
  {
    data.name <- scanMultiCovg(session.name=session.name,sample.type=sample.type,
                                      bin.size=bin.size,is.verbose=FALSE,is.plot=FALSE)
  }
  
  #=========== loading session information
  plot_output_DIR <- session.name@output.DIR$plot_output_DIR
  txt_output_DIR <- session.name@output.DIR$txt_output_DIR
  
  exon_bin_vec <- session.name@exon.info$N.exon.bin
  is_capture_vec <- as.numeric(session.name@exon.info$is.captured)
  
  if(is.null(sample.type)) # using all the samples
  {
    sample_ID_vec <- colData(data.name@eSet)$sample.name  
    
  } else {                 
    # using selected samples with given sample type
    sel_sample_idx <- which(colData(data.name@eSet)$sample.type==sample.type)
    sample_ID_vec <- colData(data.name@eSet)$sample.name[sel_sample_idx]
  }
  mch.session.sample.idx <- match(sample_ID_vec,
                                  session.name@sample.info$sample.name)
  wig_filename_vec <- session.name@sample.info$file.name[mch.session.sample.idx]
    
  N_exon <- length(exon_bin_vec)
  N_sample <- length(wig_filename_vec)
  
  
  ref_mean_wigfile <- session.name@pattern.file$avg_wig
  ref_SD_wigfile <- session.name@pattern.file$var_wig
  
  #=========== end of loading session information
  
  CNV.mtx <- mat.or.vec(N_exon,N_sample)
  colnames(CNV.mtx) <- sample_ID_vec
  
  if(is.verbose) {
  pb <- txtProgressBar(style=3,max=N_sample) }
  for (k in 1:N_sample)
  {
    
    wig_filename <- wig_filename_vec[k]
    sample_ID <- sample_ID_vec[k]
    sel.sample.name <- sample_ID_vec[k]
 
    #cat('\n',k,wig_filename,sample_ID,'\n',sep=' | ')
    if(is.verbose) {  setTxtProgressBar(pb, k) }        
	   
    cnv_res <- 
      computeSingleCNV( session.name, sel.sample.name, 
			                         ref.type=ref.type, episl=episl, small.delta=small.delta,
                               bin.size=bin.size, is.verbose=is.verbose)
    CNV.mtx[,k] <- cnv_res$CNV
  }  
  
  if(zero.median.adjust) {

           for (k in 1:N_sample)
           {
             CNV.mtx[,k] <- CNV.mtx[,k] - median(CNV.mtx[,k],na.rm=TRUE)
            }
	}

   cat('\n')	

  
  mch.sample.idx <- match(colnames(CNV.mtx),names(data.name@sample.covg.vec))
  tmp.obj <- new('PatCNVData')
  
  
  tmp.obj@sample.covg.vec <- data.name@sample.covg.vec[mch.sample.idx]
  
  #=== directly subset selected samples
  tmp.obj@eSet <- data.name@eSet[,mch.sample.idx]
  
  #tmp.obj@covg.mtx <- data.name@covg.mtx[,mch.sample.idx]
  #tmp.obj@RPKM.mtx <- data.name@RPKM.mtx[,mch.sample.idx]
  
  assay(tmp.obj@eSet,"CNV") <- CNV.mtx
  tmp.obj@datatype.vec <- union(data.name@datatype.vec,"CNV")
  
  
  tmp.obj@method.info$ref.type <- ref.type
  tmp.obj@method.info$median.covg.file <- session.name@Misc$median_RPKM_file
  tmp.obj@method.info$mean.covg.file <- session.name@Misc$mean_RPKM_file
    
  
  tmp.obj@N.exon <- data.name@N.exon
  tmp.obj@N.sample <- length(mch.sample.idx)
  #tmp.obj@sample.info <- data.frame(colData(data.name@eSet)[mch.sample.idx,],stringsAsFactors=FALSE)
  #print(colData(data.name@eSet))


  
  return(tmp.obj)
  
}

