learnPatterns <- function( session.name, refdata.name,
                                  sample.type=NULL, episl=1, bin.size=10,
                                  excluded.sample.name = NULL )
{
    
  if(class(session.name)!='PatCNVSession')
  {
    stop('input session.name should be class of "PatCNVSession" \n')
  }
  
  if(class(refdata.name)!='PatCNVData')
  {
    stop('input refdata.name should be class of "PatCNVData" \n')
  }  

  tmp_DIR <- session.name@output.DIR$txt_output_DIR
  .patCNV.create.DIR(DIR_name=tmp_DIR)
  
  #require('matrixStats')
  
  #============= summarizing total counts and mean/median RPKM
  total_count_vec <- sampleCoverage(refdata.name)
  mean_RPKM_file <- session.name@Misc$mean_RPKM_file
  median_RPKM_file <- session.name@Misc$median_RPKM_file
  
  #=========== loading configuration information
  plot_output_DIR <- session.name@output.DIR$plot_output_DIR
  txt_output_DIR <- session.name@output.DIR$txt_output_DIR
  
  exon_bin_vec <- session.name@exon.info$N.exon.bin
  is_capture_vec <- as.numeric(session.name@exon.info$is.captured)
  
  if(is.null(sample.type)) # using all the samples
  {
    sample_ID_vec <- colData(refdata.name@eSet)$sample.name  

  } else {                 
    # using selected samples with given sample type
    sel_sample_idx <- which(colData(refdata.name@eSet)$sample.type==sample.type &
                     !is.element(colData(refdata.name@eSet)$sample.name, excluded.sample.name) )
    sample_ID_vec <- colData(refdata.name@eSet)$sample.name[sel_sample_idx]
  }
  mch.session.sample.idx <- match(sample_ID_vec,
                                  session.name@sample.info$sample.name)
  wig_filename_vec <- session.name@sample.info$file.name[mch.session.sample.idx]
  
  RKM_vec <- bin.size*total_count_vec[sample_ID_vec]/1e9
  
  # computing RPKM based on total bp counts of selected samples
  # RPKM = read/ ( (kb) * (total_counts/1e6) )
  #      = read/ (bin.size/1e3 * total_counts/1e6)
  #      = read/ (bin.size*total_counts/1e9)
  
  N_exons <- length(exon_bin_vec)
  N_samples <- length(wig_filename_vec)
  cat('\n', N_samples,' samples are selected for learning patterns.\n',sep='')
  cat('\n Total BP-level coverage for selected samples: \n')
  print(total_count_vec[sample_ID_vec])
  
  
  tmp.RPKM.mtx <- coverageMatrix(obj=refdata.name,sample.name=sample_ID_vec)
  wig.mean_RPKM_vec <- apply(tmp.RPKM.mtx,1,mean)
  wig.median_RPKM_vec <- apply(tmp.RPKM.mtx,1,median)
  write.table(x=wig.mean_RPKM_vec,file=mean_RPKM_file,quote=FALSE,row.names=TRUE,col.names=FALSE,sep='\t')
  write.table(x=wig.median_RPKM_vec,file=median_RPKM_file,quote=FALSE,row.names=TRUE,col.names=FALSE,sep='\t')
  
  #=========== end of loading configuration information
  
  
  #================ initialize mean and SD wig files
  mean_vec <- mat.or.vec(N_exons,1)
  SD_vec <- mat.or.vec(N_exons,1)
  
  mean_wig_file <- session.name@pattern.file$avg_wig
  SD_wig_file <- session.name@pattern.file$var_wig
  
  
  #=== 
  z <- 1
  r_str <- readLines(wig_filename_vec[z])
  #r_int <- as.integer(readLines(wig_filename_vec[z]))
  r_int <- suppressWarnings(as.integer(readLines(wig_filename_vec[z])))
  # suppress the warning of "NAs introduced by coercion", induced by converting text header of wig to integer

  N_wig_bins <- length(r_int)
  N_exons <- length(exon_bin_vec)
  exon_header_vec <- which(is.na(r_int))
  exon_start_vec <- exon_header_vec + 1
  exon_end_vec <- c(exon_header_vec[2:N_exons]-1,length(r_int))
  
  wig_bin_mtx <- mat.or.vec(N_wig_bins,N_samples)
  wig_bin_mtx[,z] <- log2( (r_int + episl)/ RKM_vec[z])
  for(z in 2:N_samples) 
  {
    #r_int <- as.integer(readLines(wig_filename_vec[z]))
    r_int <- suppressWarnings(as.integer(readLines(wig_filename_vec[z])))
    # suppress the warning of "NAs introduced by coercion", induced by converting text header of wig to integer
    wig_bin_mtx[,z] <- log2( (r_int + episl)/ RKM_vec[z])
    cat(z,'-th sample: ', names(RKM_vec)[z],'\n')
    #print(proc.time()-ptm)
  }
  
  
  
  #avg_pattern_vec <-  apply(wig_bin_mtx,1,median)      # robust estimate

  median_vec <- rowMedians(wig_bin_mtx)
  avg_pattern_vec <- format(median_vec,digits=2,trim=TRUE)

  #var_pattern_vec <-  apply(wig_bin_mtx,1,mad)*1.4826  # robust estimate of SD
  var_pattern_vec <- format(rowMedians(wig_bin_mtx)*1.4826,digits=2,trim=TRUE)

  
  
  avg_pattern_vec[exon_header_vec] <- r_str[exon_header_vec]
  var_pattern_vec[exon_header_vec] <- r_str[exon_header_vec]
  #log2( (as.numeric(readLines(con_list[[z]],1)) + episl)/
  #        RKM_vec[z] )
  
  
  #========= writing
  
  write.table(x=avg_pattern_vec,file=mean_wig_file,quote=FALSE,row.names=FALSE,col.names=FALSE)
  cat('outputing average pattern file... \n')

  write.table(x=var_pattern_vec,file=SD_wig_file,quote=FALSE,row.names=FALSE,col.names=FALSE)
  cat('outputing variability pattern file... \n')

  
  
}

