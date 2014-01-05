createSession <- function(config.file)
{
  tmp_obj <- new("PatCNVSession")
  if(!.patCNV.file.exists(config.file)) {
    stop(paste('pattern-CNV configuration file cannot be located:',config.file,'\n'))}
  
  plot_output_DIR <- NULL
  txt_output_DIR <- NULL
  exon_key_file <- NULL
  sample_info_file <- NULL
  avg_pattern_file <- NULL
  var_pattern_file <- NULL
  median_RPKM_file <- NULL
  mean_RPKM_file <- NULL
  
  
  source(config.file,local=TRUE) 
  #"local    TRUE, FALSE or an environment, determining where the 
  #parsed expressions are evaluated. FALSE (the default) corresponds 
  #to the user's workspace (the global environment) and TRUE to 
  #the environment from which source is called."
  
#   #== debugging
#   print(avg_pattern_file)
#   print(var_pattern_file)
#   print(median_RPKM_file)
#   print(mean_RPKM_file)
#   
  
  
  #=== DIR information
  DIR_info <- list()
  DIR_info$plot_output_DIR <- .patCNV.DIR.str(plot_output_DIR)
  DIR_info$txt_output_DIR <- .patCNV.DIR.str(txt_output_DIR)
  
  
  #=== detecting existence of important files
  if(!.patCNV.file.exists(exon_key_file))
  {
    str1 <- '\n Error! Exon key file cannot be located. \n Please check \"exon_key_file\" in \"'
    str2 <- '\" and re-run the code. \n'
    stop( paste(str1,config.file,str2,sep='') )
    return()
  }
  
  if(!.patCNV.file.exists(sample_info_file))
  {
    str1 <- '\n Error! sample information file cannot be located. \n Please check \"sample_info_file\" in \"'
    str2 <- '\" and re-run the code. \n'
    stop( paste(str1,config.file,str2,sep='') )
    return()
  }
  
  #=== load exon_key
  cat('\n Loading exon information ...\n')
  exon_info <- read.delim(exon_key_file,stringsAsFactors=FALSE)
  
  if(substr(exon_info$Chr[1],1,1)!='c'){
    exon_info$Chr <- paste('chr',exon_info$Chr,sep='')   # if not starting with 'chr'
  }
  
  
  exon_info$exon_bin_vec <- as.numeric(exon_info$Bin_Count)
  exon_info$is_capture_vec <- as.numeric(exon_info$InCapture)
  cat('\n Number of exons in capture:',length(which(exon_info$is_capture_vec==1)))
  cat('\n Number of exons not in capture:',length(which(exon_info$is_capture_vec==0)))
  cat('\n')
  
  tmp.N.exon <- length(exon_info$is_capture_vec)
  
  exon_GRange <- 
    GRanges(seqnames=exon_info$Chr,ranges=IRanges(start=exon_info$Start,end=exon_info$Stop),
            gene.name=exon_info$Genes, 
            is.captured=Rle(exon_info$is_capture_vec),
            N.exon.bin=exon_info$exon_bin_vec)
              
  #=== load sample_info
  cat('\n Loading sample information ...\n')
  file_info <- read.delim(sample_info_file,stringsAsFactors=FALSE)
  
  tmp.N.unique.sample <- length(unique(file_info$sample.name))
  if(length(file_info$sample.name)!=tmp.N.unique.sample)
  {
    stop(tmp.N.unique.sample, 'unique sample names out of ', nrow(file_info), 'input samples \n
         sample name has to be unique and non-missing. \n')
  }
  
  
  file_info$ID <- paste(file_info$sample.name,'(',file_info$sample.type,')',sep='') # unique ID = name + type
  cat('\n Number of samples:', nrow(file_info))
  cat('\n Sample type distributions:')
  print(table(file_info$sample.type))
  cat('\n')
  cat('\n Locating ',nrow(file_info),'input WIG files \n')
  
  #txtpb <- txtProgressBar(min=1,max=nrow(file_info),style=3)
  for(k in 1:nrow(file_info))
  {
    #setTxtProgressBar(txtpb, k)
    cat(k,'...')
    
    if( round(k/10)*10 == k) { cat('\n') } # enter a new line for every 10 samples
    
    tmp_file <- file_info$file.name[k]
    if(!.patCNV.file.exists(tmp_file))
    {
      stop('\n Error! Input file \"',tmp_file,'\" cannot be located \n',sep='')
    }
  }
  cat('\n Done \n')
 
  
  lst.Misc <- list()
  lst.pattern <- list()
  if(!is.null(avg_pattern_file))   
  { lst.pattern$avg_wig <- avg_pattern_file } else 
  { lst.pattern$avg_wig <- 
      paste(.patCNV.DIR.str(DIR_info$txt_output_DIR),'avg_pattern.wig',sep='') }
  
  if(!is.null(var_pattern_file)) 	
  { lst.pattern$var_wig <- var_pattern_file } else 
  { lst.pattern$var_wig <- 
      paste(.patCNV.DIR.str(DIR_info$txt_output_DIR),'var_pattern.wig',sep='') }
  
  if(!is.null(median_RPKM_file)) 	
  { lst.Misc$median_RPKM_file <- median_RPKM_file } else 
  { lst.Misc$median_RPKM_file <- 
      paste(.patCNV.DIR.str(DIR_info$txt_output_DIR),'median_RPKM.txt',sep='') }
  
  if(!is.null(mean_RPKM_file)) 	
  { lst.Misc$mean_RPKM_file <- mean_RPKM_file } else 
  { lst.Misc$mean_RPKM_file <- 
      paste(.patCNV.DIR.str(DIR_info$txt_output_DIR),'mean_RPKM.txt',sep='') }
  
  tmp_obj@exon.info <- exon_GRange
  tmp_obj@sample.info <- data.frame(file_info)
  tmp_obj@output.DIR <- DIR_info
  tmp_obj@pattern.file <- lst.pattern
  tmp_obj@Misc <- lst.Misc
  tmp_obj@N.exon <- tmp.N.exon
  tmp_obj@N.sample <- tmp.N.unique.sample
  
  return(tmp_obj)
} # end of create.PatCNVSession
