makeSimulation <- function(ini.file='simu.ini',
                                  sample.info.file='sample_info.txt',
                                  with.pattern=FALSE)
{
  extdata.path <- system.file("extdata",package='patternCNV')
  
  tmp.sample.info <- system.file("extdata", "simu_sample_info.txt", package = "patternCNV")
  sample.info.data <- 
  gsub(pattern='<wig.file.path>',replacement=extdata.path,x=readLines(tmp.sample.info))
  
  writeLines(sample.info.data,sample.info.file)
  
  if(!with.pattern)
   {  tmp.simu.ini <- system.file("extdata", "simple_simu.ini", package = "patternCNV") 
   } else {
      tmp.simu.ini <- system.file("extdata", "complete_simu.ini", package = "patternCNV") 
   }
  
  simu.data0 <- gsub(pattern='<exon.key.file>',
                     replacement=system.file("extdata", "sim_exon_key.txt", package = "patternCNV"),
                     x=readLines(tmp.simu.ini))  
  simu.data  <- gsub(pattern='<sample.info.file>',
                     replacement=sample.info.file, x=simu.data0)  
  
  if(with.pattern) # filling additional pattern files
  {
    simu.data  <- gsub(pattern='<average.pattern.file>',
                       replacement=system.file("extdata", "avg_pattern.wig", package = "patternCNV"), 
                       x=simu.data)  
    simu.data  <- gsub(pattern='<variability.pattern.file>',
                       replacement=system.file("extdata", "var_pattern.wig", package = "patternCNV"), 
                       x=simu.data)
    simu.data  <- gsub(pattern='<median.RPKM.file>',
                       replacement=system.file("extdata", "median_RPKM.txt", package = "patternCNV"), 
                       x=simu.data)
    simu.data  <- gsub(pattern='<mean.RPKM.file>',
                       replacement=system.file("extdata", "mean_RPKM.txt", package = "patternCNV"), 
                       x=simu.data)
  }
  
  writeLines(simu.data,ini.file)
}

