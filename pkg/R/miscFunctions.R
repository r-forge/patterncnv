
#==================== print table with additional formatting
.patCNV.print.table <- function(tmp,N.tab=1,width=12)
{
  tmp.table <- table(tmp)
  tmp.name <- format(names(tmp.table),width=width)
  tmp.val <- format(as.character(tmp.table),width=width)
  tmp.indent <- paste(rep("\t",N.tab),collapse="")
  
  cat(tmp.indent,tmp.name,"\n",sep="\t")
  cat(tmp.indent,tmp.val,"\n",sep="\t")
}

#================== add '/' for DIR string
.patCNV.DIR.str <- function( DIR_name )
{
  tmp_str <- ""
  if (substr(DIR_name,nchar(DIR_name),nchar(DIR_name))=="/")  
   { tmp_str <- DIR_name}  else {
     tmp_str <- paste(substr(DIR_name,1,nchar(DIR_name)),"/",sep="")
   }
  return(tmp_str)
}


#================= create DIR if not existing
.patCNV.create.DIR <- function( DIR_name )
{
 if (!file.exists(DIR_name) ){
    #  cat('\"',DIR_name,'\" does not exist, creating directory...\n',sep='')
      dir.create(file.path(DIR_name),showWarnings=FALSE)
    }
}

#================ check if file exists in local path or URL
.patCNV.file.exists <- function(file_name)
{
  tmp.flag <- file.exists(file_name)
  if(tmp.flag) # file exists in local DIR
  { return(tmp.flag) } else {
    return(url.exists(file_name))  # checking file.exist as a remote file
  }
}

#===============  p-value of laplacian distribution

.patCNV.lap.pval <- function (q, location = 0, scale = 1) 
{
  if (!is.numeric(scale) | scale <=0) 
    { stop(" input scale must be positive number") }
  z.score <- (q - location)/scale
  L <- max(length(q), length(location), length(scale))
  q <- rep(q, length.out = L)
  location <- rep(location, length.out = L)
  scale <- rep(scale, length.out = L)
  p_tmp <- ifelse(q < location, 0.5 * exp(z.score), 1 - 0.5 * exp(-z.score))
  2*pmin(1-p_tmp,p_tmp)
}

 


