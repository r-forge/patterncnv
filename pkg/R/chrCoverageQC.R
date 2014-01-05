chrCoverageQC <- function(  data.name,
                                   legend.layout="topright",plot.cex=2,legend.cex=1)
{
  
    if(class(data.name)!='PatCNVData')
    {
      stop('input data.name should be class of "PatCNVData" \n')
    }  
    
   N_sample <- sampleNum(data.name)
   chr.rpkm.mtx <- mat.or.vec(22,N_sample)

   for(k in 1:22)
   {
	  tmp.chr <- paste('chr',k,sep='')
	  tmp.RPKM.mtx <- coverageMatrix(obj=data.name,chr=tmp.chr)
	  chr.rpkm.mtx[k,] <- apply(tmp.RPKM.mtx,2,sum,na.rm=TRUE)
	 }  

   samplecolor_vec <- rainbow(N_sample)
	
 	plot(-1,-1,xlim=c(1,22),ylim=c(0,max(chr.rpkm.mtx,na.rm=TRUE)*1.5),xaxp=c(1,22,21),
		xlab='Chr index',ylab='RPKM sum per Chr')
	for(j in 1:N_sample)
	{
	  lines(chr.rpkm.mtx[,j],type='b',col=samplecolor_vec[j],pch=j,lwd=1.5,cex=plot.cex)  
	}
	grid()

	  if(legend.layout!='none') 
	  {
		legend(x=legend.layout,colnames(tmp.RPKM.mtx),
	       	cex=legend.cex,col=samplecolor_vec,pch=seq(1,N_sample),lty=1,lwd=1)
	  }
 

} # end of QC function
  