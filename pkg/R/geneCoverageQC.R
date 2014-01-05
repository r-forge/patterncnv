geneCoverageQC <- function(data.name,sel.gene.name,
                                  legend.layout="topleft",legend.cex=0.5)
# if legend.layout=='none', not plotting legend
{
    
  if(class(data.name)!='PatCNVData')
  {
    stop('input data.name should be class of "PatCNVData" \n')
  }  
    
  N_sample <- sampleNum(data.name)
  sel.rpkm.mtx <- coverageMatrix(obj=data.name,gene.name=sel.gene.name,
                                 exon.ID.format="gene+pos")
  N_sel_exon <- nrow(sel.rpkm.mtx)
   
  samplecolor_vec <- rainbow(N_sample)
  comb_exon_names <- rownames(sel.rpkm.mtx)
  sample_name_vec <- colnames(sel.rpkm.mtx)
  
	par(oma = c(12, 0, 0, 0))
 	plot(-1,-1,xlim=c(1,N_sel_exon),ylim=c(0,1.6*max(sel.rpkm.mtx)),
		xlab='',ylab='RPKM per exon',xaxt="none")
	axis(1, at=seq(1,N_sel_exon),
	labels=comb_exon_names, las = 2)
	for(j in 1:N_sample)
	{
	  lines(sel.rpkm.mtx[,j],type='b',col=samplecolor_vec[j],pch=j%%26,lwd=1.5,cex=0.8)  
	}
	grid()
  if(legend.layout!='none') 
  {
    legend(x=legend.layout,sample_name_vec,
           cex=legend.cex,col=samplecolor_vec,pch=seq(1,N_sample)%%26,lty=1,lwd=1)  
  }
	
	
	par(oma = c(0, 0, 0, 0)) # restore
} # end of QC function
  