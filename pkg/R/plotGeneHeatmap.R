plotGeneHeatmap <- function( data.name, sel.gene.name,
                              capture.only=TRUE, min.ref.avgRPKM=3,ref.avg.method='median',
                              logR.cut=1.5, font.cex=0.6, heatmap.margin=c(5,15),...)
#===== ... is for plot()
{
    #require(gplots)  
    
    if(class(data.name)!='PatCNVData')
    {
      stop('input data.name should be class of "PatCNVData" \n')
    }  
  
    if(!is.element("CNV",data.name@datatype.vec))
    {
      stop('CNV results have not been computed yet \n')
    }
    
    baseline.median.file <- data.name@method.info$median.covg.file
    baseline.mean.file <- data.name@method.info$mean.covg.file
    
   if(ref.avg.method=='median')	
     {ref.avg.RPKM <- unlist(read.delim(baseline.median.file,header=TRUE)$V2)}

   if(ref.avg.method=='mean')	
     {ref.avg.RPKM <- unlist(read.delim(baseline.mean.file,header=TRUE)$V2)} 

    tmp.CNV <- cnvMatrix(obj=data.name,gene.name=sel.gene.name,exon.ID.format="gene+pos",
                         exon.score.vec=ref.avg.RPKM,min.exon.score=min.ref.avgRPKM,
                         capture.only=capture.only)  
    
    heatmap.2(
              x=tmp.CNV,density.info='none',trace='none',
              col=bluered(35),symbreaks=TRUE,
              margins=heatmap.margin, keysize=1,
	            cexRow=font.cex,cexCol=font.cex,	
              breaks=seq(-logR.cut,logR.cut,length.out=36),...)
}
  