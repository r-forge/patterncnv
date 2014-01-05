
setGeneric("exonNum", signature = "obj", function(obj) 
          {standardGeneric("exonNum")} )

setGeneric("sampleNum", signature = "obj", function(obj) 
            {standardGeneric("sampleNum")} )

setGeneric("exonInfo", signature = "obj",
           function(obj,gene.name=NULL,capture.only=FALSE,attri=NULL) 
             {standardGeneric("exonInfo")} )

setGeneric("sampleInfo",signature='obj', 
           function(obj,sample.name=NULL,attrib=c('sample.name','sample.type'))
           {standardGeneric('sampleInfo')})

setGeneric("sampleCoverage",signature='obj', 
           function(obj,sample.type=NULL)
             {standardGeneric('sampleCoverage')})

setGeneric("avgExonCoverage",signature='obj', 
           function(obj,gene.name=NULL,type=c('RPKM'), method='median') 
           {standardGeneric('avgExonCoverage')})

setGeneric("coverageMatrix",signature='obj', 
           function(obj,sample.name=NULL,
                    gene.name=NULL,chr=NULL,pos.range=NULL,
                    exon.score.vec=NULL,min.exon.score=NULL, capture.only=FALSE,
                    type=c('RPKM'), exon.ID.format='gene+pos') 
           {standardGeneric('coverageMatrix')})

setGeneric("cnvMatrix",signature='obj', 
           function(obj,sample.name=NULL,
                    gene.name=NULL,chr=NULL,pos.range=NULL,
                    exon.score.vec=NULL,min.exon.score=NULL,capture.only=FALSE,
                    exon.ID.format='gene+pos') 
           {standardGeneric('cnvMatrix')})

setGeneric("fdrMatrix",signature='obj', 
           function(obj,sample.name=NULL,
                    gene.name=NULL,chr=NULL,pos.range=NULL,
                    exon.score.vec=NULL,min.exon.score=NULL,capture.only=FALSE,
                    type='FDR', exon.ID.format='gene+pos') 
           {standardGeneric('fdrMatrix')})

