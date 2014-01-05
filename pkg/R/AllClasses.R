setClass(Class="PatCNVSession", representation=
         representation(exon.info="GRanges",
                        sample.info="data.frame",
                        output.DIR="list",
                        pattern.file='list',
                        Misc='list',
                        N.exon='numeric',
                        N.sample='numeric')         )

setClass( Class="PatCNVData", representation=
				representation( eSet="SummarizedExperiment",
								method.info='list',
								datatype.vec='character',
								sample.covg.vec='numeric',
								N.exon='numeric',
								N.sample='numeric')         )					
