rm(list=ls())
library(DNAcopy)
a=read.table('varScan/varScan.T4.cnv.called',stringsAsFactors = F,header = T) 
head(a)
CNA.object <- CNA(cbind( a$adjusted_log_ratio),
                  a$chrom,a$chr_start,
                  data.type="logratio",sampleid="T4")

smoothed.CNA.object <- smooth.CNA(CNA.object)

segment.smoothed.CNA.object <- segment(smoothed.CNA.object, verbose=1)

# plot(segment.smoothed.CNA.object, plot.type="w")
# 
# pdf('s.pdf',height = 20,width = 20)
# plot(segment.smoothed.CNA.object, plot.type="s") 
# dev.off()
# 
# pdf('p.pdf',height = 20,width = 20)
# plot(segment.smoothed.CNA.object, plot.type="p")
# dev.off()

sdundo.CNA.object <- segment(smoothed.CNA.object, 
                             undo.splits="sdundo", 
                             undo.SD=3,verbose=1)
# pdf('s2.pdf',height = 20,width = 20)
# plot(sdundo.CNA.object,plot.type="s")
# dev.off()

head(sdundo.CNA.object$output)
dim(sdundo.CNA.object$output)
dim(sdundo.CNA.object$data)
dim(sdundo.CNA.object$segRows)
head(sdundo.CNA.object$segRows)
tail(sdundo.CNA.object$segRows)

save.image(file = 'DNAcopy_T4.Rdata')





