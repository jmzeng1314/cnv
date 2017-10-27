a=read.table('wgs/GC_stat.10k.txt',fill = T)
a=na.omit(a)
a$GC = a[,4]/a[,3]
a$depth = a[,5]/a[,3]
a = a[a$depth<100,]
plot(a$GC,a$depth)

a=a[,c(1,2,6,5)]
colnames(a)=c('chrom','win.start','reads.gc','counts')
a$counts=floor(a$counts/150)
a$reads.mapq=30
library(seqCNA)  
head(a)
dim(a)
tail(a)

a$win.start=a$win.start*10000
a=a[a$chrom %in% paste0('chr',c(1:22,'X','Y')),]
a=a[a$win.start>0,]
a=a[a$counts>0,]
a=a[a$reads.gc>0,]
a=a[,c(1:3,5,4)]

## 200Kb windows to calculate the GC content and counts.
rco = readSeqsumm(tumour.data=a) 
rco = applyFilters(rco, trim.filter=1, mapq.filter=2)
rco = runSeqnorm(rco) 
rco = runGLAD(rco) 
plotCNProfile(rco) 
rco = applyThresholds(rco, seq(-0.8,4,by=0.8), 1)
plotCNProfile(rco)
summary(rco) 
head(rco@output)
writeCNProfile(rco,'./')


if(F){
  
  library(ggplot2)
  # GET EQUATION AND R-SQUARED AS STRING
  # SOURCE: http://goo.gl/K4yh
  lm_eqn <- function(x,y){
    m <- lm(y ~ x);
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,
                     list(a = format(coef(m)[1], digits = 2),
                          b = format(coef(m)[2], digits = 2),
                          r2 = format(summary(m)$r.squared, digits = 3)))
    as.character(as.expression(eq));
  }
  p=ggplot(a,aes(GC,depth)) + geom_point() +
    geom_smooth(method='lm',formula=y~x)+
    geom_text(x = 0.5, y = 100, label = lm_eqn(a$GC , a$depth), parse = TRUE)
  p=p+theme_set(theme_set(theme_bw(base_size=20)))
  p=p+theme(text=element_text(face='bold'),
            axis.text.x=element_text(angle=30,hjust=1,size =15),
            plot.title = element_text(hjust = 0.5) ,
            panel.grid = element_blank() 
  )
  pdf('GC_content_depth_10K.pdf')
  print(p)
  dev.off()
}













