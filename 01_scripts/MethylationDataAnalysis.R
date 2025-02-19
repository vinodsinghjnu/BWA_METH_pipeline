conda activate R.4.2_cloned_mapdemo_env
R
setwd('/home/vinodsingh/mutationalscanning/Workspaces/vinod/codes/bwa-meth_pipeline-master')
data.dir='/home/vinodsingh/mutationalscanning/Workspaces/vinod/codes/bwa-meth_pipeline-master/'

library(methylKit)
files1=list.files(path=paste0(data.dir,'07_methyl_dackel/Muenster/'), pattern='*_final_CpG.methylKit.gz$', full.names = T)
files2=list.files(path=paste0(data.dir,'07_methyl_dackel/PRJEB34432_NormalSperm/'), pattern='*_CpG.methylKit.gz$', full.names = T)
files=c(files1,files2)
id1=lapply(sapply(files1,basename), function(x) paste(strsplit(x,split='_')[[1]][c(3,4)], collapse='_') )
id2=lapply(sapply(files2,basename), function(x) paste0('matureSperm_',strsplit(x,split='[.]')[[1]][1], collapse='_') )
id=c(id1,id2)
names(id)=NULL

myobjDB=methRead(as.list(files),sample.id=id,assembly="hg38",
                 treatment=rep(0, length(id)),context="CpG",
                 dbtype = "tabix",dbdir = "methylDB")

print(myobjDB[[1]]@dbpath)
class(myobjDB[[1]])

Org_sampleCovSummary=lapply(seq(1,length(myobjDB)), function(x) summary(as(myobjDB[[x]],"GRanges")$coverage))

#srun --account mutationalscanning --mem 60g --pty bash

for(i in 1:length(id))
{
  
  png(paste0("methyKitPlots/", id[[i]], '_getMethylationStats.png'), width = 600, height = 750)
  par(mfrow = c(2, 1))
#png(paste0("methyKitPlots/", id[[i]], '_getMethylationStats.png'), width = 500, height = 450)
  getMethylationStats(myobjDB[[i]],plot=TRUE,both.strands=FALSE) # second samples
#dev.off()

#getMethylationStats(myobjDB[[i]],plot=FALSE,both.strands=FALSE)
#png(paste0("methyKitPlots/", id[[i]], '_getCoverageStats.png'), width = 500, height = 450)
  getCoverageStats(myobjDB[[i]],plot=TRUE,both.strands=FALSE)
dev.off()
  
}

save.image('methyKitPlots/wksp.RDATA')

# we will need to get the bases covered in all samples. The following function will merge all samples to one object for base-pair locations that are covered in all samples. 
# Setting destrand=TRUE (the default is FALSE) will merge reads on both strands of a CpG dinucleotide. This provides better coverage, but only advised when looking at CpG methylation (
  
meth=unite(myobjDB, destrand=FALSE)
head(meth)
unite(myobjDB[[1]], destrand=FALSE)
# creates a methylBase object, 
# where only CpGs covered with at least 1 sample per group will be returned


#The code below filters a methylRawList and discards bases that have coverage below 10X 
#and also discards the bases that have more than 99.9th percentile of coverage in each sample.
filtered.myobj=filterByCoverage(myobjDB,lo.count=8,lo.perc=NULL,
                                hi.count=NULL,hi.perc=99.9)


filt_sampleCovSummary=lapply(seq(1,length(filtered.myobj)), function(x) summary(as(filtered.myobj[[x]],"GRanges")$coverage))

filtered.meth=unite(filtered.myobj, destrand=FALSE)


# there were  groups defined by the treatment vector, 
# given during the creation of myobj: treatment=c(1,1,0,0)
meth.min=unite(myobj,min.per.group=1L)

jpeg(paste0("methyKitPlots/", 'MethySampleCorr.jpeg'), width = 450, height = 350)
  getCorrelation(meth,plot=TRUE)
dev.off()

png(paste0("methyKitPlots/", 'MethySampleCluster.png'), width = 600, height = 850)
  clusterSamples(meth, dist="correlation", method="ward", plot=TRUE)
dev.off()

jpeg(paste0("methyKitPlots/", 'MethySamplePCA_copmp.jpeg'), width = 450, height = 350)
  PCASamples(meth, screeplot=TRUE) # principle component
dev.off()

png(paste0("methyKitPlots/", 'MethySamplePCA.png'), width = 750, height = 750)
  PCASamples(meth)
dev.off()

methy_gr=as(meth,"GRanges")
#methy_gr=as(meth,"GRanges")

cov_df=mcols(methy_gr)[seq(1,ncol(mcols(methy_gr)), 3)]
colnames(cov_df)= unlist(id)

## test ##
tiles = tileMethylCounts(myobjDB,win.size=1000,step.size=1000,cov.bases = 10)
tiles[[11]]
tile11_gr=as(tiles[[11]],"GRanges")
methyGr_tile11=as(myobjDB[[11]],"GRanges")
temp_gr=subsetByOverlaps(methyGr_tile11,tile11_gr[1])
##

tiles2=tileMethylCounts(meth,win.size=1000,step.size=1000,cov.bases = 10)
tiles2_gr=as(tiles2,"GRanges")
tiles2_gr$methyPerc=percMethylation(tiles2)

tiles2_gr[11]
tiles2_gr[11]$methyPerc[,'PrimarySpermatocytes_179958']
a=subsetByOverlaps(methy_gr,tile11_gr[1])

tiles2_gr[which(seqnames(tiles2_gr)=='chr1' & start(tiles2_gr)==14001 & end(tiles2_gr)== 15000 ) ]

save.image('methyKitPlots/wksp.RDATA')
#load('methyKitPlots/wksp.RDATA')

# batch correction ##
sampleAnnotation=data.frame(batch_id=rep(c(1,2,3), each=4))
as=assocComp(mBase=meth,sampleAnnotation)
newObj=removeComp(meth,comp=1)
mat=percMethylation(meth)
mat[mat==100]=80
# reconstruct the methylBase from the corrected matrix
newobj=reconstruct(mat,meth)






##### CHN analysis : "CAA" "CTC" "CCT" "CAC" "CTT" "CCA" "CAT" "CTA" "CCC" ###. 




files1=list.files(path=paste0(data.dir,'07_methyl_dackel/Muenster/'), pattern='*_final_CHH.methylKit.gz$', full.names = T)
files2=list.files(path=paste0(data.dir,'07_methyl_dackel/PRJEB34432_NormalSperm/'), pattern='*_CHH.methylKit.gz$', full.names = T)
files=c(files1,files2)
id1=lapply(sapply(files1,basename), function(x) paste(strsplit(x,split='_')[[1]][c(3,4)], collapse='_') )
id2=lapply(sapply(files2,basename), function(x) paste0('matureSperm_',strsplit(x,split='[.]')[[1]][1], collapse='_') )
id=c(id1,id2)
names(id)=NULL

myobjDB_CHH=methRead(as.list(files),sample.id=id,assembly="hg38",
                 treatment=rep(0, length(id)),context="CHH",
                 dbtype = "tabix",dbdir = "methylDB_CHH")

print(myobjDB_CHH[[1]]@dbpath)
class(myobjDB_CHH[[1]])

OrgCHH_sampleCovSummary=lapply(seq(1,length(myobjDB_CHH)), function(x) summary(as(myobjDB_CHH[[x]],"GRanges")$coverage))

#srun --account mutationalscanning --mem 60g --pty bash

for(i in 1:length(id))
{
  
  png(paste0("methyKitPlots/", id[[i]], '_getMethylationStatsCHH.png'), width = 600, height = 750)
  par(mfrow = c(2, 1))
  #png(paste0("methyKitPlots/", id[[i]], '_getMethylationStats.png'), width = 500, height = 450)
  getMethylationStats(myobjDB_CHH[[i]],plot=TRUE,both.strands=FALSE) # second samples
  #dev.off()
  
  #getMethylationStats(myobjDB[[i]],plot=FALSE,both.strands=FALSE)
  #png(paste0("methyKitPlots/", id[[i]], '_getCoverageStats.png'), width = 500, height = 450)
  getCoverageStats(myobjDB_CHH[[i]],plot=TRUE,both.strands=FALSE)
  dev.off()
  
}







