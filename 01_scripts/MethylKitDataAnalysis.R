conda activate R.4.2_cloned_mapdemo_env
R
setwd('/home/vinodsingh/mutationalscanning/Workspaces/vinod/codes/bwa-meth_pipeline-master')
data.dir='/home/vinodsingh/mutationalscanning/Workspaces/vinod/codes/bwa-meth_pipeline-master/'

library(methylKit)
files1=list.files(path=paste0(data.dir,'07_methyl_dackel/Muenster/'), pattern='*_final_CpG.methylKit.gz$', full.names = T)
files2=list.files(path=paste0(data.dir,'07_methyl_dackel/PRJEB34432_NormalSperm/'), pattern='*_CpG.methylKit.gz$', full.names = T)
files3=list.files(path=paste0(data.dir,'07_methyl_dackel/PRJNA357085_Spermatogonia/'), pattern='*_CpG.methylKit.gz$', full.names = T)

files=c(files1,files2,files3)
id1=lapply(sapply(files1,basename), function(x) paste(strsplit(x,split='_')[[1]][c(3,4)], collapse='_') )
id2=lapply(sapply(files2,basename), function(x) paste0('matureSperm_',strsplit(x,split='[.]')[[1]][1], collapse='_') )
id3=lapply(sapply(files3,basename), function(x) paste0('Spgonia_',paste(strsplit(x,split='_')[[1]][c(2,3)], collapse='_') ))

id=c(id1,id2, id3)
names(id)=NULL

CpGmethy_objDB=methRead(as.list(files),sample.id=id,assembly="hg38",
                 treatment=rep(0, length(id)),context="CpG",
                 dbtype = "tabix",dbdir = "Sperm_CpGmethylDB")

print(CpGmethy_objDB[[1]]@dbpath)
class(CpGmethy_objDB[[1]])

Org_sampleCovSummary=lapply(seq(1,length(CpGmethy_objDB)), function(x) summary(as(CpGmethy_objDB[[x]],"GRanges")$coverage))

#srun --account mutationalscanning --mem 60g --pty bash
dir.create('Sperm_CpGmethyKitPlots')

for(i in 1:length(id))
{
  
  png(paste0("Sperm_CpGmethyKitPlots/", id[[i]], '_getMethylationStats.png'), width = 600, height = 750)
  par(mfrow = c(2, 1))
#png(paste0("methyKitPlots/", id[[i]], '_getMethylationStats.png'), width = 500, height = 450)
  getMethylationStats(CpGmethy_objDB[[i]],plot=TRUE,both.strands=FALSE) # second samples
#dev.off()

#getMethylationStats(myobjDB[[i]],plot=FALSE,both.strands=FALSE)
#png(paste0("methyKitPlots/", id[[i]], '_getCoverageStats.png'), width = 500, height = 450)
  getCoverageStats(CpGmethy_objDB[[i]],plot=TRUE,both.strands=FALSE)
dev.off()
  
}

#save.image('Sperm_CpGmethyKitPlots/wksp.RDATA')

# we will need to get the bases covered in all samples. The following function will merge all samples to one object for base-pair locations that are covered in all samples. 
# Setting destrand=TRUE (the default is FALSE) will merge reads on both strands of a CpG dinucleotide. This provides better coverage, but only advised when looking at CpG methylation (
  
CpGmethMerge=unite(CpGmethy_objDB, destrand=FALSE)
head(CpGmethMerge)
CpGmethMerge_gr=as(CpGmethMerge,"GRanges")
cov_df=mcols(CpGmethMerge_gr)[seq(1,ncol(mcols(CpGmethMerge_gr)), 3)]
colnames(cov_df)= unlist(id)

C_cnts_df=as.data.frame(mcols(CpGmethMerge_gr)[seq(2,ncol(mcols(CpGmethMerge_gr)), 3)])
T_cnts_df=as.data.frame(mcols(CpGmethMerge_gr)[seq(3,ncol(mcols(CpGmethMerge_gr)), 3)])
MethyDf=C_cnts_df/Reduce('+',list(C_cnts_df,T_cnts_df))
colnames(MethyDf)=id

summary_df=sapply(cov_df, summary)
summary_df[3,]

library(openxlsx)

MethyDf_M=MethyDf[,c(1,5,9,2,6,10,3,7,11,4,8,12,seq(13,17),18)]

corrDF=as.data.frame(cor(MethyDf_M))
#corrDF=corrDF[ ,order(names(corrDF))]
#corrDF=corrDF[ order(row.names(corrDF)), ]
corrDF=signif(corrDF, digits = 3)

write.xlsx(corrDF,file=paste0('Sperm_CpGmethyKitPlots','/','methycorrData.xlsx'),colNames = TRUE, rowNames = TRUE  )


#save.image('Sperm_CpGmethyKitPlots/wksp.RDATA')

#unite(CpGmethy_objDB[[1]], destrand=FALSE)
# creates a methylBase object, 
# where only CpGs covered with at least 1 sample per group will be returned


#The code below filters a methylRawList and discards bases that have coverage below 10X 
#and also discards the bases that have more than 99.9th percentile of coverage in each sample.
#filtered.CpGmethy_objDB=filterByCoverage(CpGmethy_objDB,lo.count=8,lo.perc=NULL,hi.count=NULL,hi.perc=99.9)
filtered.CpGmethy_objDB=filterByCoverage(CpGmethy_objDB,lo.count=8,lo.perc=NULL,hi.count=35,hi.perc=NULL)

filtered.CpGmethy_objDB_norm = normalizeCoverage(filtered.CpGmethy_objDB,method="median")

#filtered.CpGmethy_objDB_sperm=new("methylRawList", filtered.CpGmethy_objDB_norm[seq(1,12)], treatment = rep(1,12))
filtered.CpGmethy_objDB_sperm=new("methylRawList", filtered.CpGmethy_objDB_norm[seq(1,12)], treatment = rep(1,12))
CpGmethMerge_sperm=unite(filtered.CpGmethy_objDB_sperm, destrand=FALSE)
new_sample.ids=CpGmethMerge_sperm@sample.ids[c(1,5,9,2,6,10,3,7,11,4,8,12)]
CpGmethMerge_sperm_reorder=reorganize(CpGmethMerge_sperm,sample.ids=new_sample.ids,treatment = rep(1,12),save.db = FALSE,)


#filt_sampleCovSummary=lapply(seq(1,length(filtered.CpGmethy_objDB)), function(x) summary(as(filtered.CpGmethy_objDB[[x]],"GRanges")$coverage))


#filt.CpGmethMerge=unite(filtered.CpGmethy_objDB, destrand=FALSE)
#head(filt.CpGmethMerge)

CpGmethMerge=CpGmethMerge_sperm_reorder # replace
id2=CpGmethMerge@sample.ids

filt.CpGmethMerge_gr=as(CpGmethMerge,"GRanges")
filt.cov_df=as.data.frame(mcols(filt.CpGmethMerge_gr)[seq(1,ncol(mcols(filt.CpGmethMerge_gr)), 3)])
colnames(filt.cov_df)= unlist(id2)

filt.C_cnts_df=as.data.frame(mcols(filt.CpGmethMerge_gr)[seq(2,ncol(mcols(filt.CpGmethMerge_gr)), 3)])
filt.T_cnts_df=as.data.frame(mcols(filt.CpGmethMerge_gr)[seq(3,ncol(mcols(filt.CpGmethMerge_gr)), 3)])
filt.MethyDf=filt.C_cnts_df/Reduce('+',list(filt.C_cnts_df,filt.T_cnts_df))
colnames(filt.MethyDf)=id2

filt.summary_df=sapply(filt.cov_df, summary)
filt.summary_df[3,]

#filt.MethyDf_M=filt.MethyDf[,c(1,5,9,2,6,10,3,7,11,4,8,12,seq(13,17),18)]
#filt.MethyDf_M=filt.MethyDf[,c(1,5,9,2,6,10,3,7,11,4,8,12)]
filt.MethyDf_M=filt.MethyDf

filt.corrDF=as.data.frame(cor(filt.MethyDf_M))

filt.corrDF=signif(filt.corrDF, digits = 3)

write.xlsx(filt.corrDF,file=paste0('Sperm_CpGmethyKitPlots','/','filtNormSperm_methycorrData.xlsx'),colNames = TRUE, rowNames = TRUE  )

library(PerformanceAnalytics)
subsample1=sample(x=seq(1, nrow(filt.MethyDf_M)), size=100000, replace = FALSE, prob = NULL)
png(paste0("Sperm_CpGmethyKitPlots/", 'filtNormSperm_methyCorrData.png'), width = 850, height = 850)
chart.Correlation(filt.MethyDf_M[subsample1,], histogram=TRUE, pch=19);
dev.off()

## for coverage ##

library("tibble")
filt.cov_df_t=as_tibble(filt.cov_df)
filt.cov_corr=signif(as.data.frame(cor(filt.cov_df)), digits = 3)
upper_tri=upper.tri(filt.cov_corr, diag = FALSE)
filt.cov_corr[upper.tri(filt.cov_corr,diag = FALSE)] <- NA
write.xlsx(filt.cov_corr,file=paste0('Sperm_CpGmethyKitPlots','/','filtNormSperm_CoverageCorrData.xlsx'),colNames = TRUE, rowNames = TRUE  )

library(PerformanceAnalytics)
subsample=sample(x=seq(1, nrow(filt.cov_df_t)), size=100000, replace = FALSE, prob = NULL)
png(paste0("Sperm_CpGmethyKitPlots/", 'filtNormSperm_CoverageCorrData.png'), width = 850, height = 850)
  chart.Correlation(filt.cov_df_t[subsample,], histogram=TRUE, pch=19);
dev.off()

save.image('Sperm_CpGmethyKitPlots/wksp.RDATA')



# there were  groups defined by the treatment vector, 
# given during the creation of myobj: treatment=c(1,1,0,0)

png(paste0("Sperm_CpGmethyKitPlots/", 'filtNormSperm_MethySampleCorr.png'), width = 850, height = 850)
  getCorrelation(CpGmethMerge,plot=TRUE) # default is 2 million records
dev.off()


png(paste0("Sperm_CpGmethyKitPlots/", 'filtSperm_MethySampleCluster.png'), width = 600, height = 850)
  clusterSamples(CpGmethMerge, dist="correlation", method="ward", plot=TRUE)
dev.off()

png(paste0("Sperm_CpGmethyKitPlots/", 'filtSperm_MethySamplePCA_copmp.png'), width = 450, height = 350)
  PCASamples(CpGmethMerge, screeplot=TRUE) # principle component
dev.off()

png(paste0("Sperm_CpGmethyKitPlots/", 'filtSperm_MethySamplePCA.png'), width = 750, height = 750)
  PCASamples(CpGmethMerge)
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


data.dir='/home/vinodsingh/mutationalscanning/Workspaces/vinod/codes/bwa-meth_pipeline-master/'
library(methylKit)

CHH_files1=list.files(path=paste0(data.dir,'07_methyl_dackel/Muenster/'), pattern='*_final_CHH.methylKit.gz$', full.names = T)
CHH_files2=list.files(path=paste0(data.dir,'07_methyl_dackel/PRJEB34432_NormalSperm/'), pattern='*_CHH.methylKit.gz$', full.names = T)
CHH_files3=list.files(path=paste0(data.dir,'07_methyl_dackel/PRJNA357085_Spermatogonia/'), pattern='*_CHH.methylKit.gz$', full.names = T)
CHH_files=c(CHH_files1,CHH_files2,CHH_files3)

CHH_id1=lapply(sapply(CHH_files1,basename), function(x) paste(strsplit(x,split='_')[[1]][c(3,4)], collapse='_') )
CHH_id2=lapply(sapply(CHH_files2,basename), function(x) paste0('matureSperm_',strsplit(x,split='[.]')[[1]][1], collapse='_') )
CHH_id3=lapply(sapply(CHH_files3,basename), function(x) paste0('Spgonia_',paste(strsplit(x,split='_')[[1]][c(2,3)], collapse='_') ))
CHH_id=c(CHH_id1,CHH_id2,CHH_id3)

names(CHH_id)=NULL

CHHmethy_objDB=methRead(as.list(CHH_files),sample.id=CHH_id,assembly="hg38",
                 treatment=rep(0, length(id)),context="CHH",
                 dbtype = "tabix",dbdir = "Sperm_CHHmethylDB")

print(CHHmethy_objDB[[1]]@dbpath)
class(CHHmethy_objDB[[1]])

OrgCHH_sampleCovSummary=lapply(seq(1,length(CHHmethy_objDB)), function(x) summary(as(CHHmethy_objDB[[x]],"GRanges")$coverage))

#srun --account mutationalscanning --mem 60g --pty bash

for(i in 1:length(id))
{
  
  png(paste0("Sperm_CpGmethyKitPlots/", id[[i]], '_getMethylationStatsCHH.png'), width = 600, height = 750)
  par(mfrow = c(2, 1))
  #png(paste0("methyKitPlots/", id[[i]], '_getMethylationStats.png'), width = 500, height = 450)
  getMethylationStats(myobjDB_CHH[[i]],plot=TRUE,both.strands=FALSE) # second samples
  #dev.off()
  
  #getMethylationStats(myobjDB[[i]],plot=FALSE,both.strands=FALSE)
  #png(paste0("methyKitPlots/", id[[i]], '_getCoverageStats.png'), width = 500, height = 450)
  getCoverageStats(myobjDB_CHH[[i]],plot=TRUE,both.strands=FALSE)
  dev.off()
  
}







