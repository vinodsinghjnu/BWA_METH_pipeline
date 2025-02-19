library("optparse")

option_list = list(
  make_option(c("-D", "--BAM_Dir"), type="character", default=NULL, 
              help="BAM Files directory", metavar="character"),
  make_option(c("-P", "--BAMfiles_pattern"), type="character", default=NULL, 
              help="Pattern of Bam Files", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

BAM_filesDir=opt$BAM_Dir
BAM_filesPat=opt$BAMfiles_pattern

print('====>Parameters:')

print(paste0('BAM Files directory: ', BAM_filesDir))
print(paste0('Pattern of Bam Files: ', BAM_filesPat))

args <- commandArgs()
print(paste(c('command:', args), collapse=' '))

cat("<<======\n\n")

shhh <- suppressPackageStartupMessages # It's a library, so shhh!
library(Rsamtools); library(GenomicAlignments); shhh(library(CustomBioInfoFuctionsHumanGenome))
#setwd('/home/vinodsingh/mutationalscanning/Workspaces/vinod/codes/bwa-meth_pipeline-master')
# bam.list=list.files(path = '06_deduplicated/Muenster/', pattern = '*179944_final.bam$', full.names = T)
bam.list= list.files(path = BAM_filesDir, pattern = paste0('*',BAM_filesPat,"$"), full.names = T)

print(bam.list)

for(i in 1:length(bam.list))
{
  bfl=BamFile(bam.list[i])
  print(bfl)
  what <- c("rname", "strand", "pos")
  param <- ScanBamParam(what=what, tag = c("NM"))
  
  bam <- scanBam(bfl, param=param)
  print(largeVariables(n = 3))
  head(bam[[1]]$tag$NM)
  summary(bam[[1]]$tag$NM)
  #hist(bam[[1]]$tag$NM)
  
  library(ggplot2)
  df=data.frame(NM=bam[[1]]$tag$NM)
  p<-ggplot(df, aes(x=NM)) + 
    geom_histogram(color="black", fill="white")
  #p
  
  NM_distPlotName=sub('filt.bam$', 'NMtag_dist.png', bam.list[i]) 
  print(NM_distPlotName)
  ggsave(NM_distPlotName,p)
  
  print('==> successful')
  rm(bam)

}





