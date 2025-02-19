
library("optparse")

option_list = list(
  make_option(c("-I", "--FASTQDIR"), type="character", default=NULL, 
              help="path of fastq files", metavar="character"),
  make_option(c("-O", "--REPORTDIR"), type="character", default=NULL, 
              help="output dir for Fastqc report", metavar="character"),
  make_option(c("-T", "--THREADS"), type="interger", default=NULL, 
              help="output dir for Fastqc report", metavar="integer")
); 



library(openxlsx);library(fastqcr)


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if(is.null(opt)){
  print_help(opt_parser)
  stop("Input and Output dirs should be supplied", call.=FALSE)
}


## input ##
FASTQDIR=opt$FASTQDIR
REPORTDIR=opt$REPORTDIR
THREADS=opt$THREADS
###



print('====>Parameters:')

print(paste0('Input fastq DIR: ', FASTQDIR))
print(paste0('Output QC results Dir: ', REPORTDIR))
print(paste0('Threads: ', THREADS))


args <- commandArgs()
print(paste(c('command:', args), collapse=' '))

cat("<<======\n\n")

##

fastqc(fq.dir = FASTQDIR, qc.dir = REPORTDIR, threads = THREADS)
print(list.files(REPORTDIR))


#http://www.sthda.com/english/wiki/fastqcr-an-r-package-facilitating-quality-controls-of-sequencing-data-for-large-numbers-of-samples


wb <- createWorkbook()

## aggregate ##
REPORTDIR='/Users/au734493/mount/genomedk_mutScanning/codes/bwa-meth_pipeline-master/03_raw_data/MuensterMethyFastq/fastqc_analysis/'
QC <- qc_aggregate(REPORTDIR, progressbar = TRUE)
print(QC)

addWorksheet(wb, "FastQCSummary")
writeData(wb, "FastQCSummary", QC)



library(dplyr); library(tidyverse)
QC2= QC %>% dplyr::select(sample, module, status) %>% filter(status %in% c("WARN", "FAIL")) %>% arrange(sample)

addWorksheet(wb, "FastQCSummary2")
writeData(wb, "FastQCSummary2", QC2)

QC3=summary(QC)

addWorksheet(wb, "FastQCSummary3")
writeData(wb, "FastQCSummary3", QC3)

QC_satats=qc_stats(QC)
addWorksheet(wb, "FastQCStats")
writeData(wb, "FastQCStats", QC_satats)

OUTDIR='/Users/au734493/mount/genomedk_mutScanning/codes/bwa-meth_pipeline-master/03_raw_data/MuensterMethyFastq/'
saveWorkbook(wb, file = paste0(OUTDIR,'/QC_reportSummaryTable.xlsx'))


qc_report(REPORTDIR, result.file = paste0(OUTDIR, 'MuensterMethyFastq_QCreport.html'), experiment = "MuensterMethyData")


print("==completed successfully==")

