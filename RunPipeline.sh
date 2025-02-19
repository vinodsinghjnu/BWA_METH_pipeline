#!/bin/bash


wd='/home/vinodsingh/mutationalscanning/Workspaces/vinod/codes/bwa-meth_pipeline-master'
cd $wd

source $HOME/miniforge3/etc/profile.d/conda.sh
conda activate ngs_packages
#conda activate ngs_packages


CODESDIR='/faststorage/project/mutationalscanning/Workspaces/vinod/codes/'
fastqc='/home/vinodsingh/installed_softwares/FastQC/./fastqc'
REF="02_reference/hg38/hg38.fa"



#wd=$(pwd) 
#DATADIR='/home/vinodsingh/mutationalscanning/Workspaces/vinod/codes/bwa-meth_pipeline-master/03_raw_data/MuensterMethyFastq/'
DATADIR='/home/vinodsingh/mutationalscanning/Workspaces/vinod/codes/bwa-meth_pipeline-master/03_raw_data/'
DATANAME="PRJEB34432_NormalSperm"
FASTQDIR=${DATADIR}"$DATANAME"
INPUT=03_raw_data/${DATANAME}

# List Samples for alignment #
fastqPat=".fastq.gz"; grepPat="R1${fastqPat}|1${fastqPat}"
ls "$INPUT"/*$fastqPat | grep -E $grepPat | parallel basename {} | awk -v pat="_$grepPat" -v repl="" '{ sub(pat,repl,$0); print $0 }' > $(basename $INPUT)_samples_for_alignment.txt  # for muenster published data



# Get read length 
readLengthFile=${DATADIR}/${DATANAME}/read_length.txt
cd $FASTQDIR
for f in *.fastq.gz; do
  rl=$(zcat $f | awk '{if(NR%4==2) print length($1)}' | sort -n | uniq -c)
  echo "$f $rl"
  #echo $(zcat $f|wc -l)/4|bc
done >$readLengthFile
cd $wd

### Create QC report ###

n_cpus=20
rs=$wd/01_scripts/fastQC_report.R
master_R_bash_script=${CODESDIR}/runR_masterCluster.sh
sbatch --account mutationalscanning -J fastQC_report --time 06:00:00 -c $n_cpus --mem 80G  $master_R_bash_script $rs  "-I $FASTQDIR -O ${DATADIR} -T $n_cpus"

cd $DATADIR
mkdir ${FASTQDIR}/fastqc_analysis
$fastqc  -f fastq --outdir=${FASTQDIR}/fastqc_analysis -t 4 --memory 10000 ${FASTQDIR}/*.fastq.gz
cd ${FASTQDIR}/fastqc_analysis
multiqc . --force
cd $wd




### index the reference with BSseeker2 ###
sh 01_scripts/00_index_ref.sh $REF

### Triming ###

INPUT=03_raw_data/${DATANAME}
DATANAME=$(basename $INPUT)

# make important dirs
mkdir 04_trimmed/${DATANAME}
mkdir 05_aligned/${DATANAME}
mkdir 06_deduplicated/${DATANAME}
mkdir 07_methyl_dackel/${DATANAME}
mkdir 10_logfiles/${DATANAME}

# decide read length threshold in 01_fastp_trimming from QC report #
n_cpus=20
sbatch --account mutationalscanning -J 01_fastp_trimming.sh --time 07:00:00 -c $n_cpus --mem 80G  01_scripts/01_fastp_trimming.sh $INPUT $n_cpus



### Again FastqC for trimmed reads data ###
trim_DATADIR=04_trimmed/${DATANAME}
mkdir 04_trimmed/${DATANAME}/fastqc_analysis
$fastqc  -f fastq --outdir=${trim_DATADIR}/fastqc_analysis -t 4 --memory 10000 ${trim_DATADIR}/*fastq.gz
cd ${trim_DATADIR}/fastqc_analysis
multiqc .
cd $wd


###########################
#### MWA-METH alignment ### 
###########################


#samples=$(ls -1 "$INPUT"/*_R1_001.fastq.gz | perl -pe 's/_R1_001.fastq\.gz//g' | parallel basename {}) # for munster cluster data


samples=$(ls -1 "$INPUT"/*_R1_001.fastq.gz | perl -pe 's/_R1_001.fastq\.gz//g' | parallel basename {}) # for munster cluster data

# NOTE #
# Sample path: 03_raw_data/MuensterMethyFastq/Muenster/A006200243_179937_S22_L001_R1_001.fastq.gz, 03_raw_data/MuensterMethyFastq/Muenster/A006200243_179937_S22_L001_R2_001.fastq.gz
# We want to extract smple name from this path "A006200243_179937_S22_L001" 

parallel  "sbatch --parsable --account mutationalscanning -J 02_bwa-meth.{} --time 22:00:00 -c $n_cpus --mem 70G 01_scripts/02_bwa-meth.sh $INPUT $n_cpus {}" ::: $samples


# varify run is complete or not ##
grep -E error $(ls -1 02_bwa-meth*.out) | wc -l
grep -E "==> completed" $(ls -1  02_bwa-meth*.out) | wc -l

mv $(ls -1 02_bwa-meth*.out) 10_logfiles/${DATANAME}/



## remove duplicates ##

n_cpusD=6
parallel "sbatch --parsable -J 03_remove_duplicates_{} --time 10:00:00 --cpus-per-task $n_cpusD --mem=60GB 01_scripts/03_remove_duplicates.sh  $INPUT {}" ::: $samples


grep -E "==> completed" $(ls  03_remove_duplicates_*.out)
grep -E error $(ls -1 03_remove_duplicates_*.out) | wc -l
grep -E "ERROR:" $(ls  03_remove_duplicates_*.out)
mv $(ls -1 03_remove_duplicates_*.out) 10_logfiles/${DATANAME}/


# run  script on cluster

# Extra reads filtering ## 03a_filterReads
# filter out duplicates, low quality reads, not mapped in a proper pair 

n_cpus=2
parallel "sbatch  --account mutationalscanning -J 03a_filterReads_{} --time 01:00:00 -c $n_cpus --mem 40G 01_scripts/03a_filter_reads.sh $INPUT {}" ::: $samples
ls -1 03a_filterReads_*.out | wc -l
grep -E error $(ls -1 03a_filterReads_*.out) # slurm errors
grep -E "ERROR:TRUNCATED_FILE" $(ls  03a_filterReads_*.out)
grep -E "ERROR:" $(ls  03a_filterReads_*.out)
mv $(ls -1 03a_filterReads_*.out) 10_logfiles/${DATANAME}/






## RUN methyl-dackel_mbias plot ##

rm 04_methyl-dackel_mbias.sh_*.out
n_cpus=6
parallel  "sbatch --account mutationalscanning  --time 2:00:00 -c $n_cpus --mem 80G -J 04_methyl-dackel_mbias.sh_{} 01_scripts/04_methyl-dackel_mbias.sh $INPUT $n_cpus $REF {}" ::: $(ls -1 06_deduplicated/${DATANAME}/*.dedup.filt.bam) 
grep -E "error" $(ls  04_methyl-dackel_mbias.sh_*.out)

mv $(ls -1 04_methyl-dackel_mbias.sh_*.out) 10_logfiles/${DATANAME}/



# write data and find threshold ## not necessary
n_cpus=6
parallel  "sbatch --account mutationalscanning  --time 2:00:00 -c $n_cpus --mem 80G -J 04a_methyl-dackel_mbias.sh_{1} sh 01_scripts/04a_methyl-dackel_mbias.sh $INPUT $n_cpus {} $" ::: $real_samples 
grep -E "error" $(ls  04a_methyl-dackel_mbias.sh_*.out)


## RUN methyl-dackel_extract  ##


# per base depth of bam files ##
BAM_pat=".dedup.filt.bam"
BAMreadDepthFile=06_deduplicated/${DATANAME}/BAM_ReadDepth.txt
for f in $(ls 06_deduplicated/${DATANAME}/*${BAM_pat}); do
  echo -e "\n$f"
  tot_size=$(samtools view -H $f | grep -P '^@SQ' | cut -f 3 -d ':' | awk '{sum+=$1} END {print sum}')
  samtools depth $f |  awk -v var="$tot_size" '{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/var; print "Stdev = ",sqrt(sumsq/var - (sum/var)**2)}'
done >$BAMreadDepthFile

 

n_cpus=2 # only two CPU  are enough
#BAM_pat='_final.bam'
BAM_pat=".dedup.filt.bam"
real_samples=$(ls 06_deduplicated/${DATANAME}/*$BAM_pat | parallel echo {/.} | cut -d'.' -f 1)
real_samples=${real_samples[@]}

parallel  "sbatch --account mutationalscanning  --time 04:00:00 -c $n_cpus --mem 60G -J 05_methyl-dackel_extract_up.sh_{1} 01_scripts/05_methyl-dackel_extract_up.sh $INPUT $n_cpus $BAM_pat {}" ::: $real_samples 
head $(ls -1 07_methyl_dackel/${DATANAME}/OTOB*)
grep -E "==> completed " $(ls 05_methyl-dackel_extract.sh_*.out) | wc -l
grep -E "error" $(ls 05_methyl-dackel_extract.sh_*.out) | wc -l

mv $(ls -1 05_methyl-dackel_extract.sh_*.out) 10_logfiles/${DATANAME}/

echo $f | xargs -I {} sh -c  "echo $(basename {} ${BAM_pat})_c"
#samtools merge -o 06_deduplicated/PRJNA357085_Spermatogonia/merge_SRR5099529_SRR5099530_final.bam -@ 3 $(ls -1 06_deduplicated/PRJNA357085_Spermatogonia/*dedup.filt.bam)
f=06_deduplicated/PRJNA357085_Spermatogonia/merge_SRR5099529_SRR5099530_final.bam
tot_size=$(samtools view -H $f | grep -P '^@SQ' | cut -f 3 -d ':' | awk '{sum+=$1} END {print sum}')
samtools depth $f |  awk -v var="$tot_size" '{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/var; print "Stdev = ",sqrt(sumsq/var - (sum/var)**2)}'








rscript1='01_scripts/MethylKitDataAnalysis.R'
rscript2='01_scripts/MethyFilesRead.R'






### testing codes ##









cp -r 01_scripts 10_logfiles/${DATANAME}/
cp RunPipeline_new2.sh 10_logfiles/${DATANAME}/01_scripts/


conda  activate samtools_1.16.1
GENOME="02_reference/hg38/hg38.fa"
NCPUS=7
bamFile="06_deduplicated/Muenster/179937_Sperm_UndiffSpermatogonia_179937_final.bam"
sample='test_179937'
cd "/faststorage/project/mutationalscanning/Workspaces/vinod/codes/bwa-meth_pipeline-master/"
MethylDackel mbias -@ "$NCPUS" --noCpG --CHH --CHG  "$GENOME" "$bamFile" "${sample}_chn"
MethylDackel mbias -@ $NCPUS --noCpG --CHH --CHG $GENOME $bamFile ${sample}_chn>&OTOBchn$sample.txt

MethylDackel extract $(cat OTOBchn$sample.txt | cut -d ':' -f 2) --noCpG --CHH --CHG --minOppositeDepth 1 --maxVariantFrac 0.1 "$GENOME" $bamFile \; \
gzip {.}_CHN.methylKit \; \

###### run on old data ##

oldInput='06_deduplicated/PRJEB34432_NormalSpermBam'
mkdir 07_methyl_dackel/PRJEB34432_NormalSpermBam
old_samples=('P4_NC002' 'P4_NC004' 'P4_NC005' 'P4_NC007' 'P4_NC008')
old_samples=${old_samples[@]}


n_cpus2=6
parallel  "sbatch --account mutationalscanning  --time 2:00:00 -c $n_cpus2 --mem 80G -J 04_methyl-dackel_mbias.sh_old_{1} 01_scripts/04_methyl-dackel_mbias.sh $oldInput $n_cpus {}" ::: $old_samples 
grep -E "error" $(ls  04_methyl-dackel_mbias.sh_old*.out)


parallel "sbatch --account mutationalscanning  --time 5:00:00 -c $n_cpus --mem 60G -J 05_methyl-dackel_extract.sh_old_{1} 01_scripts/05_methyl-dackel_extract.sh $oldInput $n_cpus {}" ::: $old_samples

