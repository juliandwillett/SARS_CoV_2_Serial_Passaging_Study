#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --account=def-ioannisr 
#SBATCH --job-name='Make merged vcf file'
#SBATCH --mem-per-cpu=8000M
#SBATCH --cpus-per-task=4
#SBATCH --array=1-95%20
#1-95%1   

#https://www.ebi.ac.uk/sites/ebi.ac.uk/files/content.ebi.ac.uk/materials/2014/140217_AgriOmics/dan_bolser_snp_calling.pdf
#https://hbctraining.github.io/In-depth-NGS-Data-Analysis-Course/sessionVI/lessons/02_variant-calling.html

module load StdEnv/2020 gcc/9.3.0 bwa/0.7.17 samtools/1.15.1 freebayes/1.3.6 vcftools/0.1.16

cd /lustre04/scratch/willett/ViralCellStudies/

REF_FILE=References/MN908947.3.fa
BATCH1_CURR_FOLDER=$(sed -n "${SLURM_ARRAY_TASK_ID}p" AllFQ_Batch1_forMerge)
BATCH2_CURR_FOLDER=$(sed -n "${SLURM_ARRAY_TASK_ID}p" AllFQ_Batch2_forMerge)

FILE1_B1=FastqFiles_Mouse_FirstRun/${BATCH1_CURR_FOLDER}/*_R1_*.fastq.gz
FILE2_B1=FastqFiles_Mouse_FirstRun/${BATCH1_CURR_FOLDER}/*_R2_*.fastq.gz

cat FastqFiles_Mouse_SecondRun/$BATCH2_CURR_FOLDER/*R1_001.fastq.gz > Vcfs_Merged/${BATCH2_CURR_FOLDER}_R1_B2.fastq.gz #make merged file
cat FastqFiles_Mouse_SecondRun/$BATCH2_CURR_FOLDER/*R2_001.fastq.gz > Vcfs_Merged/${BATCH2_CURR_FOLDER}_R2_B2.fastq.gz #make merged file

FILE1_B2=Vcfs_Merged/${BATCH2_CURR_FOLDER}_R1_B2.fastq.gz
FILE2_B2=Vcfs_Merged/${BATCH2_CURR_FOLDER}_R2_B2.fastq.gz

cat $FILE1_B1 $FILE1_B2 > Vcfs_Merged/${BATCH2_CURR_FOLDER}_file1.fastq.gz
cat $FILE2_B1 $FILE2_B2 > Vcfs_Merged/${BATCH2_CURR_FOLDER}_file2.fastq.gz

PIPE_FILE1=Vcfs_Merged/${BATCH2_CURR_FOLDER}_file1.fastq.gz
PIPE_FILE2=Vcfs_Merged/${BATCH2_CURR_FOLDER}_file2.fastq.gz

MAX_READS=34253784 #min num reads in SAM file after BWA MEM step for all samples

#bwa index $REF_FILE #Indexing must be done before running script
bwa mem -M -t ${SLURM_CPUS_PER_TASK} $REF_FILE $PIPE_FILE1 $PIPE_FILE2 | \
samtools fixmate -u -m - - | \
samtools sort -u -@4 -T Vcfs_Merged/$BATCH2_CURR_FOLDER - | \
samtools markdup -@4 --reference $REF_FILE - Vcfs_Merged/${BATCH2_CURR_FOLDER}.cram

SAMPLE_READS=$(samtools view -c Vcfs_Merged/${BATCH2_CURR_FOLDER}.cram)
FRACTION=$(bc <<< "scale=10 ; $MAX_READS / $SAMPLE_READS") 

echo $BATCH2_CURR_FOLDER $SAMPLE_READS $MAX_READS $FRACTION >> merging_vcfs_numbers_fractions.txt

samtools view -s $FRACTION -@4 Vcfs_Merged/${BATCH2_CURR_FOLDER}.cram -o Vcfs_Merged/$BATCH2_CURR_FOLDER.bam
samtools faidx $REF_FILE
freebayes -f $REF_FILE Vcfs_Merged/$BATCH2_CURR_FOLDER.bam | \
vcftools --vcf - --minQ 20 --minDP 10 --recode --recode-INFO-all --out Vcfs_Merged/$BATCH2_CURR_FOLDER
rm Vcfs_Merged/${BATCH2_CURR_FOLDER}.cram
rm Vcfs_Merged/$BATCH2_CURR_FOLDER.bam
rm $PIPE_FILE1
rm $PIPE_FILE2
rm $FILE1_B2
rm $FILE2_B2

echo Done
