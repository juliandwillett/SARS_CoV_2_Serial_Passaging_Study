#!/bin/bash
#SBATCH --time=1:00:00
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

cat FastqFiles_Mouse_SecondRun/$BATCH2_CURR_FOLDER/*R1_001.fastq.gz > tmp/${BATCH2_CURR_FOLDER}_R1_B2.fastq.gz #make merged file
cat FastqFiles_Mouse_SecondRun/$BATCH2_CURR_FOLDER/*R2_001.fastq.gz > tmp/${BATCH2_CURR_FOLDER}_R2_B2.fastq.gz #make merged file

FILE1_B2=tmp/${BATCH2_CURR_FOLDER}_R1_B2.fastq.gz
FILE2_B2=tmp/${BATCH2_CURR_FOLDER}_R2_B2.fastq.gz

cat $FILE1_B1 $FILE1_B2 > tmp/${BATCH2_CURR_FOLDER}_file1.fastq.gz
cat $FILE2_B1 $FILE2_B2 > tmp/${BATCH2_CURR_FOLDER}_file2.fastq.gz

PIPE_FILE1=tmp/${BATCH2_CURR_FOLDER}_file1.fastq.gz
PIPE_FILE2=tmp/${BATCH2_CURR_FOLDER}_file2.fastq.gz

#bwa index $REF_FILE #Indexing must be done before running script
bwa mem -M -t ${SLURM_CPUS_PER_TASK} $REF_FILE $PIPE_FILE1 $PIPE_FILE2 | \
samtools fixmate -u -m - - | \
samtools sort -u -@4 -T tmp/$BATCH2_CURR_FOLDER - | \
samtools markdup -@4 --reference $REF_FILE - tmp/${BATCH2_CURR_FOLDER}.cram
samtools view -c tmp/${BATCH2_CURR_FOLDER}.cram >> num_reads_merged.txt

rm $FILE1_B2
rm $FILE2_B2
rm $PIPE_FILE1
rm $PIPE_FILE2
rm tmp/$BATCH2_CURR_FOLDER.cram
echo Done