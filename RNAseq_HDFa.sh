#!/bin/bash -l
#SBATCH -A naiss2023-22-935
#SBATCH -p core -n 14
#SBATCH -t 5-00:00:00
#SBATCH -J thbs1
#SBATCH -e thbs1_SLURM_Job_id=%j.sderr
#SBATCH -o thbs1_SLURM_Job_id=%j.sdout
#SBATCH --mail-type=All
#SBATCH --mail-user=zhuang.liu@ki.se
module load bioinfo-tools FastQC/0.11.8 java/sun_jdk1.8.0_151 trimmomatic/0.36
module load bioinfo-tools BEDTools/2.27.1 ucsc-utilities/v345 HISAT2/2.1.0 samtools/1.10

for SampleName in HDFa_Ctr_1 HDFa_Ctr_2 HDFa_OE_1 HDFa_OE_2
do
cd /crex/proj/snic2021-23-156/THBS1_NC_revision/RNAseq_HDFa/00-rawData/soapnuke/clean/${SampleName}
fastqc ${SampleName}_1.fq.gz
fastqc ${SampleName}_2.fq.gz

trimmomatic PE -threads 14 \
${SampleName}_1.fq.gz ${SampleName}_2.fq.gz \
ILLUMINACLIP:/proj/snic2019-8-262/humanGenome/Trimadapters/TruSeq3-PE-2.fa:2:30:10 \
SLIDINGWINDOW:4:15 LEADING:3 TRAILING:3 MINLEN:36 \
-baseout ${SampleName}.fq.gz

fastqc ${SampleName}_1P.fq.gz
fastqc ${SampleName}_2P.fq.gz

# hisat2 mapping the read to reference genome
##-- pay attention to the parameter of rna-strandness according to different library types
##-- which can be inferred by the RseQC https://chipster.csc.fi/manual/rseqc_infer_rnaseq_experiment.html
hisat2 -p 8 \
--rg-id=${SampleName} --rg SM:${SampleName} \
--rg LB:${SampleName}_DNBseq --rg PL:DNBSEQ \
-x /proj/snic2019-8-262/humanGenome/hisat2_index/genomeGRCh38_tran \
--dta --rna-strandness RF \
--known-splicesite-infile /proj/snic2019-8-262/humanGenome/gencodev34annotation.ss \
--no-softclip \
-1 ${SampleName}_1P.fq.gz \
-2 ${SampleName}_2P.fq.gz \
--summary-file ${SampleName}_mapping.log \
-S ${SampleName}.sam

samtools view -bS ${SampleName}.sam \
> ${SampleName}.bam
samtools sort -@ 8 ${SampleName}.bam \
-o ${SampleName}_sorted.bam
samtools index ${SampleName}_sorted.bam
rm ${SampleName}.sam
rm ${SampleName}.bam
done


# Bowtie2 together with RseQC infer the library type
GenomeIndex=/crex/proj/snic2019-8-262/humanGenome/bowtie2_index/genomeGRCh38Bowtie2Ver
bowtie2 -p 10 -q \
--very-sensitive \
-X 1000 \
-x $GenomeIndex \
-1 ${SampleName}_1P.fq.gz \
-2 ${SampleName}_2P.fq.gz \
-S ${SampleName}.sam \
2> ${SampleName}_mapping.log
##transform the sam to bam and sort, index
samtools view -bS ${SampleName}.sam > ${SampleName}.bam
samtools sort -@ 4 ${SampleName}.bam -o ${SampleName}_sorted.bam
samtools index ${SampleName}_sorted.bam
rm ${SampleName}.sam
rm ${SampleName}.bam

module load bioinfo-tools rseqc/2.6.4
infer_experiment.py \
-i /proj/snic2021-23-156/pro_radiated_KC_lncRNAseq/s1rawData/Clean/${SampleName}/${SampleName}_sorted.bam \
-r /crex/proj/snic2021-23-156/human_genome/hg38_Gencode_V28_RSeQC.bed \
> ./infer_experiment/${SampleName}.infer_experiment.txt


