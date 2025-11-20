#!/bin/bash  
#SBATCH --partition=cpu
#SBATCH -n 40             
#SBATCH --ntasks-per-node=40 

#sample="____"
#cellranger_path="/lustre/home/acct-biown/biown-mkbai/project/24clone/PMID37386030/${sample}/outs"
bam_path="/lustre/home/acct-biown/biown-mkbai/project/24clone/PMID32579974/bam/${sample}_scRNA.bam"
demo_path="/lustre/home/acct-biown/biown-mkbai/project/24clone/PMID32579974/${sample}"
ref=/lustre/home/acct-biown/biown-mkbai/reference/refdata-gex-hg19-3.0.0/fasta/genome.fa
annovar=/lustre/home/acct-biown/biown-mkbai/software/annovar/table_annovar.pl
strelka=/lustre/home/acct-biown/biown-mkbai/software/strelka-2.9.2/bin/configureStrelkaGermlineWorkflow.py

# filter bam
mkdir -p $demo_path
cd $demo_path

module load samtools
dos2unix ../${sample}_barcode.txt
samtools view -H $bam_path > ${sample}_SAM_header
samtools view $bam_path | LC_ALL=C grep -F -f ../${sample}_barcode.txt > ${sample}_SAM_body
cat ${sample}_SAM_header ${sample}_SAM_body > ${sample}.sam
samtools view -b ${sample}.sam > ${sample}.bam
samtools sort -t CB ${sample}.sam > ${sample}_sort.bam
module unload samtools

# split bam
source conda activate py3
mkdir -p $demo_path/${sample}_split
python ../run_split.py $demo_path/${sample}_sort.bam $demo_path/${sample}_split/

exit 0
