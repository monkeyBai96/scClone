#!/bin/bash  
#SBATCH --partition=cpu
#SBATCH -n 10             
#SBATCH --ntasks-per-node=10

#cellranger_path="/lustre/home/acct-biown/biown-mkbai/project/24clone/PMID37386030/${sample}/outs"
bam_path="/lustre/home/acct-biown/biown-mkbai/project/24clone/PMID32579974/bam/${sample}_scRNA.bam"
demo_path="/lustre/home/acct-biown/biown-mkbai/project/24clone/PMID32579974/${sample}"
ref=/lustre/home/acct-biown/biown-mkbai/reference/refdata-gex-hg19-3.0.0/fasta/genome.fa
annovar=/lustre/home/acct-biown/biown-mkbai/software/annovar/table_annovar.pl
strelka=/lustre/home/acct-biown/biown-mkbai/software/strelka-2.9.2/bin/configureStrelkaGermlineWorkflow.py


cd $demo_path
# array to 10 cells
eval IFS='+' read -ra cell_array <<< "$array"
length=${#cell_array[@]}

# call variants
for ((i=0; i<$length; i++ )); do
	cell=${cell_array[i]}
	module load samtools

	# change chr name
	#module load bcftools
	#samtools view -H $demo_path/${sample}_split/${cell}.bam > $demo_path/${sample}_split/${cell}.header.sam
	#sed '/^@SQ/{s/SN:\([0-9XY]\)/SN:chr\1/g;s/SN:MT/SN:chrM/g}' $demo_path/${sample}_split/${cell}.header.sam > $demo_path/${sample}_split/${cell}.header.chr.sam
	#samtools reheader $demo_path/${sample}_split/${cell}.header.chr.sam $demo_path/${sample}_split/${cell}.bam > $demo_path/${sample}_split/${cell}.chr.bam
	#module unload bcftools
	#
	samtools index $demo_path/${sample}_split/${cell}.bam
	module unload samtools
	mkdir $demo_path/strelka_tem_${cell}
	$strelka --bam ./${sample}_split/${cell}.bam \
	 	--referenceFasta $ref \
	 	--runDir $demo_path/strelka_tem_${cell}
	$demo_path/strelka_tem_${cell}/runWorkflow.py -m local -j 40

	# move variant file
	mkdir $demo_path/variants
	cp $demo_path/strelka_tem_${cell}/results/variants/variants.vcf.gz $demo_path/variants/${cell}.vcf.gz
	rm -rf $demo_path/strelka_tem_${cell}
done

exit 0
