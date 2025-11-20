#!/bin/bash

#SBATCH --job-name=annovar
#SBATCH --partition=cpu
#SBATCH -n 40
#SBATCH --ntasks-per-node=40

path=/lustre/home/acct-biown/biown-mkbai/project/24clone/PMID32579974
ref_fa=/lustre/home/acct-biown/biown-mkbai/reference/hg19/hg19.fa
annovar=/lustre/home/acct-biown/biown-mkbai/software/annovar


mkdir -p $path/${sample}/annovar

# annovar
cd $path/${sample}/variants
for i in *.vcf.gz
	do
	a=$(basename $i .vcf.gz) &&
	perl ${annovar}/table_annovar.pl ${a}.vcf.gz ${annovar}/humandb/ -buildver hg19 -out ../annovar/${a} -remove \
	 -protocol refGene,ALL.sites.2015_08,EAS.sites.2015_08,avsnp147,clinvar_20160302,cosmic70,esp6500siv2_ea,exac03,gnomad_exome \
	 -operation g,f,f,f,f,f,f,f,f -nastring . -vcfinput --otherinfo
done

exit 0
