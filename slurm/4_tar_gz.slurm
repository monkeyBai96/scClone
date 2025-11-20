#!/bin/bash

#SBATCH --job-name=tar
#SBATCH --partition=cpu
#SBATCH -n 20
#SBATCH --ntasks-per-node=40

# user
path=/lustre/home/acct-biown/biown-mkbai/project/24clone/PMID32579974

sample="P4_Normal"
cd ${path}/${sample}
tar -czvf annovar.tar.gz ./annovar/*.txt
tar -czvf variants.tar.gz ./variants

sample="P6_Normal"
cd ${path}/${sample}
tar -czvf annovar.tar.gz ./annovar/*.txt
tar -czvf variants.tar.gz ./variants

sample="P4_Tumor"
cd ${path}/${sample}
tar -czvf annovar.tar.gz ./annovar/*.txt
tar -czvf variants.tar.gz ./variants

sample="P6_Tumor"
cd ${path}/${sample}
tar -czvf annovar.tar.gz ./annovar/*.txt
tar -czvf variants.tar.gz ./variants
