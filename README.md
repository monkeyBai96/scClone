# scClone: de novo mutation detection and clonal inference tool from single-cell transcriptome sequencing

scClone is a comprehensive tool designed for **de novo mutation calling** and **clonal structure inference** using single-cell transcriptome sequencing data. It supports end-to-end analysis from raw sequencing data processing to clonal visualization, with optimized workflows for both high-performance computing (HPC/Linux) and local RStudio environments.

## GitHub Directory Structure

```bash

scClone/
├── package/          # Core R package of scClone
├── benchmark/        # Test code for myeloma cell line validation
├── case_1-4/         # Case study code (HCC, pancreatic cancer, ovarian cancer, cSCC)
├── slurm/            # Bash command run on HPC
├── document/         # Files required in the process
└── README.md         # User guide (this file)
```

## Installation

First install the scClone R package via `devtools`:

```r

# Install devtools if not already installed
if (!require("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# Install scClone from GitHub
devtools::install_github("monkeyBai96/scClone", subdir = "package")
```

## Part I: Mutation Calling (Recommended: HPC/Linux Server)

### 0. Preprocessing with Cell Ranger/Spaceranger

First generate BAM files using 10x Genomics Cell Ranger (for scRNA-seq) or Spaceranger (for spatial transcriptomics):

```bash

# Example Cell Ranger command (adjust parameters for your data)
cellranger count \
  --id=sample_id \
  --transcriptome=path/to/ref_genome \
  --fastqs=path/to/fastq_files \
  --sample=sample_name \
  --localcores=16 \
  --localmem=64
```

### 1. Split BAM by Cell Barcode

Extract and split BAM files using valid cell barcodes:

```bash

# Load required module
module load samtools

# Define input/output paths (modify as needed)
input_bam="path/to/cellranger/outs/possorted_genome_bam.bam"
sample_name="your_sample"
barcode_file="path/to/valid_barcodes.txt"
output_dir="path/to/split_bam_output"
mkdir -p ${output_dir}

# Extract SAM header and valid barcode reads
samtools view -H ${input_bam} > ${output_dir}/${sample_name}_header.sam
samtools view ${input_bam} | LC_ALL=C grep -F -f ${barcode_file} > ${output_dir}/${sample_name}_body.sam

# Reconstruct and sort SAM/BAM
cat ${output_dir}/${sample_name}_header.sam ${output_dir}/${sample_name}_body.sam > ${output_dir}/${sample_name}.sam
samtools view -b ${output_dir}/${sample_name}.sam > ${output_dir}/${sample_name}.bam
samtools sort -t CB ${output_dir}/${sample_name}.sam > ${output_dir}/${sample_name}_sorted.bam

# Split sorted BAM into single-cell BAMs (using Python script)
conda activate py3  # Activate environment with required dependencies
python ../run_split.py \
  --input_bam=${output_dir}/${sample_name}_sorted.bam \
  --output_dir=${output_dir}/${sample_name}_split/
```

### 2. Call Mutations with Strelka

Call somatic mutations for each single-cell BAM using Strelka:

```bash

# Load required module
module load samtools

# Define core parameters (modify as needed)
sample_dir="path/to/sample_output"
ref_genome="path/to/hg19.fa"  # Reference genome (hg19 recommended)
strelka_path="path/to/strelka-2.9.10/bin/configureStrelkaSomaticWorkflow.py"

# Iterate over split single-cell BAMs
for cell_bam in ${sample_dir}/${sample_name}_split/*.bam; do
  # Get cell ID from BAM filename
  cell_id=$(basename ${cell_bam} .bam)
  
  # Index BAM file
  samtools index ${cell_bam}
  
  # Create temporary Strelka directory
  strelka_temp_dir=${sample_dir}/strelka_temp_${cell_id}
  mkdir -p ${strelka_temp_dir}
  
  # Configure and run Strelka
  ${strelka_path} \
    --bam ${cell_bam} \
    --referenceFasta ${ref_genome} \
    --runDir ${strelka_temp_dir}
  
  # Execute workflow (40 threads, adjust based on server resources)
  ${strelka_temp_dir}/runWorkflow.py -m local -j 40
  
  # Move and organize variant files
  variant_output_dir=${sample_dir}/variants
  mkdir -p ${variant_output_dir}
  cp ${strelka_temp_dir}/results/variants/variants.vcf.gz ${variant_output_dir}/${cell_id}.vcf.gz
  
  # Clean up temporary files
  rm -rf ${strelka_temp_dir}
done
```

### 3. Variant Annotation with ANNOVAR

Annotate VCF files using ANNOVAR (functional, population, and clinical annotations):

```bash

# Define ANNOVAR paths (modify as needed)
annovar_path="path/to/annovar"
humandb_path="${annovar_path}/humandb"
variant_dir="path/to/variants"
annotation_output_dir="path/to/annovar_annotations"
mkdir -p ${annotation_output_dir}

# Iterate over VCF files for annotation
for vcf_file in ${variant_dir}/*.vcf.gz; do
  # Get cell ID from VCF filename
  cell_id=$(basename ${vcf_file} .vcf.gz)
  
  # Run ANNOVAR (hg19 build, multiple annotation databases)
  ${annovar_path}/table_annovar.pl \
    ${vcf_file} \
    ${humandb_path}/ \
    -buildver hg19 \
    -out ${annotation_output_dir}/${cell_id} \
    -remove \
    -protocol refGene,ALL.sites.2015_08,EAS.sites.2015_08,avsnp147,clinvar_20160302,cosmic70,esp6500siv2_ea,exac03,gnomad_exome \
    -operation g,f,f,f,f,f,f,f,f \
    -nastring . \
    -vcfinput \
    --otherinfo
done

# Compress results for storage
tar -czvf ${annotation_output_dir}/annovar_annotations.tar.gz ${annotation_output_dir}/*.txt
tar -czvf ${variant_dir}/raw_variants.tar.gz ${variant_dir}/*.vcf.gz
```

## Part II: Genotype Inference & Clonal Structure Analysis (Recommended: Local RStudio, RAM > 16G)

### 0. Single-Cell Expression Preprocessing with Seurat

Process scRNA-seq expression data to get cell clusters and metadata:
```r

obj <- CreateSeuratObject(counts = raw, min.cells = 5, meta.data = data.frame(row.names = meta$cell, sample = meta$sample, tissue = meta$tissue))
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)
obj <- FindNeighbors(obj, dims = 1:10)
obj <- FindClusters(obj, resolution = 1)
obj <- RunUMAP(obj, dims = 1:10)
obj <- RunTSNE(obj, dims = 1:10)
DimPlot(object = obj, reduction = 'umap',label = TRUE, group.by = "seurat_clusters", pt.size = 3,raster = TRUE)
saveRDS(obj, "demo.rds")
```

### 1. Extract raw mutation & statistic

Extract raw mutation information and VAF from cell variants and statistics:
```r
col_name <- c("Chr","Pos","Ref","Alt","Func.refGene","Gene.refGene","GeneDetail.refGene","ExonicFunc.refGene","AAChange.refGene","Qual","Filter","Info","Format","VAF","cell")
data <- read_annovar(path = mut_path, col_keep = c(1,2,4,5,6,7,8,9,10,11:37,46:50))
colnames(data)[c(1:9,37:42)] <- col_name
sta <- data[,c(1:6,8,12,42)]
data_VAF <- split_FORMAT_col(dat = data, split = 41, name = 40, sep = ":") # slow
data_all <- cbind(data[,c(-40,-41)], data_VAF)

p <- statistic_single_mut(sta = sta, ref = "hg19", bar = "sd")
ggplot2::ggsave("basic_statistic_sd.pdf", p, width = 23, height = 10)
```

### 2. SVM

Using SVM to filter mutation:
```r
single_svm$mark <- ifelse(single_snp$edit=="no" & single_snp$dbsnp1=="yes", "True", "Test")
single_svm$mark <- ifelse(single_snp$edit=="yes" & single_snp$dbsnp1=="no", "False", single_svm$mark)
single_svm$mark <- ifelse(nchar(single_raw$Ref)==nchar(single_raw$Alt) & nchar(single_raw$Alt)>1, "False", single_svm$mark)
train <- single_svm[which(single_svm$mark %in% c("True","False")),]
test <- single_svm[which(single_svm$mark %in% c("Test")),]
# Fitting model
fit <- svm(mark ~ ., data = train)
# Predict Output
pred <- predict(fit, test[,c(-13)])
```

### 3. Infer genotype

From VAF to infer genotype:
```r
single <- mutation_filter(mut = single, min_celltype = 10, min_mut_per_cell = 10, min_cell_one_mut = 15, meta = meta)
matrix_REF <- col1_to_rowname(as.matrix(dcast(single[,c("cell","var","AD_1")], var~cell, value.var = "AD_1", fun.aggregate=sum)))
matrix_ALT <- col1_to_rowname(as.matrix(dcast(single[,c("cell","var","AD_2")], var~cell, value.var = "AD_2", fun.aggregate=sum)))
matrix_VAF <- make_VAF_matrix(matrix_REF = matrix_REF, matrix_ALT = matrix_ALT, fast_cutoff = 0.9)
```

### 4. Smooth genotype

Borrow information from neighbor cells:
```r
matrix_smooth <- make_smooth_matrix(matrix_VAF = matrix_VAF, matrix_ALT = matrix_ALT,matrix_REF = matrix_REF,meta = meta_select, sig = sig,
                                    ref = "hg19", 
                                    method = "manhattan",
                                    leader_min_depth = 20,
                                    neighbor_min_dis = 300,
                                    K_smooth = 0.2)
```

### 5. RPCA

Cell-mutation matrix denoising:
```r

matrix_rpca <- rpca(as.matrix(matrix_VAF), max.iter = 500)
```

### 6. Randomforest

Using RF to get clonal mutation:
```r
demo <- run_randomforest(matrix = matrix_VAF, meta = meta, train_ratio = 0.7, output_path = paste0(work_path,"/7_rf"))
importance <- demo[["importance"]]
var_keep <- rownames(importance)[which(importance$MeanDecreaseAccuracy>0&importance$MeanDecreaseGini>0)]
matrix_rf <- matrix_VAF[var_keep,]
```

### 7. Visualization

Visualization of heatmap and clonal structure:
```r

meta_select <- read.csv(paste0(work_path,"/1_basic/meta_select.csv"), row.names = 1)
matrix_VAF <- read.csv(paste0(work_path,"/7_rf/matrix_rf.csv"), row.names = 1)
cluster <- get_hclust(matrix = matrix_VAF, method = "manhattan", cut = 300)

p1 <- pheatmap(matrix_VAF[,meta_select$cell], annotation_col = meta_select[,c("cluster_new","sample","tissue")], show_rownames = F, show_colnames = F, 
               cluster_cols = F, cluster_rows = F, color = colorRampPalette(brewer.pal(9,"GnBu"))(50),
               clustering_distance_cols = "manhattan", annotation_colors = heatmap_color)

df <- meta_select %>% make_long(tissue, sample, cluster_new) 
df$node <- factor(df$node, levels = c(unique(meta_select$tissue),
                                      unique(meta_select$sample),
                                      unique(meta_select$cluster_new)))
p2 <- ggplot(df, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) +
  geom_sankey(flow.alpha = 0.7, node.color = "gray70", smooth = 5, width = 0.2) +
  geom_sankey_label(aes(fill = factor(node)), size = 5, color = "white") +
  scale_fill_manual(values = sankey_color)+
  theme_sankey(base_size = 18) +
  labs(x = NULL) +
  theme(legend.position = "none", plot.title = element_text(hjust = .5)) +
  ggtitle("tissue-sample-cluster")

table <- as.data.frame(table(meta_select$cluster_new, meta_select$sample))
p3 <- ggplot(table, aes(area = Freq, fill = Var1, label = Var1, subgroup = Var2)) +
  geom_treemap() +
  geom_treemap_text(fontface = "italic", colour = "white", place = "centre",grow = TRUE,alpha=.6)+
  geom_treemap_subgroup_border(color='white')+
  scale_fill_manual(values = sankey_color)

# save
graphics.off()
p <- plot_grid(p1$gtable, p2, p3, ncol = 3)
ggsave("demo.pdf",p, width = 30, height = 10)
```

## Notes
**Environment Requirements:**
**Local RStudio**: R ≥ 4.0, Seurat, scClone, e1071, randomForest, ggplot2, etc. (install via BiocManager/install.packages)
