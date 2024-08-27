# RARE pipelines
va**R**iant p**A**thogenicity p**RE**diction pipelines

A set easy-to-use, open source pipelines for prioritizing variants in Mendelian Diseases.

## Features
This repository offers a comprehensive pipeline designed to annotate and prioritize variants associated with Mendelian diseases. The pipeline leverages two key components: [Maverick](https://github.com/ZuchnerLab/Maverick) and [CADA](https://github.com/Chengyao-Peng/CADA). The pipeline only **support only GRCh38** genome.

Maverick (Mendelian Approach to Variant Effect pRedICtion built in Keras): Maverick is a state-of-the-art variant effect prediction model built using Keras. It employs transformer architectures to process a diverse set of input features, enabling it to accurately classify variants as benign, dominant pathogenic, or recessive pathogenic.

CADA (Case Annotations and Disease Dnnotations): is a powerful tool for identifying disease-causing genes in rare syndromes. It integrates disease-level annotations from the Human Phenotype Ontology (HPO) with clinical case-level data to construct a gene-phenotype association network. By applying network representation learning techniques, CADA effectively prioritizes genes based on their likelihood of involvement in the disease.

The whole pipeline is into a single Singularity container named **rare_pipes**.

## Installation
Prerequisite: Due to licensing restrictions, Annovar could not be directly included in the Singularity image. To use Maverick, you'll need to setup Annovar separately with the steps you find below.

```bash
# Step 1: Download Annovar
# Annovar executables are bound by both academic and commercial licenses, but can be downloaded (after registration) from https://www.openbioinformatics.org/annovar/annovar_download_form.php
# Users must obtain their own licenses for use of Annovar and then move the release archive to the a directory of ypur choise. 

# Step 2: Extract Annovar to a new annovar/ subdirectory.
tar -zxvf annovar.tar.gz

# 3: Configure the Gencode annotations for use with Annovar (GRCh38)
cd annovar/
./annotate_variation.pl -downdb -build hg38 seq humandb/hg38_seq/
# The following line of code may not be necessary with newer versions of Annovar. Use as needed. 
mv humandb/hg38_seq/chroms/* humandb/hg38_seq/
./annotate_variation.pl --build hg38 --downdb wgEncodeGencodeBasicV33 humandb/
./retrieve_seq_from_fasta.pl -format genericGene -seqdir humandb/hg38_seq/ -outfile humandb/hg38_wgEncodeGencodeBasicV33Mrna.fa humandb/hg38_wgEncodeGencodeBasicV33.txt 
cd ../
```

Next, to install the rare pipes Docker/Singularity image, run the following commands:
```bash
# 1. Install with Singularity and test it
singularity pull docker://gambalab/rare_pipes:1.0.0

# 2. Install with Docker and test it
docker pull gambalab/rare_pipes:1.0.0
 ```

