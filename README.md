# RARE pipelines
va**R**iant p**A**thogenicity p**RE**diction pipelines

A set easy-to-use, open source pipelines for prioritizing variants in Mendelian Diseases.

## Features
This repository offers a comprehensive pipeline designed to annotate and prioritize variants associated with Mendelian diseases. The pipeline leverages two key components: [Maverick (Mendelian Approach to Variant Effect pRedICtion built in Keras)](https://github.com/ZuchnerLab/Maverick) and [CADA (Case Annotations and Disease Annotations)](https://github.com/Chengyao-Peng/CADA). The pipeline **support only GRCh38** genome.

* **Maverick:** Maverick is a state-of-the-art variant effect prediction model built using Keras. It employs transformer architectures to process a diverse set of input features, enabling it to accurately classify variants as benign, dominant pathogenic, or recessive pathogenic.

* **CADA**: is a powerful tool for identifying disease-causing genes in rare syndromes. It integrates disease-level annotations from the Human Phenotype Ontology (HPO) with clinical case-level data to construct a gene-phenotype association network. By applying network representation learning techniques, CADA effectively prioritizes genes based on their likelihood of involvement in the disease.

The whole pipeline is into a single Singularity container named **rare_pipes**.

## Installation
Prerequisite: Due to licensing restrictions, Annovar could not be directly included in the Singularity image. To use Maverick, you'll need to setup Annovar separately with the steps you find below.

```bash
# Step 1: Download Annovar
# Annovar executables are bound by both academic and commercial licenses,
# but can be downloaded (after registration) from https://www.openbioinformatics.org/annovar/annovar_download_form.php
#
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

## Pipeline Overview
The pipeline begins by leveraging [snpEFF](https://pcingola.github.io/SnpEff/) and [dbSNP 4.9](https://sites.google.com/site/jpopgen/dbNSFP) to identify a preliminary set of noteworthy mutations. Specifically, it focuses on the following criteria:

* **SNPs:** Only coding SNPs with an allele frequency (AF) less than 0.05 and a damaging prediction score greater than 0.5 for at least one of the **Alpha Missense**, **EVE**, or **REVEL** predictors, or with an **ESM1b** prediction score less than -7.5, are retained.
* **INDELs:** All INDELs are included.
* **Stop Gain, Stop Loss, and Start Loss Variants:** These variants are also retained.

After filtering based on these criteria, the pipeline proceeds to analyze the selected mutations using Maverick. Finally, the results from Maverick are integrated with CADA results (if HPO codes are provided) and data collected from [OMIM](https://www.omim.org/) and [Orphanet](https://www.orpha.net/) databases.

## Example of use

* **Input:**
	* **Not annotated** VCF file produced by **Deepvariant**. Althoug any other variant caller shuold also work.
	* The absolute path to output folder where store the results.
	* Sample Name.
	* The absolute path to the reference FASTA file.
	* The absolute path to Annovar.
	* Optional bed file on which cut mutatins.

* **Output:**
	* Filtered vcf file used as input of Maverick.
	* TSV file with scored variants.

The pipeline is coded in the ```run_rare.sh``` script and can be run in the following way:

```bash
# Let's define a honey_pipe_exec variable to excec the several commands 
RARE_exec="singularity exec --bind /usr/lib/locale/ path/to/rare_pipes_1.0.0.sif"

# Let's see the help
${RARE_exec} run_rare.sh -h
```
```
Run RARE pipeline for GRCh38.
Syntax: run_rare.sh [h|i|o|s|r|a|m|b|p|d|t]
options:
-h     Print this Help.
-i     Input vcf file. (Required)
-o     Output directory. (Required)
-s     Sample name. (Required)
-r     Absolute path to the reference FASTA file. (Required)
-a     Absolute path to Annovar. (Required)
-m     Model type, e.g. full or light. Default full. (Optional)
-p     A string of comma-separated HPO terms describing the patient, e.g. HP:0000573,HP:0001102,HP:0003115 (Optional)
-d     Vcf file containing mutation to discard. Must be indexed. (Optional)
-t     Trio vcf file annotated with rtg tool. Must be indexed and sample name be contained into it. (Optional)
```

So a typical case of use will be something like this:
```bash
${RARE_exec} \
    run_rare.sh \
    -i /path/to/input/folder/with/vcf.gz \
    -o /path/to/output_folder \
    -s sample_name \
    -r /path/to/reference.fa \
    -p HP:0002354,HP:0000544
```

## Pipeline Output Files:
After the pipeline execution, results will be written to the /path/to/output_folder. You will find three primary files:

* ```${sample_name}.UD_AN001_P.ontarget.annotated.dedup.filtered.vcf.gz```: The primary VCF file containing filtered, annotated, and deduplicated variants.
* ```${sample_name}.UD_AN001_P.ontarget.annotated.dedup.filtered.vcf.gz.tbi```: The tabix index for the VCF file, enabling efficient querying.
* ```${sample_name}.ranked.mutations.tsv```: A tab-separated file containing mutations ranked by pathogenicity score.

**Mutation Ranking** in ```${sample_name}.ranked.mutations.tsv```:
*	If HPO terms are provided: Mutations will be ranked according to the ```RARE.score``` column, which is the average of the ```Maverick.Score``` and ```CADA.score``` columns.
*	If HPO terms are not provided: Mutations will be ranked solely by the ```Maverick.Score```.

## Explanation of ```ranked.mutations.tsv``` file
Below you can find an explations of the most important columns in the ```${sample_name}.ranked.mutations.tsv```

**Core VCF Information**:
*	**CHROM**: Chromosome on which the variant is located.
*	**POS**: Base-pair position of the variant on the chromosome (1-based indexing).
*	**REF**: Reference allele at the variant position.
*	**ALT**: Alternative allele(s) at the variant position.
*	**QUAL**: Phred-scaled quality score for the variant call.
*	**VarType**: Type of variant (e.g., missense, stop gain ecc.).


**Genotype Information**:
*	**Symbol**: Gene symbol associated with the variant (if applicable).
*	**AD**: Allelic depth for each allele (separated by commas).
*	**DP**: Total read depth at the variant position.
*	**GQ**: Genotype quality score.
*	**GT**: Genotype call for the sample (e.g., 0/1, 1/1).
*	**PL**: Phred-scaled genotype likelihoods.
*	**VAF**: Variant Allele Frequency


**Functional Impact Prediction**:
*	**AlphaMissense**: Score for the predicted missense effect of the variant.
*	**ESM1b**: Score for the predicted effect of with ESM.
*	**EVE**: Evolutionary model of Variant Effect score.
*	**REVEL**: Rare Exome Variant Ensemble Learner score.


**Disease Association Annotations**:
*	**Orphanet.id**: Orphanet identifier for a disease associated with the variant (if applicable).
*	**Orphanet.Phenotype**: Orphanet disease phenotype associated with the variant (if applicable).
*	**Orphanet.inheritance**: Orphanet inheritance mode for the associated disease (if applicable).
*	**OMIM.id**: OMIM identifier for a Mendelian disease associated with the variant (if applicable).
*	**OMIM.title**: Title of the associated OMIM disease.
*	**OMIM.Phenotype**: OMIM disease phenotype associated with the variant (if applicable).
*	**clinvar_MedGen_id**: MedGen identifier from ClinVar.
*	**clinvar_OMIM_id**: OMIM identifier from ClinVar.
*	**clinvar_Orphanet_id**: Orphanet identifier from ClinVar.
*	**clinvar_clnsig**: Clinical significance classification from ClinVar.
*	**clinvar_hgvs**: HGVS notation for the variant from ClinVar.
*	**clinvar_id**: ClinVar accession identifier.
*	**clinvar_review**: Review status of the variant in ClinVar.
*	**clinvar_trait**: Associated trait(s) from ClinVar.
*	**clinvar_var_source**: Source of the variant information in ClinVar.


**GenomeAD annotation**:
*	**gnomAD_exomes_AF**: Allele frequency in gnomAD exomes dataset.
*	**gnomAD_genomes_AF**: Allele frequency in gnomAD genomes dataset.


## Acknowledgements
Rare pipes happily makes use of many open source packages. We would like to specifically call out a few key ones:
*	[bcftools](https://github.com/samtools/bcftools)
*   [bedtools](https://github.com/arq5x/bedtools2)
*	[tabix](https://github.com/samtools/tabix)
*	[Maverick](https://github.com/ZuchnerLab/Maverick)
*	[CADA](https://github.com/Chengyao-Peng/CADA)
*	[vcf2tsvpy](https://github.com/sigven/vcf2tsvpy)
*	[snpEFF](https://pcingola.github.io/SnpEff/)
*	[dbNSFP](https://sites.google.com/site/jpopgen/dbNSFP)

We thank all of the developers and contributors to these packages for their work.
