#!/bin/bash

#############################################
# Color Definitions                         #
#############################################
# Reset
 Color_Off=$'\033[0m'       # Text Reset

 # Regular Colors
 Black=$'\033[0;30m'        # Black
 Red=$'\033[0;31m'          # Red
 Green=$'\033[0;32m'        # Green
 Yellow=$'\033[0;33m'       # Yellow
 Blue=$'\033[0;34m'         # Blue
 Purple=$'\033[0;35m'       # Purple
 Cyan=$'\033[0;36m'         # Cyan
 White=$'\033[0;37m'        # White

print_info(){
 dt=$(date '+%d/%m/%Y %H:%M:%S')
 echo "[${Cyan}${dt}${Color_Off}] [${Green}info${Color_Off}] ${1}"
}

print_error(){
 dt=$(date '+%d/%m/%Y %H:%M:%S')
 echo "[${Cyan}${dt}${Color_Off}] [${Red}error${Color_Off}] ${1}"
}

################################################################################
# Help                                                                         #
################################################################################
Help()
{
   # Display Help
    echo "Run variant filtering."
    echo "Syntax: filter_variants.sh [h|i|o|s|r|b]"
        echo "options:"
        echo "-h     Print this Help."
        echo "-i     Path to the input vcf file. (Required)"
        echo "-o     Output directory. (Required)"
        echo "-s     Sample name. (Required)"
        echo "-r     Absolute path to the reference FASTA file. (Required)"
        echo "-b     Bed file of interested regions. (Optional)"
        echo
}


OPTIND=1         # Reset in case getopts has been used previously in the shell.
BED_file=""
while getopts ":hi:o:s:a:b:r:" option; do
   case $option in
      h) # display Help
         Help
         exit;;
      i)
         INPUT_VCF=${OPTARG}
         ((count++))
         ;;
      o)
         OUT_folder=${OPTARG}
         ((count++))
         ;;
      s)
         SAMPLE=${OPTARG}
         ((count++))
         ;;
      b)
         BED_file=${OPTARG}
         ;;
      r)
         REF_fasta=${OPTARG}
         ((count++))
         ;;
      :)
         print_error "Option -${OPTARG} requires an argument."
         exit 1
         ;;
     \?) # incorrect option
         print_error "Invalid Input option -${OPTARG}"
         exit;;
   esac
done


# Check the number of input args correspond
if [[ ${count} == 0 ]]; then
   print_error "No arguments in input, please see the help (-h)"
   exit
elif [[ ${count} -lt 4 ]]; then
      print_error "Missing some input arguments, please see the help (-h)"
      exit
fi

##################################################
# Set vars
eval "$(conda shell.bash hook)"
conda activate bio
snpSIFT="/opt/snpEff/SnpSift.jar"
snpEFF="/opt/snpEff/snpEff.jar"
dbNSFP=/opt/snpEff/data/dbNSFP4.9_small.txt.gz
##################################################

mkdir -p ${OUT_folder}
TMP_DIR=$(mktemp -d --tmpdir=${OUT_folder})


if [ "${BED_file}" == "" ]; then
	print_info "Starting Annotation (BED file not provided).."
	bcftools view -f PASS ${INPUT_VCF} |
	bcftools norm -m-both | \
	bcftools norm -f ${REF_fasta} | \
	java -Xmx8g -jar ${snpEFF} ann \
		-nodownload \
		-noStats \
		-no-downstream \
		-no-intergenic \
		-no-upstream \
		-no-utr \
		GRCh38.p14 - | \
	java -jar ${snpSIFT} filter "( ANN[*].IMPACT = 'HIGH' ) | (ANN[*].IMPACT = 'MODERATE')" | \
	bgzip -c  > ${TMP_DIR}/${SAMPLE}.ontarget.vcf.gz
	tabix -p vcf ${TMP_DIR}/${SAMPLE}.ontarget.vcf.gz
else
	if [ -f ${BED_file} ]; then
		print_info "Starting Annotation (BED file provided).."
		bcftools view -f PASS ${INPUT_VCF} | \
		bedtools intersect -header -a - -b ${BED_file} | \
		bcftools norm -m-both | \
		bcftools norm -f ${REF_fasta} | \
		java -Xmx8g -jar ${snpEFF} ann \
			-nodownload \
			-noStats \
			-no-downstream \
			-no-intergenic \
			-no-upstream \
			-no-utr \
			GRCh38.p14 - | \
		java -jar ${snpSIFT} filter "( ANN[*].IMPACT = 'HIGH' ) | (ANN[*].IMPACT = 'MODERATE')" | \
		bgzip -c  > ${TMP_DIR}/${SAMPLE}.ontarget.vcf.gz
		tabix -p vcf ${TMP_DIR}/${SAMPLE}.ontarget.vcf.gz
	else
		print_error "Input BED file not found!"
	fi
fi

print_info "dbNSFP annotation.."
java -jar ${snpSIFT} dbnsfp \
			-db ${dbNSFP} \
			-v ${TMP_DIR}/${SAMPLE}.ontarget.vcf.gz \
			-f AlphaMissense_score,EVE_score,ESM1b_score,REVEL_score,gnomAD_exomes_AF,gnomAD_genomes_AF,clinvar_id,clinvar_clnsig,clinvar_trait,clinvar_review,clinvar_hgvs,clinvar_var_source,clinvar_MedGen_id,clinvar_OMIM_id,clinvar_Orphanet_id \
			| java -jar ${snpSIFT} varType - \
			| bgzip -c > ${TMP_DIR}/${SAMPLE}.ontarget.annotated.dedup.vcf.gz
tabix -f -p vcf ${TMP_DIR}/${SAMPLE}.ontarget.annotated.dedup.vcf.gz


print_info "Variant Filtering.."
zcat ${TMP_DIR}/${SAMPLE}.ontarget.annotated.dedup.vcf.gz | \
	java -jar ${snpSIFT} filter "(VARTYPE = 'SNP')" | \
	java -jar ${snpSIFT} filter "(dbNSFP_gnomAD_genomes_AF < 0.05 | dbNSFP_gnomAD_exomes_AF < 0.05 | (na dbNSFP_gnomAD_genomes_AF) | (na dbNSFP_gnomAD_exomes_AF) )" | \
	java -jar ${snpSIFT} filter "(dbNSFP_AlphaMissense_score > 0.5 | dbNSFP_EVE_score > 0.5 | dbNSFP_REVEL_score > 0.5 | dbNSFP_ESM1b_score < -7.5 )" | \
	bgzip -c > ${TMP_DIR}/${SAMPLE}.ontarget.annotated.dedup.SNPs.filtered.vcf.gz
tabix -p vcf ${TMP_DIR}/${SAMPLE}.ontarget.annotated.dedup.SNPs.filtered.vcf.gz

zcat ${TMP_DIR}/${SAMPLE}.ontarget.annotated.dedup.vcf.gz | \
	java -jar ${snpSIFT} filter "(VARTYPE != 'SNP')" | \
	bgzip -c > ${TMP_DIR}/${SAMPLE}.ontarget.annotated.dedup.INDELs.filtered.vcf.gz
tabix -p vcf ${TMP_DIR}/${SAMPLE}.ontarget.annotated.dedup.INDELs.filtered.vcf.gz

zcat ${TMP_DIR}/${SAMPLE}.ontarget.annotated.dedup.vcf.gz | \
	java -jar ${snpSIFT} filter "ANN[*].EFFECT has 'stop_lost'" | \
	bgzip -c > ${TMP_DIR}/${SAMPLE}.ontarget.annotated.dedup.SL.filtered.vcf.gz
tabix -p vcf ${TMP_DIR}/${SAMPLE}.ontarget.annotated.dedup.SL.filtered.vcf.gz

zcat ${TMP_DIR}/${SAMPLE}.ontarget.annotated.dedup.vcf.gz | \
	java -jar ${snpSIFT} filter "ANN[*].EFFECT has 'stop_gained'" | \
	bgzip -c > ${TMP_DIR}/${SAMPLE}.ontarget.annotated.dedup.SG.filtered.vcf.gz
tabix -p vcf ${TMP_DIR}/${SAMPLE}.ontarget.annotated.dedup.SG.filtered.vcf.gz

zcat ${TMP_DIR}/${SAMPLE}.ontarget.annotated.dedup.vcf.gz | \
	java -jar ${snpSIFT} filter "ANN[*].EFFECT has 'start_lost'" | \
	bgzip -c > ${TMP_DIR}/${SAMPLE}.ontarget.annotated.dedup.SAL.filtered.vcf.gz
tabix -p vcf ${TMP_DIR}/${SAMPLE}.ontarget.annotated.dedup.SAL.filtered.vcf.gz

print_info "Merge files.."
ionice -c 3 bcftools concat -a \
        --threads 4 \
        "${TMP_DIR}/${SAMPLE}.ontarget.annotated.dedup.SNPs.filtered.vcf.gz" \
        "${TMP_DIR}/${SAMPLE}.ontarget.annotated.dedup.INDELs.filtered.vcf.gz" \
        "${TMP_DIR}/${SAMPLE}.ontarget.annotated.dedup.SL.filtered.vcf.gz" \
        "${TMP_DIR}/${SAMPLE}.ontarget.annotated.dedup.SG.filtered.vcf.gz" \
        "${TMP_DIR}/${SAMPLE}.ontarget.annotated.dedup.SAL.filtered.vcf.gz" | \
        bcftools sort - -o - | \
        bcftools norm -d all | \
        bgzip -c > "${OUT_folder}/${SAMPLE}.ontarget.annotated.dedup.filtered.vcf.gz"
tabix -p vcf "${OUT_folder}/${SAMPLE}.ontarget.annotated.dedup.filtered.vcf.gz"


#zcat ${TMP_DIR}/${SAMPLE}.ontarget.annotated.dedup.vcf.gz | \
#	java -jar ${snpSIFT} filter "!(na dbNSFP_clinvar_id)" | \
#	bgzip -c > "${OUT_folder}/${SAMPLE}.ontarget.annotated.dedup.filtered.Clinvar.vcf.gz"
#tabix -p vcf "${OUT_folder}/${SAMPLE}.ontarget.annotated.dedup.filtered.Clinvar.vcf.gz"

rm -rf ${TMP_DIR}
print_info "DONE!!"