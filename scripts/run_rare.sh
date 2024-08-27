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
    echo "Run RARE pipeline for GRCh38."
    echo "Syntax: run_rare.sh [h|i|o|s|m|a]"
        echo "options:"
        echo "-h     Print this Help."
        echo "-i     Input vcf file. (Required)"
        echo "-o     Output directory. (Required)"
        echo "-s     Sample name. (Required)"
        echo "-r     Absolute path to the reference FASTA file. (Required)"
        echo "-a     Absolute path to Annovar. (Required)"
        echo "-m     Model type, e.g. full or light. Default full. (Optional)"
        echo "-b     Bed file of interested regions. (Optional)"
        echo
}

# A POSIX variable
OPTIND=1         # Reset in case getopts has been used previously in the shell.
MODEL="full"
BED_file=""
while getopts ":hi:o:s:a:r:b:" option; do
   case $option in
      h) # display Help
         Help
         exit;;
      i)
         fileName=${OPTARG}
         ((count++))
         ;;
      o)
         OUT_DIR=${OPTARG}
         ((count++))
         ;;
      s)
         SAMPLE=${OPTARG}
         ((count++))
         ;;
      m)
         MODEL=${OPTARG}
         ;;
      a)
         ANNOVAR_PATH=${OPTARG}
         ((count++))
         ;;
      r)
         REF_fasta=${OPTARG}
         ((count++))
         ;;
      b)
         BED_file=${OPTARG}
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
elif [[ ${count} -lt 5 ]]; then
      print_error "Missing some input arguments, please see the help (-h)"
      exit
fi

##################################################
# Set vars
eval "$(conda shell.bash hook)"
conda activate maverick
MAVERICK_ROOT="/opt/maverick"
cd ${MAVERICK_ROOT}
export PATH=${ANNOVAR_PATH}:$PATH
#################################################


mkdir -p ${OUT_DIR}
TMP_DIR=$(mktemp -d --tmpdir=${OUT_DIR})

# Filter variants
print_info "Filter Variants Step: Get putative mutations.."
if [ "${BED_file}" != "" ]; then
   opt_args="-b ${BED_file}" 
fi
filter_variants.sh \
        -i ${fileName} \
        -o ${OUT_DIR} \
        -s ${SAMPLE} \
        -r ${REF_fasta} ${opt_args}

# update vcf filname
fileName=${OUT_DIR}/${SAMPLE}.ontarget.annotated.dedup.filtered.vcf.gz

# process variants with annovar
print_info "Maverick Step 1/3: Get coding changes with Annovar.."
# remove chr prefix if present
if (file ${fileName} | grep -q gzip ) ; then
    zcat ${fileName} | sed 's/^chr//' > "${TMP_DIR}/${SAMPLE}.vcf"
 else
 	sed 's/^chr//' ${fileName} > "${TMP_DIR}/${SAMPLE}.vcf"
fi


convert2annovar.pl -format vcf4 --keepindelref "${TMP_DIR}/${SAMPLE}.vcf" > ${TMP_DIR}/${SAMPLE}.avinput

annotate_variation.pl -dbtype wgEncodeGencodeBasicV33 \
					  -buildver hg38 \
					  --exonicsplicing ${TMP_DIR}/${SAMPLE}.avinput \
					  ${ANNOVAR_PATH}/humandb/

# if there are no scorable variants, end early
SCORABLEVARIANTS=$(cat ${TMP_DIR}/${SAMPLE}.avinput.exonic_variant_function | wc -l || true)
if [[ ${SCORABLEVARIANTS} -eq 0 ]]; then print_error "No scorable variants found"; exit 0; fi

coding_change.pl ${TMP_DIR}/${SAMPLE}.avinput.exonic_variant_function \
				 ${ANNOVAR_PATH}/humandb/hg38_wgEncodeGencodeBasicV33.txt \
				 ${ANNOVAR_PATH}/humandb/hg38_wgEncodeGencodeBasicV33Mrna.fa \
				 --includesnp --onlyAltering --alltranscript \
				 > "${TMP_DIR}/${SAMPLE}.coding_changes.txt"

# select transcript
print_info "Maverick Step 2/4: Select transcript.."
python Maverick/InferenceScripts/groomAnnovarOutput_GRCh38.py --inputBase=${TMP_DIR}/${SAMPLE}

# if there are no scorable variants, end early
SCORABLEVARIANTS=$(cat ${TMP_DIR}/${SAMPLE}.groomed.txt | wc -l || true)
if [[ ${SCORABLEVARIANTS} -lt 2 ]]; then print_error "No scorable variants found"; exit 0; fi

# add on annotations
print_info "Maverick Step 3/4: Merge on annotations.."
python Maverick/InferenceScripts/annotateVariants_GRCh38.py --inputBase=${TMP_DIR}/${SAMPLE}

# run variants through each of the models
if [ "${MODEL}" == "full" ]; then
	print_info "Maverick Step 4/4: Score variants with the 8 models .."
	python Maverick/InferenceScripts/runModels_GRCh38.py --inputBase=${TMP_DIR}/${SAMPLE}
	python Maverick/InferenceScripts/rankVariants_GRCh38.py --inputBase=${TMP_DIR}/${SAMPLE}
else
	print_info "Maverick Step 4/4: Score variants with architecture 1 model 1 only .."
	python Maverick/InferenceScripts/runLiteModel_GRCh38.py --inputBase=${TMP_DIR}/${SAMPLE}
	python Maverick/InferenceScripts/rankVariantsLite_GRCh38.py --inputBase=${TMP_DIR}/${SAMPLE}
fi

print_info "Post process Maverick results.."
conda activate bio
vcf2tsvpy --input_vcf ${fileName} --out_tsv ${TMP_DIR}/${SAMPLE}.vcf.tsv
postprocess_results.R  ${TMP_DIR}/${SAMPLE}.vcf.tsv ${TMP_DIR}/${SAMPLE}.finalScores.txt ${OUT_DIR}/${SAMPLE}.ranked.mutations.tsv

rm -rf ${TMP_DIR}
print_info "Done"
