#!/bin/bash

## Example run: ./run_happy_vcfeval.sh UDP10618 UDP10618 /data/markellocj/ t146
#   COHORT_NAME="UDP10618"
#   SAMPLE_NAME="UDP11320"
#   WORK_DIR="/data/markellocj/"
#   VG_VERSION="t146"
#   PACKAGE_DIR="/data/markellocj/vg_workaround_pipeline_pkg"
#
#   cd ${WORK_DIR}/${COHORT_NAME}_run/output_vg_${VG_VERSION}/output_${SAMPLE_NAME}_vg_${VG_VERSION}
#   rm -fr ${SAMPLE_NAME}_vcf_compare
#   sbatch --partition=norm --cpus-per-task=32 --mem=100g --time=8:00:00 ${PACKAGE_DIR}/run_happy_vcfeval.sh ${COHORT_NAME} ${SAMPLE_NAME} ${WORK_DIR} ${VG_VERSION}


COHORT_NAME=$1
SAMPLE_NAME=$2
WORK_DIR=$3
VG_VERSION=$4
SAMPLE_VCF=$5

SAMPLE_VCF_FILENAME=${SAMPLE_VCF##*/}
SAMPLE_VCF_BASENAME=${SAMPLE_VCF_FILENAME%.*.*}
VCF_ANALYSIS_DIR="${WORK_DIR}/${COHORT_NAME}_run/output_vg_${VG_VERSION}/output_${SAMPLE_NAME}_vg_${VG_VERSION}/${SAMPLE_NAME}_vcf_compare"

SDF_DIR="${HOME}/sdf_references/hs37d5.fa.sdf"
UDN_BASELINE_FILE="${HOME}/${COHORT_NAME}_run/output_vg_${VG_VERSION}/output_${SAMPLE_NAME}_vg_${VG_VERSION}/${SAMPLE_NAME}_vcf_compare/${UDN_VCF_NAME}.new.vcf.gz"
SNP_BASELINE_FILE="${HOME}/${COHORT_NAME}_run/output_vg_${VG_VERSION}/output_${SAMPLE_NAME}_vg_${VG_VERSION}/${SAMPLE_NAME}_vcf_compare/${SAMPLE_NAME}_only_SNPchip.new.vcf.gz"
VCFEVAL_UDN_vs_SNP_OUTPUT="${HOME}/${COHORT_NAME}_run/output_vg_${VG_VERSION}/output_${SAMPLE_NAME}_vg_${VG_VERSION}/${SAMPLE_NAME}_vcf_compare/vcfeval_output_SAMPLE_vcf_TRUTH"

module load hap.py
module load samtools

# Preprocess UDN vcf
cp ${SAMPLE_VCF} ${VCF_ANALYSIS_DIR}
cd ${VCF_ANALYSIS_DIR}
bgzip -d ${SAMPLE_VCF_FILENAME}
awk '{if($0 !~ /^#|^chr/) print "chr"$0; else print $0}' ${SAMPLE_VCF_BASENAME}.vcf > ${SAMPLE_VCF_BASENAME}.new.vcf
bgzip ${SAMPLE_VCF_BASENAME}.new.vcf
tabix -f -p vcf ${SAMPLE_VCF_BASENAME}.new.vcf.gz
bgzip ${SAMPLE_VCF_BASENAME}.vcf
tabix -f -p vcf ${SAMPLE_VCF_BASENAME}.vcf.gz


## Run hap.py SAMPLE VCF against TRUTH VCF
cd ${VCF_ANALYSIS_DIR}
hap.py ${SAMPLE_NAME}_only_SNPchip.new.vcf.gz ${SAMPLE_VCF_BASENAME}.new.vcf.gz -o happy_${SAMPLE_NAME}_UDN_snpchip --threads 32 2> happy_output.txt


module load singularity

## Run vcfeval SAMPLE VCF against TRUTH VCF
cd ${WORK_DIR}
rm -fr "${VCF_ANALYSIS_DIR}/vcfeval_output_SAMPLE_vcf_TRUTH"
singularity exec -H ${PWD}:${HOME}  docker://realtimegenomics/rtg-tools rtg vcfeval -T 32 -b ${SNP_BASELINE_FILE} -c ${UDN_BASELINE_FILE} -t ${SDF_DIR} -o ${VCFEVAL_UDN_vs_SNP_OUTPUT}

