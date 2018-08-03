#!/bin/bash

COHORT_NAME=$1
COHORT_WORK_DIR=$2
UDPBINFO_DIR=$3

# Example INPUT:
#COHORT_NAME="UDP10618"
#COHORT_WORK_DIR="/data/markellocj/UDP_PEDIGREE_EXPERIMENT_1/UDP10618_run"
#UDPBINFO_DIR="Udpbinfo/usr/markellocj"

# Extract cohort sample names
COHORT_SAMPLE_NAMES=($(ls ${COHORT_WORK_DIR} | grep '^UDP' | awk -F'_' '{print $1}' | sort | uniq))

# COPY Cohort BAMs to UDPBINFO Workdirectory
for sample_name in ${COHORT_SAMPLE_NAMES[@]}
do
    cp -r ${COHORT_WORK_DIR}/${sample_name}_workdir/${sample_name}_surjected_bams/ /data/${UDPBINFO_DIR}
    SAMPLE_DRAGEN_WORK_DIR="/staging/markellocj/${sample_name}"
    TMP_DIR="/staging/markellocj/tmp"
    ssh -t markellocj@helix.nih.gov ssh udpdragen01.nhgri.nih.gov "mkdir -p ${SAMPLE_DRAGEN_WORK_DIR}" && \
    ssh -t markellocj@helix.nih.gov ssh udpdragen01.nhgri.nih.gov "mkdir -p ${TMP_DIR}" && \
    ssh -t markellocj@helix.nih.gov ssh udpdragen01.nhgri.nih.gov "dragen -f -r /staging/hs37d5_v7 -b /staging/helix/${UDPBINFO_DIR}/${sample_name}_surjected_bams/${sample_name}_merged.sorted.dupmarked.reordered.bam --verbose --bin_memory=50000000000 --enable-map-align false --enable-variant-caller true --pair-by-name=true --vc-emit-ref-confidence GVCF --vc-sample-name ${sample_name} --intermediate-results-dir ${TMP_DIR} --output-directory ${SAMPLE_DRAGEN_WORK_DIR} --output-file-prefix ${sample_name}_dragen_genotypecall_wgs > ${COHORT_WORK_DIR}/${sample_name}_workdir/dragen_call_${sample_name}_wgs.gvcf.stdout 2> ${COHORT_WORK_DIR}/${sample_name}_workdir/dragen_call_${sample_name}_wgs.gvcf.stderr"
done


# RUN Dragen joint genotyping using cohort GVCFs
JOINT_GENOTYPE_DRAGEN_WORK_DIR="/staging/markellocj/output_cohort_joint_call_${COHORT_NAME}"
ssh -t markellocj@helix.nih.gov ssh udpdragen01.nhgri.nih.gov "mkdir -p ${JOINT_GENOTYPE_DRAGEN_WORK_DIR}/sample_gvcfs" && \
ssh -t markellocj@helix.nih.gov ssh udpdragen01.nhgri.nih.gov "mkdir -p ${JOINT_GENOTYPE_DRAGEN_WORK_DIR}/${COHORT_NAME}_cohort_joint_genotyped_output" && \
ssh -t markellocj@helix.nih.gov ssh udpdragen01.nhgri.nih.gov "touch ${JOINT_GENOTYPE_DRAGEN_WORK_DIR}/gvcf_list.txt"

for sample_name in ${COHORT_SAMPLE_NAMES[@]}
do
    SAMPLE_DRAGEN_WORK_DIR="/staging/markellocj/${sample_name}"
    ssh -t markellocj@helix.nih.gov ssh udpdragen01.nhgri.nih.gov "cp ${SAMPLE_DRAGEN_WORK_DIR}/${sample_name}_dragen_genotypecall_wgs.gvcf.gz ${JOINT_GENOTYPE_DRAGEN_WORK_DIR}/sample_gvcfs" && \
    ssh -t markellocj@helix.nih.gov ssh udpdragen01.nhgri.nih.gov "echo \"${JOINT_GENOTYPE_DRAGEN_WORK_DIR}/sample_gvcfs/${sample_name}_dragen_genotypecall_wgs.gvcf.gz\" \>\> ${JOINT_GENOTYPE_DRAGEN_WORK_DIR}/gvcf_list.txt"
done

ssh -t markellocj@helix.nih.gov ssh udpdragen01.nhgri.nih.gov "dragen -f -r /staging/hs37d5_v7 --enable-joint-genotyping true --intermediate-results-dir ${TMP_DIR} --output-directory ${JOINT_GENOTYPE_DRAGEN_WORK_DIR}/${COHORT_NAME}_cohort_joint_genotyped_output --output-file-prefix cohort_joint_genotyped_${COHORT_NAME} --variant-list ${JOINT_GENOTYPE_DRAGEN_WORK_DIR}/gvcf_list.txt > ${COHORT_WORK_DIR}/dragen_cohort_joint_call_${COHORT_NAME}.stdout 2> ${COHORT_WORK_DIR}/dragen_cohort_joint_call_${COHORT_NAME}.stderr"


# Copy final results to helix directory and do work directory cleanup
mkdir /data/${UDPBINFO_DIR}/output_cohort_joint_call_${COHORT_NAME} && chmod ug+rw -R /data/${UDPBINFO_DIR}/output_cohort_joint_call_${COHORT_NAME} && \
ssh -t markellocj@helix.nih.gov ssh udpdragen01.nhgri.nih.gov "cp -R ${JOINT_GENOTYPE_DRAGEN_WORK_DIR}/${COHORT_NAME}_cohort_joint_genotyped_output/. /staging/helix/${UDPBINFO_DIR}/output_cohort_joint_call_${COHORT_NAME}" && \
ssh -t markellocj@helix.nih.gov ssh udpdragen01.nhgri.nih.gov "rm -fr ${JOINT_GENOTYPE_DRAGEN_WORK_DIR}/"
for sample_name in ${COHORT_SAMPLE_NAMES[@]}
do
    SAMPLE_DRAGEN_WORK_DIR="/staging/markellocj/${sample_name}"
    ssh -t markellocj@helix.nih.gov ssh udpdragen01.nhgri.nih.gov "rm -fr ${SAMPLE_DRAGEN_WORK_DIR}"
done
for sample_name in ${COHORT_SAMPLE_NAMES[@]}
do
    rm -fr /data/${UDPBINFO_DIR}/${sample_name}_surjected_bams/
done

# Run snpEff on the final joint-genotyped vcf
module load snpEff
mv /data/${UDPBINFO_DIR}/output_cohort_joint_call_${COHORT_NAME} ${COHORT_WORK_DIR}/output_cohort_joint_call_${COHORT_NAME}
time snpEff -Xmx20g -i VCF -o VCF -noLof -noHgvs -formatEff -classic GRCh37.75 ${COHORT_WORK_DIR}/output_cohort_joint_call_${COHORT_NAME}/cohort_joint_genotyped_${COHORT_NAME}.vcf.gz > ${COHORT_WORK_DIR}/output_cohort_joint_call_${COHORT_NAME}/cohort_joint_genotyped_${COHORT_NAME}.snpeff.vcf


