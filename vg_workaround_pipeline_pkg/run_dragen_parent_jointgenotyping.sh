#!/bin/bash

COHORT_NAME=$1
PARENT_1_SAMPLE_NAME=$2
PARENT_2_SAMPLE_NAME=$3
COHORT_WORK_DIR=$4
UDPBINFO_DIR=$5

# Example INPUT:
#COHORT_NAME="UDP10618"
#PARENT_1_SAMPLE_NAME="UDP11319"
#PARENT_2_SAMPLE_NAME="UDP11320"
#COHORT_WORK_DIR="/data/markellocj/UDP_PEDIGREE_EXPERIMENT_1/UDP10618_run"
#UDPBINFO_DIR="Udpbinfo/usr/markellocj"

# COPY Parental BAMs to UDPBINFO Workdirectory
cp -r ${COHORT_WORK_DIR}/${PARENT_1_SAMPLE_NAME}_workdir/${PARENT_1_SAMPLE_NAME}_surjected_bams/ /data/${UDPBINFO_DIR}
cp -r ${COHORT_WORK_DIR}/${PARENT_2_SAMPLE_NAME}_workdir/${PARENT_2_SAMPLE_NAME}_surjected_bams/ /data/${UDPBINFO_DIR}

PARENT_1_DRAGEN_WORK_DIR="/staging/markellocj/${PARENT_1_SAMPLE_NAME}"
PARENT_2_DRAGEN_WORK_DIR="/staging/markellocj/${PARENT_2_SAMPLE_NAME}"
TMP_DIR="/staging/markellocj/tmp"

# RUN GVCF genotyping on parents through Dragen
ssh -t markellocj@helix.nih.gov ssh udpdragen01.nhgri.nih.gov "mkdir -p ${PARENT_1_DRAGEN_WORK_DIR}" && \
ssh -t markellocj@helix.nih.gov ssh udpdragen01.nhgri.nih.gov "mkdir -p ${PARENT_2_DRAGEN_WORK_DIR}" && \
ssh -t markellocj@helix.nih.gov ssh udpdragen01.nhgri.nih.gov "mkdir -p ${TMP_DIR}" && \
ssh -t markellocj@helix.nih.gov ssh udpdragen01.nhgri.nih.gov "dragen -f -r /staging/hs37d5_v7 -b /staging/helix/${UDPBINFO_DIR}/${PARENT_1_SAMPLE_NAME}_surjected_bams/${PARENT_1_SAMPLE_NAME}_merged.sorted.dupmarked.reordered.bam --verbose --bin_memory=50000000000 --enable-map-align false --enable-variant-caller true --pair-by-name=true --vc-emit-ref-confidence GVCF --vc-sample-name ${PARENT_1_SAMPLE_NAME} --intermediate-results-dir ${TMP_DIR} --output-directory ${PARENT_1_DRAGEN_WORK_DIR} --output-file-prefix ${PARENT_1_SAMPLE_NAME}_dragen_genotypecall_wgs > ${COHORT_WORK_DIR}/${PARENT_1_SAMPLE_NAME}_workdir/dragen_call_${PARENT_1_SAMPLE_NAME}_wgs.gvcf.stdout 2> ${COHORT_WORK_DIR}/${PARENT_1_SAMPLE_NAME}_workdir/dragen_call_${PARENT_1_SAMPLE_NAME}_wgs.gvcf.stderr" && \
ssh -t markellocj@helix.nih.gov ssh udpdragen01.nhgri.nih.gov "dragen -f -r /staging/hs37d5_v7 -b /staging/helix/${UDPBINFO_DIR}/${PARENT_2_SAMPLE_NAME}_surjected_bams/${PARENT_2_SAMPLE_NAME}_merged.sorted.dupmarked.reordered.bam --verbose --bin_memory=50000000000 --enable-map-align false --enable-variant-caller true --pair-by-name=true --vc-emit-ref-confidence GVCF --vc-sample-name ${PARENT_2_SAMPLE_NAME} --intermediate-results-dir ${TMP_DIR} --output-directory ${PARENT_2_DRAGEN_WORK_DIR} --output-file-prefix ${PARENT_2_SAMPLE_NAME}_dragen_genotypecall_wgs > ${COHORT_WORK_DIR}/${PARENT_2_SAMPLE_NAME}_workdir/dragen_call_${PARENT_2_SAMPLE_NAME}_wgs.gvcf.stdout 2> ${COHORT_WORK_DIR}/${PARENT_2_SAMPLE_NAME}_workdir/dragen_call_${PARENT_2_SAMPLE_NAME}_wgs.gvcf.stderr"
sh -t markellocj@helix.nih.gov ssh udpdragen01.nhgri.nih.gov "mkdir -p ${JOINT_GENOTYPE_DRAGEN_WORK_DIR}"

# RUN Dragen joint genotyping using parental GVCFs
JOINT_GENOTYPE_DRAGEN_WORK_DIR="/staging/markellocj/output_parental_joint_call_${COHORT_NAME}"
ssh -t markellocj@helix.nih.gov ssh udpdragen01.nhgri.nih.gov "mkdir -p ${JOINT_GENOTYPE_DRAGEN_WORK_DIR}/sample_gvcfs" && \
ssh -t markellocj@helix.nih.gov ssh udpdragen01.nhgri.nih.gov "cp ${PARENT_1_DRAGEN_WORK_DIR}/${PARENT_1_SAMPLE_NAME}_dragen_genotypecall_wgs.gvcf.gz ${PARENT_2_DRAGEN_WORK_DIR}/${PARENT_2_SAMPLE_NAME}_dragen_genotypecall_wgs.gvcf.gz ${JOINT_GENOTYPE_DRAGEN_WORK_DIR}/sample_gvcfs" && \
ssh -t markellocj@helix.nih.gov ssh udpdragen01.nhgri.nih.gov "touch ${JOINT_GENOTYPE_DRAGEN_WORK_DIR}/gvcf_list.txt" && \
ssh -t markellocj@helix.nih.gov ssh udpdragen01.nhgri.nih.gov "echo \"${JOINT_GENOTYPE_DRAGEN_WORK_DIR}/sample_gvcfs/${PARENT_1_SAMPLE_NAME}_dragen_genotypecall_wgs.gvcf.gz\" \>\> ${JOINT_GENOTYPE_DRAGEN_WORK_DIR}/gvcf_list.txt" && \
ssh -t markellocj@helix.nih.gov ssh udpdragen01.nhgri.nih.gov "echo \"${JOINT_GENOTYPE_DRAGEN_WORK_DIR}/sample_gvcfs/${PARENT_2_SAMPLE_NAME}_dragen_genotypecall_wgs.gvcf.gz\" \>\> ${JOINT_GENOTYPE_DRAGEN_WORK_DIR}/gvcf_list.txt" && \
ssh -t markellocj@helix.nih.gov ssh udpdragen01.nhgri.nih.gov "dragen -f -r /staging/hs37d5_v7 --enable-joint-genotyping true --intermediate-results-dir ${TMP_DIR} --output-directory ${JOINT_GENOTYPE_DRAGEN_WORK_DIR} --output-file-prefix parental_joint_genotyped_${COHORT_NAME} --variant-list ${JOINT_GENOTYPE_DRAGEN_WORK_DIR}/gvcf_list.txt > ${COHORT_WORK_DIR}/dragen_parental_joint_call_${COHORT_NAME}.stdout 2> ${COHORT_WORK_DIR}/dragen_parental_joint_call_${COHORT_NAME}.stderr"


# Copy final results to helix directory
mkdir /data/${UDPBINFO_DIR}/output_parental_joint_call_${COHORT_NAME} && chmod ug+rw -R /data/${UDPBINFO_DIR}/output_parental_joint_call_${COHORT_NAME} && \
ssh -t markellocj@helix.nih.gov ssh udpdragen01.nhgri.nih.gov "cp -R ${JOINT_GENOTYPE_DRAGEN_WORK_DIR}/. /staging/helix/${UDPBINFO_DIR}/output_parental_joint_call_${COHORT_NAME}" && \
ssh -t markellocj@helix.nih.gov ssh udpdragen01.nhgri.nih.gov "rm -fr ${JOINT_GENOTYPE_DRAGEN_WORK_DIR}/"

mv /data/${UDPBINFO_DIR}/output_parental_joint_call_${COHORT_NAME} ${COHORT_WORK_DIR}/output_parental_joint_call_${COHORT_NAME}
rm -fr /data/${UDPBINFO_DIR}/${PARENT_1_SAMPLE_NAME}_surjected_bams/ /data/${UDPBINFO_DIR}/${PARENT_2_SAMPLE_NAME}_surjected_bams/

