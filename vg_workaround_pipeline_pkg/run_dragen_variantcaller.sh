#!/bin/bash

SAMPLE_NAME=$1
WORK_DIR=$2
UDPBINFO_DIR=$3

# Example INPUT:
#SAMPLE_NAME="NA24385"
#WORK_DIR="/data/markellocj/HG002_run/output_vg_t206/output_NA24385_vg_t206_poscontrol_ref_mpmap"
#UDPBINFO_DIR="Udpbinfo/usr/markellocj"

cp -r ${WORK_DIR}/${SAMPLE_NAME}_surjected_bams/ /data/${UDPBINFO_DIR}

DRAGEN_WORK_DIR="/staging/markellocj/${SAMPLE_NAME}"
TMP_DIR="/staging/markellocj/tmp"

ssh -t markellocj@helix.nih.gov ssh udpdragen01.nhgri.nih.gov "mkdir -p ${DRAGEN_WORK_DIR}" && \
ssh -t markellocj@helix.nih.gov ssh udpdragen01.nhgri.nih.gov "mkdir -p ${TMP_DIR}" && \
ssh -t markellocj@helix.nih.gov ssh udpdragen01.nhgri.nih.gov "dragen -f -r /staging/hs37d5_v7 -b /staging/helix/${UDPBINFO_DIR}/${SAMPLE_NAME}_surjected_bams/${SAMPLE_NAME}_merged.sorted.dupmarked.reordered.bam --verbose --bin_memory=50000000000 --enable-map-align false --enable-variant-caller true --pair-by-name=true --vc-sample-name ${SAMPLE_NAME} --intermediate-results-dir /staging/markellocj/tmp --output-directory ${DRAGEN_WORK_DIR} --output-file-prefix ${SAMPLE_NAME}_dragen_genotypecall_wgs > ${WORK_DIR}/dragen_call_${SAMPLE_NAME}_wgs.stdout 2> ${WORK_DIR}/dragen_call_${SAMPLE_NAME}_wgs.stderr"

mkdir /data/${UDPBINFO_DIR}/${SAMPLE_NAME}_dragen_genotyper && chmod ug+rw -R /data/${UDPBINFO_DIR}/${SAMPLE_NAME}_dragen_genotyper && \
ssh -t markellocj@helix.nih.gov ssh udpdragen01.nhgri.nih.gov "cp -R ${DRAGEN_WORK_DIR}/. /staging/helix/${UDPBINFO_DIR}/${SAMPLE_NAME}_dragen_genotyper" && \
ssh -t markellocj@helix.nih.gov ssh udpdragen01.nhgri.nih.gov "rm -fr ${DRAGEN_WORK_DIR}/"

mv /data/${UDPBINFO_DIR}/${SAMPLE_NAME}_dragen_genotyper ${WORK_DIR}/${SAMPLE_NAME}_dragen_genotyper
rm -fr /data/${UDPBINFO_DIR}/${SAMPLE_NAME}_surjected_bams/

