#!/bin/bash

SAMPLE_NAME=$1
WORK_DIR=$2
PROCESS_BAM_CORES=32

cd ${WORK_DIR} 

# Merge the chromosomal BAM files into a single sample wgs BAM file
module load samtools
SWARM_MAP_OUTPUT_DIR_NAME="${SAMPLE_NAME}_surjected_bams"
SWARM_MAP_OUTPUT_DIR_PATH="${WORK_DIR}/${SWARM_MAP_OUTPUT_DIR_NAME}"
time samtools merge -f --threads ${PROCESS_BAM_CORES} ${SWARM_MAP_OUTPUT_DIR_PATH}/${SAMPLE_NAME}_merged.dupmarked.reordered.bam ${SWARM_MAP_OUTPUT_DIR_NAME}/*.sorted.dupmarked.reordered.bam
time samtools sort -T ${SWARM_MAP_OUTPUT_DIR_NAME}/tmpSort --threads 32 ${SWARM_MAP_OUTPUT_DIR_NAME}/${SAMPLE_NAME}_merged.dupmarked.reordered.bam > ${SWARM_MAP_OUTPUT_DIR_NAME}/${SAMPLE_NAME}_merged.sorted.dupmarked.reordered.bam
# Index the final merged and reordered BAM file
time samtools index -@ ${PROCESS_BAM_CORES} ${SWARM_MAP_OUTPUT_DIR_PATH}/${SAMPLE_NAME}_merged.sorted.dupmarked.reordered.bam


