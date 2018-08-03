#!/bin/bash

## Assumes xg, gcsa, id_ranges, and cid_ranges files names

SAMPLE_NAME=$1
WORK_DIR=$2
GRAPH_FILES_DIR_PATH=$3
VG_CONTAINER=$4
PACKAGE_DIR=$5

FASTQ_DIR_NAME="split_fastqs_${SAMPLE_NAME}"
FASTQ_CHUNK_WORKDIR_PATH="${WORK_DIR}/${FASTQ_DIR_NAME}"
SWARM_MAP_OUTPUT_DIR_NAME="${SAMPLE_NAME}_output_gams"
BED_DIR_NAME="gam_refactor_bed_${SAMPLE_NAME}"
BED_WORKDIR_PATH="${WORK_DIR}/${BED_DIR_NAME}"
BED_CONTAINER_PATH="${HOME}/${BED_DIR_NAME}"
SWARM_MAP_OUTPUT_WORKDIR_PATH="${WORK_DIR}/${SWARM_MAP_OUTPUT_DIR_NAME}"
SWARM_MAP_OUTPUT_CONTAINER_PATH="${HOME}/${SWARM_MAP_OUTPUT_DIR_NAME}"
MAP_SWARMFILE_NAME="${WORK_DIR}/map_swarmfile_${SAMPLE_NAME}"
ALIGNMENT_CORES=32

GAM_DIR_NAME="${SAMPLE_NAME}_output_gams"
CHR_OUTPUT_BAM_DIR_NAME="${SAMPLE_NAME}_surjected_bams"
CHR_OUTPUT_BAM_DIR_PATH="${WORK_DIR}/${CHR_OUTPUT_BAM_DIR_NAME}"
SURJECT_GAMS_SWARMFILE_NAME="${WORK_DIR}/surject_gams_swarmfile_${SAMPLE_NAME}"
SURJECT_CORES=32

SWARM_MAP_OUTPUT_DIR_NAME="${SAMPLE_NAME}_surjected_bams"
SWARM_MAP_OUTPUT_DIR_PATH="${WORK_DIR}/${SWARM_MAP_OUTPUT_DIR_NAME}"
PROCESS_BAMS_SWARMFILE_NAME="${WORK_DIR}/process_bams_swarmfile_${SAMPLE_NAME}"
PROCESS_BAM_CORES=32

cd ${WORK_DIR}

rm ${MAP_SWARMFILE_NAME} ${SURJECT_GAMS_SWARMFILE_NAME} ${PROCESS_BAMS_SWARMFILE_NAME}

if [ ! -d "${SWARM_MAP_OUTPUT_WORKDIR_PATH}" ]; then
    mkdir -p ${SWARM_MAP_OUTPUT_WORKDIR_PATH}
    chmod 2770 ${SWARM_MAP_OUTPUT_WORKDIR_PATH}
fi

if [ ! -d "${CHR_OUTPUT_BAM_DIR_PATH}" ]; then
    mkdir -p ${CHR_OUTPUT_BAM_DIR_PATH}
    chmod 2770 ${CHR_OUTPUT_BAM_DIR_PATH}
fi

# Extract chunk ids (e.g. fq_chunk_1.part.aa.fq.gz and fq_chunk_2.part.aa.fq.gz --> aa)
read_chunk_ids=($(ls -l ${FASTQ_CHUNK_WORKDIR_PATH} | awk -F'.' '{print $3}' | sort | uniq | xargs))
XG_FILE_NAME=($(ls ${GRAPH_FILES_DIR_PATH} | grep 'xg$'))
GCSA_FILE_NAME=($(ls ${GRAPH_FILES_DIR_PATH} | grep 'gcsa$'))

# Iterate for each chunk id and run vg map for each chunked fastq
for i in ${read_chunk_ids[@]}
do
    echo $i
    READ_FILE_1_CONTAINER_PATH=${HOME}/${FASTQ_DIR_NAME}/fq_chunk_1.part.${i}.fq.gz
    READ_FILE_2_CONTAINER_PATH=${HOME}/${FASTQ_DIR_NAME}/fq_chunk_2.part.${i}.fq.gz
    map_command="module load singularity; export TMPDIR=/lscratch/${SLURM_JOB_ID}; cd ${WORK_DIR}; time singularity -q exec -H ${WORK_DIR}:${HOME} --pwd ${HOME} -B ${GRAPH_FILES_DIR_PATH}:/mnt docker://${VG_CONTAINER} vg mpmap -S -f ${READ_FILE_1_CONTAINER_PATH} -f ${READ_FILE_2_CONTAINER_PATH} -x /mnt/${XG_FILE_NAME} -g /mnt/${GCSA_FILE_NAME} -t ${ALIGNMENT_CORES} > ${SWARM_MAP_OUTPUT_WORKDIR_PATH}/${SAMPLE_NAME}_${i}.gam"
    echo $map_command >> ${MAP_SWARMFILE_NAME}

    GAM_FILE="${HOME}/${GAM_DIR_NAME}/${SAMPLE_NAME}_${i}.gam"
    surject_command="module load singularity; export TMPDIR=/lscratch/${SLURM_JOB_ID}; cd ${WORK_DIR}; time singularity -q exec -H ${WORK_DIR}:${HOME} --pwd ${HOME} -B ${GRAPH_FILES_DIR_PATH}:/mnt docker://${VG_CONTAINER} vg surject -i -x /mnt/${XG_FILE_NAME} -t ${SURJECT_CORES} -b ${GAM_FILE} > ${CHR_OUTPUT_BAM_DIR_PATH}/${SAMPLE_NAME}_${i}.bam"
    echo $surject_command >> ${SURJECT_GAMS_SWARMFILE_NAME}

    process_bam_command="module load samtools picard && cd ${WORK_DIR} && time samtools sort -T ${SWARM_MAP_OUTPUT_DIR_PATH}/tmpSort_${i} --threads ${PROCESS_BAM_CORES} ${SWARM_MAP_OUTPUT_DIR_PATH}/${SAMPLE_NAME}_${i}.bam > ${SWARM_MAP_OUTPUT_DIR_PATH}/${SAMPLE_NAME}_${i}.sorted.bam && time java -Xmx4g -XX:ParallelGCThreads=5 -jar \$PICARDJARPATH/picard.jar MarkDuplicates VALIDATION_STRINGENCY=LENIENT I=${SWARM_MAP_OUTPUT_DIR_PATH}/${SAMPLE_NAME}_${i}.sorted.bam O=${SWARM_MAP_OUTPUT_DIR_PATH}/${SAMPLE_NAME}_${i}.sorted.dupmarked.bam M=${SWARM_MAP_OUTPUT_DIR_PATH}/marked_dup_metrics_chr${i}.txt && rm -f ${SWARM_MAP_OUTPUT_DIR_PATH}/${SAMPLE_NAME}_${i}.sorted.bam && time java -Xmx20g -XX:ParallelGCThreads=32 -jar \$PICARDJARPATH/picard.jar ReorderSam VALIDATION_STRINGENCY=LENIENT INPUT=${SWARM_MAP_OUTPUT_DIR_PATH}/${SAMPLE_NAME}_${i}.sorted.dupmarked.bam OUTPUT=${SWARM_MAP_OUTPUT_DIR_PATH}/${SAMPLE_NAME}_${i}.sorted.dupmarked.reordered.bam REFERENCE=/data/markellocj/fasta_references/hs37d5_reference/hs37d5.fa && rm -f ${SWARM_MAP_OUTPUT_DIR_PATH}/${SAMPLE_NAME}_${i}.sorted.dupmarked.bam"
    echo $process_bam_command >> ${PROCESS_BAMS_SWARMFILE_NAME}    

done


