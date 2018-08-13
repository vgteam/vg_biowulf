#!/bin/bash

# Example run: ./run_bamqc.sh UDP10927 /data/markellocj/UDP_PEDIGREE_EXPERIMENT_1/UDP10927_run/UDP10927_workdir/UDP10927_surjected_bams

SAMPLE_NAME=$1
WORK_DIR=$2

## Run bamqc on merged surjected BAM file
module load singularity
cd ${WORK_DIR}
mkdir -p bamqc_out
singularity exec -H ${PWD}:${HOME} --pwd ${HOME} docker://quay.io/cmarkello/bamqc bamqc --outdir=/home/markellocj/bamqc_out/ ${SAMPLE_NAME}_merged.sorted.dupmarked.reordered.bam

## Run mosdepth on surjected BAM file
module load mosdepth
cd ${WORK_DIR}
mkdir -p mosdepth_output
mosdepth -t 16 "${WORK_DIR}/mosdepth_output/${SAMPLE_NAME}_VG" ${SAMPLE_NAME}_merged.sorted.dupmarked.reordered.bam

