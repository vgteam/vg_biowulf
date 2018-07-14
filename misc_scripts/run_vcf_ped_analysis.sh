#!/bin/bash

# Example run: ./run_vcf_ped_analysis.sh UDP10618 t146

COHORT_NAME=$1
VG_VERSION=$2
COHORT_WORKDIR="${HOME}/${COHORT_NAME}_run"
JOINT_VCF="${COHORT_WORKDIR}/output_vg_${VG_VERSION}/${COHORT_NAME}_dragen_call_from_vg_surject_vg_${VG_VERSION}/joint_genotyped_${COHORT_NAME}.vcf.gz"
OUTPUT="${COHORT_WORKDIR}/output_vg_${VG_VERSION}/mendelian_analysis_vg_${VG_VERSION}/${COHORT_NAME}_mendelian_inconsistent.vcf.gz"
PED_FILE="${HOME}/${COHORT_NAME}_run/${COHORT_NAME}.ped"
WORK_DIR="/data/markellocj/"
SDF_DIR="${HOME}/sdf_references/hs37d5.fa.sdf"

module load singularity

## Run pedigree analysis on vg surjected and Dragen called variants
cd ${WORK_DIR}
mkdir -p ${COHORT_NAME}_run/output_vg_${VG_VERSION}/mendelian_analysis_vg_${VG_VERSION}
singularity exec -H ${PWD}:${HOME}  docker://realtimegenomics/rtg-tools rtg mendelian -i ${JOINT_VCF} --output-inconsistent ${OUTPUT} --pedigree ${PED_FILE} -l -t ${SDF_DIR}


