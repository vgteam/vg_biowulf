#!/bin/bash
#################################################################################################
##
##  Script to run an entire UDP cohort through the full VG workaround pipeline
##      (up to surjected BAMs) for Whole Genome data.
##
##  Inputs:
##
##  Assumptions:
##      - Reads are located in a different directory than in the working directory of where
##          the argument WORK_DIR is set.
##          - This is to get around the problem of binding two different directories when
##              running Singularity containers.
##
##  Last modified:
##  Last modified by: Charles Markello
##
#################################################################################################

## Create help statement
usage(){
cat << EOF

This script runs the VG pipeline up to GAM surjection to BAM files (no variant calling is made).

Inputs:
    -i Sample ID (in format UDP####)

Outputs:

Assumptions:

EOF

}

## Check whether script is being run on Biowulf
if [ "$HOSTNAME" == helix.nih.gov ]; then
    usage
    exit 1
fi

## Check number of arguments
if [ $# -lt 4 ] || [[ $@ != -* ]]; then
    usage
    exit 1
fi

PACKAGE_DIR="/data/markellocj/vg_workaround_pipeline_pkg_cleanup"

VG_CONTAINER="quay.io/vgteam/vg:v1.9.0-0-g3286e131-t206-run"
## Parse through arguments
while getopts "i:g:w:c:m:s:h" OPTION; do
    case $OPTION in
        i)
            COHORT_NAME=$OPTARG
        ;;
        g)
            GRAPH_FILES_DIR=$OPTARG
        ;;
        w)
            COHORT_SET_WORK_DIR=$OPTARG
        ;;
        c)
            VG_CONTAINER=$OPTARG
        ;;
        m)
            MAP_ALGORITHM=$OPTARG
        ;;
        s)
            READS_PER_CHUNK=$OPTARG
        ;;
        h)
            usage
            exit 1
        ;;
        \?)
            usage
            exit 1
        ;;
    esac
done

## STEP1: COLLECT READS. Generate sample work directories and create read gathering swarm script
COHORT_WORK_DIR="${COHORT_SET_WORK_DIR}/${COHORT_NAME}_run"
COHORT_DATA_DIR="/data/Udpdata/Families/${COHORT_NAME}"
COHORT_NAMES_LIST=($(ls $COHORT_DATA_DIR/SnpChip/MostRecent | grep 'UDP' | awk -F'_' '{print $1}' | awk -F'.' '{print $1}' | uniq))
CAT_READS_SWARMFILE_PATH="${COHORT_WORK_DIR}/cat_reads_swarmfile_${COHORT_NAME}"

mkdir -p ${COHORT_WORK_DIR}

echo "Cohort UDP Sample List: ${COHORT_NAMES_LIST[@]}"

for SAMPLE_NAME in ${COHORT_NAMES_LIST[@]}
do
    ## Make sample work directory and sample readfile directories
    SAMPLE_WORK_DIR="${COHORT_WORK_DIR}/${SAMPLE_NAME}_workdir"
    mkdir -p ${SAMPLE_WORK_DIR}
    mkdir ${SAMPLE_WORK_DIR}/${SAMPLE_NAME}_reads
    
    ## Concatenate reads and put them in the sample-specific readfile directory
    INDIVIDUAL_DATA_DIR="/data/Udpdata/Individuals/${SAMPLE_NAME}"
    PAIR_1_READS=()
    PAIR_2_READS=()
    LANE_NUMS=($(ls ${INDIVIDUAL_DATA_DIR}/WGS/Rawreads/HudsonAlpha_1 | awk -F'-' '{print $2}'| awk -F'_' '{print $1"_"$2}' | sort | uniq | xargs))
    for LANE_NUM in ${LANE_NUMS[@]}
    do
        PAIR_1_READS+=(${INDIVIDUAL_DATA_DIR}/WGS/Rawreads/HudsonAlpha_1/"$(ls ${INDIVIDUAL_DATA_DIR}/WGS/Rawreads/HudsonAlpha_1 | grep "${LANE_NUM}_1")")
        PAIR_2_READS+=(${INDIVIDUAL_DATA_DIR}/WGS/Rawreads/HudsonAlpha_1/"$(ls ${INDIVIDUAL_DATA_DIR}/WGS/Rawreads/HudsonAlpha_1 | grep "${LANE_NUM}_2")")
    done

    echo "cat ${PAIR_1_READS[@]} > ${SAMPLE_WORK_DIR}/${SAMPLE_NAME}_reads/${SAMPLE_NAME}_read_pair_1.fq.gz" >> ${CAT_READS_SWARMFILE_PATH}
    echo "cat ${PAIR_2_READS[@]} > ${SAMPLE_WORK_DIR}/${SAMPLE_NAME}_reads/${SAMPLE_NAME}_read_pair_2.fq.gz" >> ${CAT_READS_SWARMFILE_PATH}
done

COLLECT_READS_JOBID=$(swarm -f ${CAT_READS_SWARMFILE_PATH} --time=4:00:00)
echo "Running read collection swarm script. Jobid:${COLLECT_READS_JOBID}"

cd ${COHORT_WORK_DIR}

## STEP2: RUN MAPPING PIPELINE. Run the graph mapping pipeline per sample.
for SAMPLE_NAME in ${COHORT_NAMES_LIST[@]}
do
    SAMPLE_WORK_DIR="${COHORT_WORK_DIR}/${SAMPLE_NAME}_workdir"
    READ_FILE_1="${SAMPLE_WORK_DIR}/${SAMPLE_NAME}_reads/${SAMPLE_NAME}_read_pair_1.fq.gz"
    READ_FILE_2="${SAMPLE_WORK_DIR}/${SAMPLE_NAME}_reads/${SAMPLE_NAME}_read_pair_2.fq.gz"
    SAMPLE_MAP_RUN_JOBID=$(sbatch --cpus-per-task=32 --mem=100g --time=4:00:00 --dependency=afterok:${COLLECT_READS_JOBID} ${PACKAGE_DIR}/sample_master_script.sh -i ${SAMPLE_NAME} -r ${READ_FILE_1} -r ${READ_FILE_2} -g ${GRAPH_FILES_DIR} -w ${SAMPLE_WORK_DIR} -c ${VG_CONTAINER} -m ${MAP_ALGORITHM} -s ${READS_PER_CHUNK} > ${SAMPLE_NAME}_master_script.stdout)
    echo "Running sample ${SAMPLE_NAME} through the graph alignment pipeline. Jobid:${SAMPLE_MAP_RUN_JOBID}"
done


exit

