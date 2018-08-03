#!/bin/bash
#################################################################################################
##
##  Script to run the full VG workaround pipeline (up to surjected BAMs) for Whole Genome data
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
if [ $# -lt 5 ] || [[ $@ != -* ]]; then
    usage
    exit 1
fi

PACKAGE_DIR="/data/markellocj/vg_workaround_pipeline_pkg_cleanup"

VG_CONTAINER="quay.io/vgteam/vg:v1.9.0-0-g3286e131-t206-run"
READS_PER_CHUNK=10000000

## Parse through arguments
while getopts "i:r:g:w:c:m:s:h" OPTION; do
    case $OPTION in
        i)
            SAMPLE_NAME=$OPTARG
        ;;
        r)
            FASTQ_FILES+=($OPTARG)
        ;;
        g)
            GRAPH_FILES_DIR=$OPTARG
        ;;
        w)
            WORK_DIR=$OPTARG
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

## STEP1: SPLIT READS. Generate and run read splitting swarm jobscript.
INPUT_READ_FILE_1=${FASTQ_FILES[0]}
INPUT_READ_FILE_2=${FASTQ_FILES[1]}

echo "SAMPLE_NAME: ${SAMPLE_NAME}"
echo "INPUT_READ_FILE_1: ${INPUT_READ_FILE_1}"
echo "INPUT_READ_FILE_2: ${INPUT_READ_FILE_2}"
echo "GRAPH_FILES_DIR: ${GRAPH_FILES_DIR}"
echo "WORK_DIR: ${WORK_DIR}"
echo "VG_CONTAINER: ${VG_CONTAINER}"

if [ ! -d "${WORK_DIR}" ]; then
    mkdir -p ${WORK_DIR}
    chmod 2770 ${WORK_DIR}
fi

echo "Splitting reads files into chunks of ${READS_PER_CHUNK} reads."
${PACKAGE_DIR}/generate_split_reads_swarm_script.sh ${INPUT_READ_FILE_1} ${INPUT_READ_FILE_2} ${SAMPLE_NAME} ${WORK_DIR} ${VG_CONTAINER} ${READS_PER_CHUNK}
echo "Finished splitting reads."


## RUN EITHER VG MAP OR VG MPMAP ALGORITHM PIPELINES.
if [ "$MAP_ALGORITHM" = "vg_mpmap" ]; then
    ## STEP2: GENERATE SURJECT SWARM FILES. Generate swarm script to run vg surject on each chunked GAM into final BAM files.
    ${PACKAGE_DIR}/generate_vg_mpmap_pipeline_swarm_scripts.sh ${SAMPLE_NAME} ${WORK_DIR} ${GRAPH_FILES_DIR} ${VG_CONTAINER} ${PACKAGE_DIR}
    MAP_SWARMFILE_NAME="${WORK_DIR}/map_swarmfile_${SAMPLE_NAME}"
    SURJECT_GAMS_SWARMFILE_NAME="${WORK_DIR}/surject_gams_swarmfile_${SAMPLE_NAME}"    
    PROCESS_BAMS_SWARMFILE_NAME="${WORK_DIR}/process_bams_swarmfile_${SAMPLE_NAME}"   
    echo "Generated swarm scripts for vg mpmap alignment pipeline."
 
    ## STEP3: CHUNK ALIGNMENT and RUN SURJECTION.
    CHUNK_ALIGNMENT_JOBID=$(swarm -f ${MAP_SWARMFILE_NAME} -g 120 -t 32 --time 6:00:00 --gres=lscratch:200 --maxrunning 10)
    echo "Running swarm VG mpmap graph alignment. Jobid:${CHUNK_ALIGNMENT_JOBID}"

    SURJECT_GAMS_JOBID=$(swarm -f ${SURJECT_GAMS_SWARMFILE_NAME} --dependency=afterok:${CHUNK_ALIGNMENT_JOBID} -g 100 -t 32 --time 12:00:00 --gres=lscratch:200)
    echo "Running surject swarm script. Jobid:${SURJECT_GAMS_JOBID}"

else
    ## STEP2: GENERATE SURJECT SWARM FILES. Generate swarm script to run vg map on each fastq chunk and process the final BAM files.
    ${PACKAGE_DIR}/generate_vg_map_surject_swarm_scripts.sh ${SAMPLE_NAME} ${WORK_DIR} ${GRAPH_FILES_DIR} ${VG_CONTAINER} ${PACKAGE_DIR}
    MAP_SWARMFILE_NAME="${WORK_DIR}/map_swarmfile_${SAMPLE_NAME}"
    PROCESS_BAMS_SWARMFILE_NAME="${WORK_DIR}/process_bams_swarmfile_${SAMPLE_NAME}"
    echo "Generated swarm scripts for vg map alignment pipeline."
    
    ## STEP3: CHUNK ALIGNMENT WITH VG SURJECTION.
    SURJECT_GAMS_JOBID=$(swarm -f ${MAP_SWARMFILE_NAME} -g 100 -t 32 --time 8:00:00 --gres=lscratch:200 --maxrunning 15)
    echo "Running swarm VG map graph alignment. Jobid:${SURJECT_GAMS_JOBID}"

fi
    
## STEP4: SORT, MARKDUPLICATES and REORDER BAMs.
PROCESS_BAMS_JOBID=$(swarm -f ${PROCESS_BAMS_SWARMFILE_NAME} --dependency=afterok:${SURJECT_GAMS_JOBID} -g 100 -t 32 --time 12:00:00 --maxrunning 20)
echo "Running process BAMs swarm script. Jobid:${PROCESS_BAMS_JOBID}"

## STEP5: MERGE BAM FILES.
MERGE_BAMS_JOBID=$(sbatch --cpus-per-task=32 --mem=100g --time=12:00:00 --dependency=afterok:${PROCESS_BAMS_JOBID} ${PACKAGE_DIR}/run_merge_processed_bam_script.sh ${SAMPLE_NAME} ${WORK_DIR})
echo "Running merge of final BAMs. Jobid:${MERGE_BAMS_JOBID}"

## STEP6: CLEAN UP WORK DIRECTORY.
CLEAN_WORKDIR_JOBID=$(sbatch --time=1:00:00 --dependency=afterok:${MERGE_BAMS_JOBID} ${PACKAGE_DIR}/run_workdir_cleanup.sh ${SAMPLE_NAME} ${WORK_DIR} ${PACKAGE_DIR})
echo "Running work directory cleanup script. Jobid:${CLEAN_WORKDIR_JOBID}"


exit

