#!/usr/bin/env bash

# Author: Julian Q. Zhou
# https://github.com/julianqz
# Date:   2024-02-28
#
# Run `cellranger vdj` for a list of BCR or TCR samples
#
# Prereqs:  
# 1) The following must be in ${PROJ_ID}/aux/
#    - a sample list: "cr_list_${RECEPTOR}_${PROJ_ID}.txt"
#      each row is semi-colon-separated
#
#      EITHER (if using a centralized fastq dir for all samples)
#             [sample];[comma-separated fastq id(s)]
# 
#      OR     (if using a different fastq dir for each sample)
#             [sample];[comma-separated fastq id(s)];[(a single) fastq dir]


# Print usage
usage () {
    echo -e "Usage: `basename $0` [OPTIONS]"
    echo -e "  -J  Project ID."
    echo -e "  -K  Receptor type. Either 'bcr' or 'tcr'."                        
    echo -e "  -T  Path to the top-level working dir." 
    echo -e "  -L  Whether to use a centralized fastq dir. Default is true." 
    echo -e "  -F  Path to the centralized fastq dir." 
    echo -e "  -R  Path to the reference dir." 
    echo -e "  -Y  Number of cores for cellranger."    
    echo -e "  -Z  Amount of memory for cellranger."
    echo -e "  -W  Delete .bam* files. Default is true."
    echo -e "  -Q  Whether to specify --chain. Default is false."
    echo -e "  -h  This message."
}

PROJ_ID_SET=false
RECEPTOR_SET=false
PATH_ROOT_SET=false
BOOL_ONE_FASTQ_DIR_SET=false
PATH_FASTQ_SET=false
PATH_REF_SET=false
BOOL_SPECIFY_CHAIN_SET=false
BOOL_DEL_BAM_SET=false

# Get commandline arguments
while getopts "J:K:T:L:F:R:Y:Z:W:Q:h" OPT; do
    case "$OPT" in
    J)  PROJ_ID=$OPTARG
        PROJ_ID_SET=true
        ;;
    K)  RECEPTOR=$OPTARG
        RECEPTOR_SET=true
        ;;
    T)  PATH_ROOT=$(realpath $OPTARG)
        PATH_ROOT_SET=true
        ;;
    L)  BOOL_ONE_FASTQ_DIR=$OPTARG
        BOOL_ONE_FASTQ_DIR_SET=true
        ;;
    F)  PATH_FASTQ=$(realpath $OPTARG)
        PATH_FASTQ_SET=true
        ;;
    R)  PATH_REF=$(realpath $OPTARG)
        PATH_REF_SET=true
        ;;
    Y)  CR_N=$OPTARG
        ;;
    Z)  CR_M=$OPTARG
        ;;
    W)  BOOL_DEL_BAM=$OPTARG
        BOOL_DEL_BAM_SET=true
        ;;
    Q)  BOOL_SPECIFY_CHAIN=$OPTARG
        BOOL_SPECIFY_CHAIN_SET=true
        ;;
    h)  usage
        exit
        ;;
    \?) echo "Invalid option: -$OPTARG" >&2
        exit 1
        ;;
    :)  echo "Option -$OPTARG requires an argument" >&2
        exit 1
        ;;
    esac
done

# Exit if no project ID provided
if ! $PROJ_ID_SET; then
    echo "You must specify a project ID via the -J option" >&2
    exit 1
fi

# Exit if no receptor type provided
if ! $RECEPTOR_SET; then
    echo "You must specify a receptor type via the -K option" >&2
    exit 1
fi

# Exit if no top-level working directory provided
if ! $PATH_ROOT_SET; then
    echo "You must specify a top-level working directory that contains the project folder via the -T option" >&2
    exit 1
fi


# Set BOOL_ONE_FASTQ_DIR to true if no -L specified
if ! $BOOL_ONE_FASTQ_DIR_SET; then
    BOOL_ONE_FASTQ_DIR=true
fi

# If BOOL_ONE_FASTQ_DIR true, exit if no centralized fastq directory provided
if [[ $BOOL_ONE_FASTQ_DIR && ! $PATH_FASTQ_SET ]]; then
    echo "You must specify a centralized fastq directory that contains all the fastq files via the -F option" >&2
    exit 1
fi

# Exit if no reference directory provided
if ! $PATH_REF_SET; then
    echo "You must specify a reference directory that contains the genome references via the -R option" >&2
    exit 1
fi


# Set BOOL_DEL_BAM to true if no -W specified
if ! $BOOL_DEL_BAM_SET; then
    BOOL_DEL_BAM=true
fi


# Set BOOL_SPECIFY_CHAIN to false if no -Q specified
if ! $BOOL_SPECIFY_CHAIN_SET; then
    BOOL_SPECIFY_CHAIN=false
fi

# --chain

# If -Q/BOOL_SPECIFY_CHAIN true, manually specify --chain based on $RECEPTOR
#    bcr >>> IG
#    tcr >>> TR
# Else   >>> auto

if $BOOL_SPECIFY_CHAIN; then
    if [[ ${RECEPTOR} == "bcr" ]]; then
        CHAIN="IG"
    elif [[ ${RECEPTOR} == "tcr" ]]; then
        CHAIN="TR"
    else
        echo "-K/RECEPTOR should be either bcr or tcr" >&2
        exit 1
    fi
else
    CHAIN="auto"
fi

# paths

PATH_PROJ="${PATH_ROOT}/${PROJ_ID}"
# no error if existing
mkdir -p "${PATH_PROJ}"

# csv, txt, logs
PATH_AUX="${PATH_PROJ}/aux/"
mkdir -p "${PATH_AUX}"

# cellranger vdj outputs
PATH_OUTPUT="${PATH_PROJ}/cr_${RECEPTOR}/"
mkdir -p "${PATH_OUTPUT}"

PATH_LOG="${PATH_AUX}log_cr_${RECEPTOR}_$(date '+%m%d%Y_%H%M%S').log"

NAME_LIST="cr_list_${RECEPTOR}_${PROJ_ID}.txt"
PATH_LIST="${PATH_AUX}${NAME_LIST}" 


cellranger --version &> "${PATH_LOG}"
echo "--localcores=${CR_N}; --localmem=${CR_M}" &>> "${PATH_LOG}"
echo "--chain=${CHAIN}" &>> "${PATH_LOG}"
echo "Whethe to use a single centralized FASTQ folder: ${BOOL_ONE_FASTQ_DIR}" &>> "${PATH_LOG}"
echo "Sample list: ${NAME_LIST}" &>> "${PATH_LOG}"

N_LINES=$(wc -l < "${PATH_LIST}")
echo "N_LINES: ${N_LINES}" &>> "${PATH_LOG}"


cd "${PATH_OUTPUT}"

for ((IDX=1; IDX<=${N_LINES}; IDX++)); do

    # read current line
    CUR_LINE=$(sed "${IDX}q;d" "${PATH_LIST}") 

    # important for .txt to be semi-colon-separated
    # this allows for multiple fastq ids to be comma-separated for the same sample, if applicable
    IFS=";"
    read -a strarr <<< "${CUR_LINE}"

	# sample ID
	CUR_ID=${strarr[0]}

    # fastq id(s) (this is what cellranger calls `sample` -- confusing eh?)
    CUR_FASTQ_IDS=${strarr[1]}

    # fastq dir
    if $BOOL_ONE_FASTQ_DIR; then
        # using a single centralized fastq dir
        CUR_PATH_FASTQ="${PATH_FASTQ}"
    else
        # using individualized fastq dir
        CUR_PATH_FASTQ=${strarr[2]}
    fi


	echo "IDX: ${IDX}; CUR_ID: ${CUR_ID}" &>> "${PATH_LOG}"
    echo "   FASTQ dir: ${CUR_PATH_FASTQ}" &>> "${PATH_LOG}"

	# sample-specific log
	PATH_LOG_ID="${PATH_AUX}log_cr_${RECEPTOR}_${IDX}_${CUR_ID}_$(date '+%m%d%Y_%H%M%S').log"

    # vdj does not have --nosecondary option
	cellranger vdj \
		--id "${CUR_ID}" \
		--fastqs "${CUR_PATH_FASTQ}" \
        --sample "${CUR_FASTQ_IDS}" \
        --reference "${PATH_REF}" \
		--localcores "${CR_N}" \
		--localmem "${CR_M}" \
        --chain "${CHAIN}" \
		&> "${PATH_LOG_ID}"

    # remove .bam, .bambi, etc.

    if $BOOL_DEL_BAM; then
        rm "${PATH_OUTPUT}${CUR_ID}"/outs/*.bam*
    fi

    rm "${PATH_OUTPUT}${CUR_ID}"/_*
    rm -r "${PATH_OUTPUT}${CUR_ID}/SC_VDJ_ASSEMBLER_CS"
    rm "${PATH_OUTPUT}${CUR_ID}/${CUR_ID}.mri.tgz"

done

echo "ALL DONE" &>> "${PATH_LOG}"
