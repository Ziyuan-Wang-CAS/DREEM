#!/bin/bash
# ===================================================================

REF_NAME=human_mitochondrial_genome-NCBI_NC_012920.1  # 参考序列的名称
POS_START=12337                                        # 分析的起始位置
POS_END=14148                                          # 分析的终止位置
DATA_NAME=MTS_s4U_APDF20_0906_rep1                            # 分析样本的名称
FASTQ=no                                             # 是否从fastq文件开始进行分析
CTRL=no                                              # 是否为对照组样本
STRUCT=yes                                             # 是否使用RNAStructure进行结构预测
# -------------------------------------------------------------------
CHUNK_SIZE=100                                        # 每个任务的分析区间大小
# -------------------------------------------------------------------
TOTAL_LENGTH=$((POS_END - POS_START + 1))
NUM_TASKS=$(( (TOTAL_LENGTH + CHUNK_SIZE - 1) / CHUNK_SIZE ))
# -------------------------------------------------------------------
CUR_DATE=$(date +%Y%m%d)
QSUB_LOG_DIR=/sibcb1/hanshuolab1/wangziyuan/qsub_logs/DREEM/${CUR_DATE}_${DATA_NAME}_${POS_START}_${POS_END}/
mkdir -p "$QSUB_LOG_DIR"
# -------------------------------------------------------------------
qsub -q b1.q \
    -l hostname=fnode006 \
    -wd $QSUB_LOG_DIR \
    -t 1:$NUM_TASKS \
    -v REF_NAME="$REF_NAME",POS_START="$POS_START",POS_END="$POS_END",DATA_NAME="$DATA_NAME",FASTQ="$FASTQ",CTRL="$CTRL",STRUCT="$STRUCT",CHUNK_SIZE="$CHUNK_SIZE" \
    DREEM_Parallel_Run.sh