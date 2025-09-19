#!/bin/bash
# ===================================================================

eval "$(conda shell.bash hook)"
conda activate dreem || {
    echo "[error] conda environment 'dreem' not found!";
    exit 1;
}
# -------------------------------------------------------------------
REF_NAME=human_mitochondrial_genome-NCBI_NC_012920.1  # 参考序列的名称
POS_START=12337                                        # 分析的起始位置
POS_END=14148                                          # 分析的终止位置
DATA_NAME=MTS_s4U_APDF20_0906_rep1                            # 分析样本的名称
FASTQ=no                                             # 是否从fastq文件开始进行分析
CTRL=no                                              # 是否为对照组样本
STRUCT=yes                                             # 是否使用RNAStructure进行结构预测
# -------------------------------------------------------------------
WORK_DIR=/sibcb1/hanshuolab1/wangziyuan/DREEM/code/
CUR_DATE=$(date +%Y%m%d)
REF_FILE=${REF_NAME}.fasta
RESULTS_DIR=../results_${CUR_DATE}_${DATA_NAME}_${POS_START}_${POS_END}/
LOG_DIR=../logs/${CUR_DATE}_${DATA_NAME}_${POS_START}_${POS_END}/
# -------------------------------------------------------------------
if [ "$FASTQ" = "yes" ]; then
    FASTQ_FLAG="--fastq"
else
    FASTQ_FLAG=""
fi
if [ "$CTRL" = "yes" ]; then
    CTRL_FLAG="--ctrl"
else
    CTRL_FLAG=""
fi
if [ "$STRUCT" = "yes" ]; then
    STRUCT_FLAG="--struct"
else
    STRUCT_FLAG=""
fi
# -------------------------------------------------------------------
cd $WORK_DIR
mkdir -p "$LOG_DIR"
cd /sibcb1/hanshuolab1/wangziyuan/DREEM/data/$DATA_NAME
if [ "$FASTQ" = "yes" ]; then
    if [ -f "$REF_FILE" ]; then
        bowtie2-build "$REF_FILE" "$REF_NAME"
        echo "[INFO] Bowtie2 index builded for $REF_FILE" >> "${WORK_DIR%code/}${LOG_DIR#../}${CUR_DATE}_${DATA_NAME}_${POS_START}_${POS_END}.log"
    else
        echo "[error] Reference file $REF_FILE not found!"
        exit 1
    fi
fi
cd $WORK_DIR
echo "[INFO] Current working directory: $WORK_DIR"
echo "[INFO] Executing: python Run_DREEM.py ../data/$DATA_NAME $RESULTS_DIR $DATA_NAME $REF_NAME $POS_START $POS_END $FASTQ_FLAG $STRUCT_FLAG $CTRL_FLAG" >> "${LOG_DIR}${CUR_DATE}_${DATA_NAME}_${POS_START}_${POS_END}.log"
nohup stdbuf -oL python Run_DREEM.py ../data/$DATA_NAME $RESULTS_DIR $DATA_NAME $REF_NAME $POS_START $POS_END $FASTQ_FLAG $STRUCT_FLAG $CTRL_FLAG >> "${LOG_DIR}${CUR_DATE}_${DATA_NAME}_${POS_START}_${POS_END}.log" 2>&1 &
# -------------------------------------------------------------------
PID=$!
if ps -p $PID > /dev/null; then
    echo "[INFO] Executing: python Run_DREEM.py ../data/$DATA_NAME $RESULTS_DIR $DATA_NAME $REF_NAME $POS_START $POS_END $FASTQ_FLAG $STRUCT_FLAG $CTRL_FLAG"
    echo "[INFO] Run started at [$(date)] with PID: $PID"
    echo "[INFO] Run started at [$(date)] with PID: $PID" >> "${LOG_DIR}${CUR_DATE}_${DATA_NAME}_${POS_START}_${POS_END}.log"
    echo "      ( If you want to stop the process, run: kill -9 $PID )"
    echo "      ( If you want to view terminal output in real time, run: tail -f ${WORK_DIR%code/}${LOG_DIR#../}${CUR_DATE}_${DATA_NAME}_${POS_START}_${POS_END}.log )"
    wait $PID
    PYTHON_STATUS=$?
    if [ $PYTHON_STATUS -eq 0 ]; then
        echo "[INFO] [$(date)] DREEM analysis completed successfully!"
        echo "[INFO] [$(date)] DREEM analysis completed successfully!" >> "${LOG_DIR}${CUR_DATE}_${DATA_NAME}_${POS_START}_${POS_END}.log"
    else
        echo "[error] [$(date)] DREEM analysis failed! Please check the log file in ${WORK_DIR%code/}${LOG_DIR#../} for more information."
        echo "[error] [$(date)] DREEM analysis failed! Please check the log file in ${WORK_DIR%code/}${LOG_DIR#../} for more information." >> "${LOG_DIR}${CUR_DATE}_${DATA_NAME}_${POS_START}_${POS_END}.log"
        exit 1
    fi
else
    echo "[error] [$(date)] DREEM failed to start! Please check the log file in ${WORK_DIR%code/}${LOG_DIR#../} for more information."
    echo "[error] [$(date)] DREEM failed to start! Please check the log file in ${WORK_DIR%code/}${LOG_DIR#../} for more information." >> "${LOG_DIR}${CUR_DATE}_${DATA_NAME}_${POS_START}_${POS_END}.log"
    exit 1
fi