#!/bin/bash
# ===================================================================

source /sibcb/program/install/anaconda3/etc/profile.d/conda.sh
eval "$(conda shell.bash hook)"
conda activate dreem || {
    echo "[error] conda environment 'dreem' not found!";
    exit 1;
}
# -------------------------------------------------------------------
REF_NAME=${REF_NAME}                                  # 参考序列的名称
DATA_NAME=${DATA_NAME}                                # 分析样本的名称
POS_START=${POS_START}                                # 分析的起始位置
POS_END=${POS_END}                                    # 分析的终止位置
FASTQ=${FASTQ}                                        # 是否从fastq文件开始进行分析
CTRL=${CTRL}                                          # 是否为对照组样本
STRUCT=${STRUCT}                                      # 是否使用RNAStructure进行结构预测
CHUNK_SIZE=${CHUNK_SIZE}                              # 每个任务的分析区间大小
# -------------------------------------------------------------------
if [ -z "$REF_NAME" ]; then
  echo "[error] REF_NAME not set!"
  exit 1
fi
if [ -z "$DATA_NAME" ]; then
  echo "[error] DATA_NAME not set!"
  exit 1
fi
if [ -z "$POS_START" ]; then
  echo "[error] POS_START not set!"
  exit 1
fi
if [ -z "$POS_END" ]; then
  echo "[error] POS_END not set!"
  exit 1
fi
if [ -z "$CHUNK_SIZE" ]; then
  echo "[error] CHUNK_SIZE not set!"
  exit 1
fi
# -------------------------------------------------------------------
WORK_DIR=/sibcb1/hanshuolab1/wangziyuan/DREEM/code/
CUR_DATE=$(date +%Y%m%d)
REF_FILE=${REF_NAME}.fasta
# -------------------------------------------------------------------
TASK_START=$((POS_START + (SGE_TASK_ID - 1) * CHUNK_SIZE))
TASK_END=$((TASK_START + CHUNK_SIZE - 1))
if [ "$TASK_END" -gt "$POS_END" ]; then
    TASK_END=$POS_END
fi
# -------------------------------------------------------------------
RESULTS_DIR=../results_${CUR_DATE}_${DATA_NAME}/results_${CUR_DATE}_${DATA_NAME}_${TASK_START}_${TASK_END}/
LOG_DIR=../logs/${CUR_DATE}_${DATA_NAME}_${TASK_START}_${TASK_END}/
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
        echo "[INFO] Bowtie2 index builded for $REF_FILE" >> "${WORK_DIR%code/}${LOG_DIR#../}${CUR_DATE}_${DATA_NAME}_${TASK_START}_${TASK_END}.log"
    else
        echo "[error] Reference file $REF_FILE not found!"
        exit 1
    fi
fi
cd $WORK_DIR

echo "[INFO] Current working directory: $WORK_DIR"
echo "[INFO] Task ID: $SGE_TASK_ID, POS_START: $TASK_START, POS_END: $TASK_END"
echo "[INFO] Executing: python Run_DREEM.py ../data/$DATA_NAME $RESULTS_DIR $DATA_NAME $REF_NAME $TASK_START $TASK_END $FASTQ_FLAG $STRUCT_FLAG $CTRL_FLAG" >> "${LOG_DIR}${CUR_DATE}_${DATA_NAME}_${TASK_START}_${TASK_END}.log"
nohup stdbuf -oL python Run_DREEM.py ../data/$DATA_NAME $RESULTS_DIR $DATA_NAME $REF_NAME $TASK_START $TASK_END $FASTQ_FLAG $STRUCT_FLAG $CTRL_FLAG >> "${LOG_DIR}${CUR_DATE}_${DATA_NAME}_${TASK_START}_${TASK_END}.log" 2>&1 &

PID=$!

if ps -p $PID > /dev/null; then
    echo "[INFO] Executing: python Run_DREEM.py ../data/$DATA_NAME $RESULTS_DIR $DATA_NAME $REF_NAME $TASK_START $TASK_END $FASTQ_FLAG $STRUCT_FLAG $CTRL_FLAG"
    echo "[INFO] Run started at [$(date)] with PID: $PID"
    echo "[INFO] Run started at [$(date)] with PID: $PID" >> "${LOG_DIR}${CUR_DATE}_${DATA_NAME}_${TASK_START}_${TASK_END}.log"
    echo "      ( If you want to stop the process, run: kill -9 $PID )"
    echo "      ( If you want to view terminal output in real time, run: tail -f ${WORK_DIR%code/}${LOG_DIR#../}${CUR_DATE}_${DATA_NAME}_${TASK_START}_${TASK_END}.log )"
    wait $PID
    PYTHON_STATUS=$?
    if [ $PYTHON_STATUS -eq 0 ]; then
        echo "[INFO] [$(date)] DREEM analysis completed successfully for task $SGE_TASK_ID!"
        echo "[INFO] [$(date)] DREEM analysis completed successfully for task $SGE_TASK_ID!" >> "${LOG_DIR}${CUR_DATE}_${DATA_NAME}_${TASK_START}_${TASK_END}.log"
    else
        echo "[error] [$(date)] DREEM analysis failed for task $SGE_TASK_ID! Please check the log file in ${WORK_DIR%code/}${LOG_DIR#../} for more information."
        echo "[error] [$(date)] DREEM analysis failed for task $SGE_TASK_ID! Please check the log file in ${WORK_DIR%code/}${LOG_DIR#../} for more information." >> "${LOG_DIR}${CUR_DATE}_${DATA_NAME}_${TASK_START}_${TASK_END}.log"
        exit 1
    fi
else
    echo "[error] [$(date)] DREEM failed to start for task $SGE_TASK_ID! Please check the log file in ${WORK_DIR%code/}${LOG_DIR#../} for more information."
    echo "[error] [$(date)] DREEM failed to start for task $SGE_TASK_ID! Please check the log file in ${WORK_DIR%code/}${LOG_DIR#../} for more information." >> "${LOG_DIR}${CUR_DATE}_${DATA_NAME}_${TASK_START}_${TASK_END}.log"
    exit 1
fi