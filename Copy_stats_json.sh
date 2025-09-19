#!/bin/bash
# ===================================================================

EXP_RESULTS=results_20250321_MTS_s4U_APDF_B177R
CONTROL_RESULTS=results_20250321_MTS_s4U_AP_B176R
# -------------------------------------------------------------------
EXP_DATE=$(echo "$EXP_RESULTS" | awk -F'_' '{print $2}')
EXP_NAME=$(echo "$EXP_RESULTS" | cut -d'_' -f3-)

CONTROL_DATE=$(echo "$CONTROL_RESULTS" | awk -F'_' '{print $2}')
CONTROL_NAME=$(echo "$CONTROL_RESULTS" | cut -d'_' -f3-)
# -------------------------------------------------------------------
python /sibcb1/hanshuolab1/wangziyuan/DREEM/code/Copy_stats_json.py $EXP_DATE $EXP_NAME $CONTROL_DATE $CONTROL_NAME