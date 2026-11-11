#!/usr/bin/env bash

set -uo pipefail

ROOT="/home/cloyola/PAML4.10/Guidance3_LINSI_all167"
INPUT_LIST="$ROOT/all167_inputs_G3_LINSI_by_size.txt"
RESULTS="$ROOT/results"
LOGS="$ROOT/logs"
STATUS="$ROOT/status"
ARCHIVES="$ROOT/incomplete_archives"

G3_BIN="/home/cloyola/anaconda3/envs/guidance3/bin/guidance3"
MAFFT="/home/cloyola/anaconda3/envs/guidance2/bin/mafft"

MAX_JOBS=20
PROC_NUM=4

mkdir -p "$RESULTS" "$LOGS" "$STATUS" "$ARCHIVES"

if [[ ! -x "$G3_BIN" ]]; then
    echo "ERROR: GUIDANCE3 no es ejecutable: $G3_BIN"
    exit 1
fi

if [[ ! -x "$MAFFT" ]]; then
    echo "ERROR: MAFFT no es ejecutable: $MAFFT"
    exit 1
fi

if [[ $(wc -l < "$INPUT_LIST") -ne 167 ]]; then
    echo "ERROR: la lista no contiene 167 entradas."
    exit 1
fi

while IFS= read -r input; do
    if [[ ! -s "$input" ]]; then
        echo "ERROR: entrada ausente o vacía: $input"
        exit 1
    fi
done < "$INPUT_LIST"

run_one() {
    local input="$1"
    local base og out log status_file tmp_status
    local start_iso end_iso start_epoch end_epoch elapsed
    local rc state detail child stamp archive

    base=$(basename "$input")
    og=${base%%.*}

    out="$RESULTS/${og}_guidance"
    log="$LOGS/${og}.log"
    status_file="$STATUS/${og}.tsv"

    if [[ -f "$out/ENDS_OK" ]]; then
        printf '%s\t%s\t-\t-\t0\tOK\t0\tENDS_OK preexistente\n' \
          "$og" "$input" > "$status_file"
        echo "[SKIP] $og ya estaba completo"
        return 0
    fi

    if [[ -d "$out" ]]; then
        stamp=$(date +%Y%m%d_%H%M%S)
        archive="$ARCHIVES/${og}_guidance_${stamp}"
        mv "$out" "$archive"

        if [[ -f "$log" ]]; then
            mv "$log" "$ARCHIVES/${og}_${stamp}.log"
        fi

        echo "[ARCHIVADO] $og -> $archive"
    fi

    start_iso=$(date -Is)
    start_epoch=$(date +%s)

    echo "[INICIO] $og $start_iso"

    nice -n 10 env \
      OMP_NUM_THREADS=1 \
      OPENBLAS_NUM_THREADS=1 \
      MKL_NUM_THREADS=1 \
      NUMEXPR_NUM_THREADS=1 \
      TF_NUM_INTRAOP_THREADS=1 \
      TF_NUM_INTEROP_THREADS=1 \
      "$G3_BIN" \
        --seqFile "$input" \
        --msaProgram MAFFT_LINSI \
        --seqType aa \
        --outDir "$out" \
        --program GUIDANCE3 \
        --bootstraps 100 \
        --disableConvergence \
        --colCutoff 0.93 \
        --outOrder as_input \
        --proc_num "$PROC_NUM" \
        --mafft "$MAFFT" \
        > "$log" 2>&1 &

    child=$!

    trap '
        kill -TERM "$child" 2>/dev/null || true
        wait "$child" 2>/dev/null || true
        exit 143
    ' TERM INT

    wait "$child"
    rc=$?

    trap - TERM INT

    end_epoch=$(date +%s)
    end_iso=$(date -Is)
    elapsed=$((end_epoch-start_epoch))

    if [[ $rc -eq 0 && -f "$out/ENDS_OK" ]]; then
        state="OK"
        detail="GUIDANCE3 completo"
        echo "[OK] $og duración=${elapsed}s"
    else
        state="ERROR"
        detail="exit=$rc; revisar $log"
        echo "[ERROR] $og exit=$rc duración=${elapsed}s"
    fi

    tmp_status="${status_file}.tmp.$$"

    printf '%s\t%s\t%s\t%s\t%d\t%s\t%d\t%s\n' \
      "$og" "$input" "$start_iso" "$end_iso" "$elapsed" \
      "$state" "$rc" "$detail" > "$tmp_status"

    mv "$tmp_status" "$status_file"
}

stop_all() {
    local -a children=()

    trap - TERM INT
    mapfile -t children < <(jobs -pr)

    if ((${#children[@]} > 0)); then
        echo "[DETENCIÓN] terminando trabajos activos..."
        kill -TERM "${children[@]}" 2>/dev/null || true
        wait
    fi

    exit 143
}

trap stop_all TERM INT

echo "[LOTE INICIADO] $(date -Is); simultáneos=$MAX_JOBS; proc_num=$PROC_NUM"

running=0

while IFS= read -r input; do
    run_one "$input" &
    ((running+=1))

    if ((running >= MAX_JOBS)); then
        wait -n || true
        ((running-=1))
    fi
done < "$INPUT_LIST"

wait
trap - TERM INT

SUMMARY="$ROOT/timing_summary.tsv"

{
    printf 'OG\tentrada\tinicio\tfin\tsegundos\testado\texit_code\tdetalle\n'
    cat "$STATUS"/*.tsv 2>/dev/null | sort -V
} > "$SUMMARY"

ok=$(awk -F '\t' 'NR>1 && $6=="OK" {n++} END {print n+0}' "$SUMMARY")
errors=$(awk -F '\t' 'NR>1 && $6=="ERROR" {n++} END {print n+0}' "$SUMMARY")

echo "[LOTE TERMINADO] $(date -Is); completas=$ok; errores=$errors"
