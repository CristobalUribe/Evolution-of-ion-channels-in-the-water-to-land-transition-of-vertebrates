#!/usr/bin/env bash
# Portable adaptation: setup generalized; scientific GUIDANCE3 options unchanged.

set -uo pipefail

ROOT="${ROOT:-${GUIDANCE_ROOT:-$(pwd)/guidance3_analysis}}"
INPUT_LIST="${INPUT_LIST:-$ROOT/og_protein_files.txt}"
RESULTS="$ROOT/results"
LOGS="$ROOT/logs"
STATUS="$ROOT/status"
ARCHIVES="$ROOT/incomplete_archives"

G3_BIN="${G3_BIN:-$(command -v guidance3 || true)}"
MAFFT="${MAFFT:-$(command -v mafft || true)}"

MAX_JOBS="${MAX_JOBS:-2}"
PROC_NUM="${PROC_NUM:-4}"
EXPECTED_OGS="${EXPECTED_OGS:-0}"

mkdir -p "$RESULTS" "$LOGS" "$STATUS" "$ARCHIVES"

if [[ ! -x "$G3_BIN" ]]; then
    echo "ERROR: GUIDANCE3 no es ejecutable: $G3_BIN"
    exit 1
fi

if [[ ! -x "$MAFFT" ]]; then
    echo "ERROR: MAFFT no es ejecutable: $MAFFT"
    exit 1
fi

[[ -f "$INPUT_LIST" ]] || { echo "ERROR: lista ausente: $INPUT_LIST" >&2; exit 1; }
[[ "$MAX_JOBS" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: MAX_JOBS debe ser positivo" >&2; exit 1; }
[[ "$PROC_NUM" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: PROC_NUM debe ser positivo" >&2; exit 1; }
[[ "$EXPECTED_OGS" =~ ^[0-9]+$ ]] || { echo "ERROR: EXPECTED_OGS debe ser no negativo" >&2; exit 1; }
count=$(wc -l < "$INPUT_LIST")
(( count > 0 )) || { echo "ERROR: la lista está vacía" >&2; exit 1; }
if (( EXPECTED_OGS > 0 && count != EXPECTED_OGS )); then
    echo "ERROR: se esperaban $EXPECTED_OGS OG y se encontraron $count" >&2
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
