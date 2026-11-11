#!/usr/bin/env bash

# Lanzador reanudable para los modelos branch-site de PAML/codeml.
# Ejecuta alternative y null de cada OG, limita la concurrencia y valida
# cada resultado antes de crear el marcador CODEML_RUN_OK.

set -uo pipefail

PAML_ROOT="${PAML_ROOT:-}"
CODEML="${CODEML:-/home/cloyola/PAML4.10/paml-4.10.10-linux-x86_64/bin/codeml}"
MAX_JOBS="${MAX_JOBS:-80}"
EXPECTED_OGS="${EXPECTED_OGS:-167}"
NICE_LEVEL="${NICE_LEVEL:-10}"

usage() {
    cat <<'EOF'
Uso:
  bash run_guidance3_paml.sh --paml-root DIR [opciones]

Opciones:
  --paml-root DIR      Directorio que contiene alternative/ y null/.
  --codeml ARCHIVO     Ejecutable codeml.
  --max-jobs N         Máximo de procesos codeml simultáneos (default: 80).
  --expected-ogs N     Número esperado de OG (default: 167; 0 desactiva el control).
  --nice N             Prioridad nice de cada codeml (default: 10).
  -h, --help           Mostrar esta ayuda.

También se aceptan las variables PAML_ROOT, CODEML, MAX_JOBS,
EXPECTED_OGS y NICE_LEVEL.
EOF
}

die() {
    printf 'ERROR: %s\n' "$*" >&2
    exit 2
}

while (( $# > 0 )); do
    case "$1" in
        --paml-root)
            (( $# >= 2 )) || die "falta el valor de --paml-root"
            PAML_ROOT="$2"
            shift 2
            ;;
        --codeml)
            (( $# >= 2 )) || die "falta el valor de --codeml"
            CODEML="$2"
            shift 2
            ;;
        --max-jobs)
            (( $# >= 2 )) || die "falta el valor de --max-jobs"
            MAX_JOBS="$2"
            shift 2
            ;;
        --expected-ogs)
            (( $# >= 2 )) || die "falta el valor de --expected-ogs"
            EXPECTED_OGS="$2"
            shift 2
            ;;
        --nice)
            (( $# >= 2 )) || die "falta el valor de --nice"
            NICE_LEVEL="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            die "opción desconocida: $1"
            ;;
    esac
done

[[ -n "$PAML_ROOT" ]] || die "debes indicar --paml-root"
PAML_ROOT="${PAML_ROOT%/}"

[[ "$MAX_JOBS" =~ ^[1-9][0-9]*$ ]] || die "--max-jobs debe ser un entero positivo"
[[ "$EXPECTED_OGS" =~ ^[0-9]+$ ]] || die "--expected-ogs debe ser un entero no negativo"
[[ "$NICE_LEVEL" =~ ^-?[0-9]+$ ]] || die "--nice debe ser un entero"
(( NICE_LEVEL >= -20 && NICE_LEVEL <= 19 )) || die "--nice debe estar entre -20 y 19"

for command_name in awk date find flock grep mkdir mv nice sed sort tail wc; do
    command -v "$command_name" >/dev/null 2>&1 || die "no se encontró el comando: $command_name"
done

[[ -d "$PAML_ROOT/alternative" ]] || die "no existe $PAML_ROOT/alternative"
[[ -d "$PAML_ROOT/null" ]] || die "no existe $PAML_ROOT/null"
[[ -x "$CODEML" ]] || die "codeml no existe o no es ejecutable: $CODEML"
[[ -f "$PAML_ROOT/PAML_INPUTS_READY" ]] || die "falta $PAML_ROOT/PAML_INPUTS_READY"

LAUNCH_LOG="$PAML_ROOT/guidance3_codeml_launcher.log"
LOG_LOCK="$PAML_ROOT/guidance3_codeml_launcher.log.lock"
SUMMARY="$PAML_ROOT/guidance3_codeml_run_summary.tsv"
BATCH_MARKER="$PAML_ROOT/CODEML_BATCH_COMPLETE"

exec 9>"$PAML_ROOT/guidance3_codeml_launcher.lock"
if ! flock -n 9; then
    die "ya hay otro lanzador activo para $PAML_ROOT"
fi

log_message() {
    local timestamp line
    timestamp="$(date --iso-8601=seconds)"
    line="[$timestamp] $*"
    (
        exec 200>"$LOG_LOCK"
        flock -x 200
        printf '%s\n' "$line" >> "$LAUNCH_LOG"
    )
    printf '%s\n' "$line"
}

ctl_value() {
    local ctl_file="$1"
    local wanted_key="$2"
    awk -F '=' -v wanted_key="$wanted_key" '
        {
            key=$1
            gsub(/[[:space:]]/, "", key)
            if (key == wanted_key) {
                value=$2
                sub(/#.*/, "", value)
                gsub(/^[[:space:]]+|[[:space:]]+$/, "", value)
                gsub(/^"|"$/, "", value)
                gsub(/^\047|\047$/, "", value)
                print value
                exit
            }
        }
    ' "$ctl_file"
}

output_path_for() {
    local run_dir="$1"
    local ctl_file="$run_dir/codeml.ctl"
    local output_name
    output_name="$(ctl_value "$ctl_file" outfile)"
    [[ -n "$output_name" ]] || return 1
    [[ "$output_name" != /* && "$output_name" != *".."* ]] || return 1
    printf '%s/%s\n' "$run_dir" "$output_name"
}

last_lnl_line() {
    local output_file="$1"
    grep 'lnL(ntime:' "$output_file" 2>/dev/null | tail -n 1
}

extract_stats() {
    local output_file="$1"
    local line ntime np lnl
    line="$(last_lnl_line "$output_file")"
    [[ -n "$line" ]] || return 1

    ntime="$(sed -nE 's/.*lnL\(ntime:[[:space:]]*([0-9]+)[[:space:]]+np:.*/\1/p' <<< "$line")"
    np="$(sed -nE 's/.*np:[[:space:]]*([0-9]+)\):.*/\1/p' <<< "$line")"
    lnl="$(sed -nE 's/.*\):[[:space:]]*([-+0-9.eE]+).*/\1/p' <<< "$line")"

    [[ "$lnl" =~ ^[-+]?[0-9]+([.][0-9]+)?([eE][-+]?[0-9]+)?$ ]] || return 1
    [[ "$ntime" =~ ^[0-9]+$ ]] || ntime="NA"
    [[ "$np" =~ ^[0-9]+$ ]] || np="NA"
    printf '%s\t%s\t%s\n' "$ntime" "$np" "$lnl"
}

is_complete() {
    local run_dir="$1"
    local output_file stdout_file
    output_file="$(output_path_for "$run_dir")" || return 1
    stdout_file="$run_dir/codeml.stdout"

    [[ -s "$output_file" ]] || return 1
    extract_stats "$output_file" >/dev/null || return 1

    if grep -q 'Time used:' "$output_file" 2>/dev/null; then
        return 0
    fi
    [[ -s "$stdout_file" ]] && grep -q 'Time used:' "$stdout_file" 2>/dev/null
}

write_run_record() {
    local run_dir="$1"
    local model="$2"
    local og="$3"
    local status="$4"
    local started_at="$5"
    local finished_at="$6"
    local duration_s="$7"
    local exit_code="$8"
    local output_file="$9"
    local stats ntime np lnl temp_file

    stats="$(extract_stats "$output_file" 2>/dev/null || true)"
    if [[ -n "$stats" ]]; then
        IFS=$'\t' read -r ntime np lnl <<< "$stats"
    else
        ntime="NA"
        np="NA"
        lnl="NA"
    fi

    temp_file="$run_dir/codeml_run.tsv.tmp.$$"
    {
        printf 'model\tog\tstatus\tstarted_at\tfinished_at\tduration_s\texit_code\tntime\tnp\tlnL\toutput\tstdout\n'
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$model" "$og" "$status" "$started_at" "$finished_at" \
            "$duration_s" "$exit_code" "$ntime" "$np" "$lnl" \
            "$output_file" "$run_dir/codeml.stdout"
    } > "$temp_file"
    mv "$temp_file" "$run_dir/codeml_run.tsv"
}

mark_success() {
    local run_dir="$1"
    local model="$2"
    local og="$3"
    local status="$4"
    local started_at="$5"
    local finished_at="$6"
    local duration_s="$7"
    local exit_code="$8"
    local output_file stats

    output_file="$(output_path_for "$run_dir")" || return 1
    if [[ -f "$run_dir/CODEML_RUN_ERROR" ]]; then
        mv "$run_dir/CODEML_RUN_ERROR" \
            "$run_dir/CODEML_RUN_ERROR.superseded.$(date +%Y%m%dT%H%M%S).$$"
    fi
    write_run_record "$run_dir" "$model" "$og" "$status" "$started_at" \
        "$finished_at" "$duration_s" "$exit_code" "$output_file"
    stats="$(extract_stats "$output_file")" || return 1
    printf 'status=%s\tmodel=%s\tog=%s\tntime_np_lnL=%s\tfinished_at=%s\n' \
        "$status" "$model" "$og" "$stats" "$finished_at" \
        > "$run_dir/CODEML_RUN_OK"
    return 0
}

archive_incomplete_attempt() {
    local run_dir="$1"
    local output_file archive_dir path moved=0
    local -a possible_files

    output_file="$(output_path_for "$run_dir" 2>/dev/null || true)"
    possible_files=(
        "$output_file"
        "$run_dir/codeml.stdout"
        "$run_dir/codeml.pid"
        "$run_dir/codeml_run.tsv"
        "$run_dir/CODEML_RUN_OK"
        "$run_dir/CODEML_RUN_ERROR"
        "$run_dir/rst"
        "$run_dir/rst1"
        "$run_dir/rub"
        "$run_dir/lnf"
        "$run_dir/2NG.dN"
        "$run_dir/2NG.dS"
        "$run_dir/2NG.t"
        "$run_dir/4fold.nuc"
    )

    archive_dir="$run_dir/previous_attempt_$(date +%Y%m%dT%H%M%S)_$$"
    for path in "${possible_files[@]}"; do
        [[ -n "$path" && -e "$path" ]] || continue
        if (( moved == 0 )); then
            mkdir -p "$archive_dir"
            moved=1
        fi
        mv "$path" "$archive_dir/"
    done

    if (( moved == 1 )); then
        printf '%s\n' "$archive_dir"
    fi
}

validate_inputs() {
    local model og run_dir ctl_file seqfile treefile
    local ctl_model nssites fix_omega omega cleandata
    local errors=0

    for og in "${OGS[@]}"; do
        for model in alternative null; do
            run_dir="$PAML_ROOT/$model/$og"
            ctl_file="$run_dir/codeml.ctl"
            if [[ ! -d "$run_dir" || ! -f "$run_dir/PAML_INPUT_OK" || ! -s "$ctl_file" ]]; then
                printf 'Entrada incompleta: %s/%s\n' "$model" "$og" >&2
                (( errors += 1 ))
                continue
            fi

            seqfile="$(ctl_value "$ctl_file" seqfile)"
            treefile="$(ctl_value "$ctl_file" treefile)"
            if [[ -z "$seqfile" || "$seqfile" == /* || "$seqfile" == *".."* || ! -s "$run_dir/$seqfile" ]]; then
                printf 'seqfile inválido: %s/%s (%s)\n' "$model" "$og" "$seqfile" >&2
                (( errors += 1 ))
            fi
            if [[ -z "$treefile" || "$treefile" == /* || "$treefile" == *".."* || ! -s "$run_dir/$treefile" ]]; then
                printf 'treefile inválido: %s/%s (%s)\n' "$model" "$og" "$treefile" >&2
                (( errors += 1 ))
            fi
            output_path_for "$run_dir" >/dev/null || {
                printf 'outfile inválido: %s/%s\n' "$model" "$og" >&2
                (( errors += 1 ))
            }

            ctl_model="$(ctl_value "$ctl_file" model)"
            nssites="$(ctl_value "$ctl_file" NSsites)"
            fix_omega="$(ctl_value "$ctl_file" fix_omega)"
            omega="$(ctl_value "$ctl_file" omega)"
            cleandata="$(ctl_value "$ctl_file" cleandata)"
            if [[ "$ctl_model" != "2" || "$nssites" != "2" || "$cleandata" != "1" ]]; then
                printf 'Controles incompatibles: %s/%s (model=%s NSsites=%s cleandata=%s)\n' \
                    "$model" "$og" "$ctl_model" "$nssites" "$cleandata" >&2
                (( errors += 1 ))
            fi
            if [[ "$model" == "alternative" ]]; then
                if [[ "$fix_omega" != "0" || ! "$omega" =~ ^0?\.4(0*)?$ ]]; then
                    printf 'Control alternativo inválido: %s (fix_omega=%s omega=%s)\n' \
                        "$og" "$fix_omega" "$omega" >&2
                    (( errors += 1 ))
                fi
            elif [[ "$fix_omega" != "1" || ! "$omega" =~ ^1(\.0*)?$ ]]; then
                printf 'Control nulo inválido: %s (fix_omega=%s omega=%s)\n' \
                    "$og" "$fix_omega" "$omega" >&2
                (( errors += 1 ))
            fi
        done
    done

    (( errors == 0 )) || die "se detectaron $errors problemas en los inputs; no se lanzó codeml"
}

run_one() {
    local model="$1"
    local og="$2"
    local run_dir="$PAML_ROOT/$model/$og"
    local output_file started_at finished_at start_epoch end_epoch duration_s rc
    local codeml_pid=""
    local archived=""

    output_file="$(output_path_for "$run_dir")" || {
        log_message "[ERROR] $model/$og: no se pudo determinar outfile"
        return 0
    }

    if is_complete "$run_dir"; then
        finished_at="$(date --iso-8601=seconds)"
        if [[ -f "$run_dir/CODEML_RUN_OK" ]]; then
            if [[ ! -s "$run_dir/codeml_run.tsv" ]]; then
                mark_success "$run_dir" "$model" "$og" "SKIPPED_VALID" \
                    "NA" "$finished_at" "NA" "0" || true
            fi
            log_message "[SKIP] $model/$og ya estaba completo"
        else
            mark_success "$run_dir" "$model" "$og" "RECOVERED" \
                "NA" "$finished_at" "NA" "0" || true
            log_message "[RECOVERED] $model/$og salida completa preexistente"
        fi
        return 0
    fi

    archived="$(archive_incomplete_attempt "$run_dir")"
    if [[ -n "$archived" ]]; then
        log_message "[ARCHIVE] $model/$og intento incompleto -> $archived"
    fi

    started_at="$(date --iso-8601=seconds)"
    start_epoch="$(date +%s)"
    log_message "[START] $model/$og"

    cd "$run_dir" || {
        log_message "[ERROR] $model/$og: no se pudo entrar en $run_dir"
        return 0
    }

    terminate_worker() {
        trap - TERM INT
        if [[ -n "$codeml_pid" ]] && kill -0 "$codeml_pid" 2>/dev/null; then
            kill -TERM "$codeml_pid" 2>/dev/null || true
            wait "$codeml_pid" 2>/dev/null || true
        fi
        finished_at="$(date --iso-8601=seconds)"
        end_epoch="$(date +%s)"
        duration_s=$(( end_epoch - start_epoch ))
        printf 'status=INTERRUPTED\tmodel=%s\tog=%s\tfinished_at=%s\n' \
            "$model" "$og" "$finished_at" > "$run_dir/CODEML_RUN_ERROR"
        log_message "[INTERRUPTED] $model/$og duración=${duration_s}s"
        exit 143
    }
    trap terminate_worker TERM INT

    env \
        OMP_NUM_THREADS=1 \
        OPENBLAS_NUM_THREADS=1 \
        MKL_NUM_THREADS=1 \
        NUMEXPR_NUM_THREADS=1 \
        TF_NUM_INTRAOP_THREADS=1 \
        TF_NUM_INTEROP_THREADS=1 \
        nice -n "$NICE_LEVEL" "$CODEML" codeml.ctl \
        > "$run_dir/codeml.stdout" 2>&1 &
    codeml_pid=$!
    printf '%s\n' "$codeml_pid" > "$run_dir/codeml.pid"

    wait "$codeml_pid"
    rc=$?
    trap - TERM INT

    finished_at="$(date --iso-8601=seconds)"
    end_epoch="$(date +%s)"
    duration_s=$(( end_epoch - start_epoch ))

    if (( rc == 0 )) && is_complete "$run_dir" && \
        mark_success "$run_dir" "$model" "$og" "OK" "$started_at" \
            "$finished_at" "$duration_s" "$rc"; then
        mv "$run_dir/codeml.pid" "$run_dir/codeml.last.pid"
        log_message "[OK] $model/$og duración=${duration_s}s"
    else
        write_run_record "$run_dir" "$model" "$og" "ERROR" "$started_at" \
            "$finished_at" "$duration_s" "$rc" "$output_file"
        printf 'status=ERROR\tmodel=%s\tog=%s\texit_code=%s\tfinished_at=%s\n' \
            "$model" "$og" "$rc" "$finished_at" > "$run_dir/CODEML_RUN_ERROR"
        mv "$run_dir/codeml.pid" "$run_dir/codeml.last.pid"
        log_message "[ERROR] $model/$og exit=$rc duración=${duration_s}s"
    fi

    return 0
}

write_global_summary() {
    local temp_file="$SUMMARY.tmp.$$"
    local og model record
    printf 'model\tog\tstatus\tstarted_at\tfinished_at\tduration_s\texit_code\tntime\tnp\tlnL\toutput\tstdout\n' \
        > "$temp_file"
    for og in "${OGS[@]}"; do
        for model in alternative null; do
            record="$PAML_ROOT/$model/$og/codeml_run.tsv"
            if [[ -s "$record" ]]; then
                tail -n +2 "$record" >> "$temp_file"
            fi
        done
    done
    mv "$temp_file" "$SUMMARY"
}

mapfile -t OGS < <(
    find "$PAML_ROOT/alternative" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' \
        | sort -V
)

(( ${#OGS[@]} > 0 )) || die "no se encontraron directorios OG en alternative/"
if (( EXPECTED_OGS > 0 && ${#OGS[@]} != EXPECTED_OGS )); then
    die "se esperaban $EXPECTED_OGS OG y se encontraron ${#OGS[@]}"
fi

null_og_count="$(find "$PAML_ROOT/null" -mindepth 1 -maxdepth 1 -type d -printf '.' | wc -c)"
if (( null_og_count != ${#OGS[@]} )); then
    die "alternative contiene ${#OGS[@]} OG, pero null contiene $null_og_count"
fi

validate_inputs

stop_all() {
    local -a worker_pids
    trap - TERM INT
    log_message "[STOP] señal recibida; deteniendo trabajadores"
    mapfile -t worker_pids < <(jobs -pr)
    if (( ${#worker_pids[@]} > 0 )); then
        kill -TERM "${worker_pids[@]}" 2>/dev/null || true
        wait "${worker_pids[@]}" 2>/dev/null || true
    fi
    log_message "[STOPPED] lanzador detenido de forma controlada"
    exit 130
}
trap stop_all TERM INT

expected_models=$(( ${#OGS[@]} * 2 ))
log_message "[LOTE INICIADO] OGs=${#OGS[@]} modelos=$expected_models max_jobs=$MAX_JOBS cleandata=1"

active_jobs=0
for og in "${OGS[@]}"; do
    for model in alternative null; do
        run_one "$model" "$og" &
        (( active_jobs += 1 ))
        if (( active_jobs >= MAX_JOBS )); then
            wait -n || true
            (( active_jobs -= 1 ))
        fi
    done
done

while (( active_jobs > 0 )); do
    wait -n || true
    (( active_jobs -= 1 ))
done

trap - TERM INT
write_global_summary

complete_count="$(find "$PAML_ROOT/alternative" "$PAML_ROOT/null" \
    -mindepth 2 -maxdepth 2 -type f -name 'CODEML_RUN_OK' | wc -l)"
error_count="$(find "$PAML_ROOT/alternative" "$PAML_ROOT/null" \
    -mindepth 2 -maxdepth 2 -type f -name 'CODEML_RUN_ERROR' | wc -l)"

if (( complete_count == expected_models && error_count == 0 )); then
    printf 'READY\tOGs=%s\tmodels=%s\tfinished_at=%s\n' \
        "${#OGS[@]}" "$complete_count" "$(date --iso-8601=seconds)" > "$BATCH_MARKER"
    log_message "[LOTE TERMINADO] completas=$complete_count errores=0 resumen=$SUMMARY"
    exit 0
fi

log_message "[LOTE INCOMPLETO] completas=$complete_count/$expected_models errores=$error_count resumen=$SUMMARY"
exit 1
