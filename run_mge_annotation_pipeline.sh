#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<USAGE
Usage:
  $(basename "$0") -i <input_fasta_or_dir> -o <outdir> [-t threads]

Description:
  Run MobileElementFinder, IntegronFinder, and ISEScan on one FASTA
  file or on all FASTA files in a directory, then normalize outputs into TSV.

Required:
  -i   Input FASTA file or directory containing FASTA files
  -o   Output directory

Optional:
  -t   Threads per job (default: 4)
  -m   MobileElementFinder conda env name (default: mefinder_env)
  -n   IntegronFinder conda env name (default: integronfinder_env)
  -s   ISEScan conda env name (default: isescan_env)

Examples:
  $(basename "$0") -i GI285.fasta -o mge_out -t 4
  $(basename "$0") -i fasta_dir/ -o mge_out -t 8
USAGE
}

INPUT=""
OUTDIR=""
THREADS=4
ME_ENV="mefinder_env"
IF_ENV="integronfinder_env"
IS_ENV="isescan_env"

while getopts ":i:o:t:m:n:s:h" opt; do
    case ${opt} in
        i) INPUT="$OPTARG" ;;
        o) OUTDIR="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        m) ME_ENV="$OPTARG" ;;
        n) IF_ENV="$OPTARG" ;;
        s) IS_ENV="$OPTARG" ;;
        h) usage; exit 0 ;;
        \?) echo "[ERROR] Invalid option: -$OPTARG" >&2; usage; exit 1 ;;
        :) echo "[ERROR] Option -$OPTARG requires an argument." >&2; usage; exit 1 ;;
    esac
done

if [[ -z "$INPUT" || -z "$OUTDIR" ]]; then
    usage
    exit 1
fi

if ! command -v conda >/dev/null 2>&1; then
    echo "[ERROR] conda not found in PATH." >&2
    exit 1
fi

mkdir -p "$OUTDIR"
OUTDIR=$(cd "$OUTDIR" && pwd)

PY_MERGER="$(cd "$(dirname "$0")" && pwd)/merge_mge_annotations.py"
if [[ ! -f "$PY_MERGER" ]]; then
    echo "[ERROR] merge_mge_annotations.py not found next to this script." >&2
    exit 1
fi

collect_fastas() {
    local p="$1"
    if [[ -f "$p" ]]; then
        echo "$p"
    elif [[ -d "$p" ]]; then
        find "$p" -maxdepth 1 -type f \( \
            -iname "*.fa" -o -iname "*.fna" -o -iname "*.fasta" -o -iname "*.fas" -o -iname "*.fsa" \
        \) | sort
    else
        return 1
    fi
}

mapfile -t FASTAS < <(collect_fastas "$INPUT")
if [[ ${#FASTAS[@]} -eq 0 ]]; then
    echo "[ERROR] No FASTA files found." >&2
    exit 1
fi

for fasta in "${FASTAS[@]}"; do
    fasta=$(cd "$(dirname "$fasta")" && pwd)/$(basename "$fasta")
    fname=$(basename "$fasta")
    sample=${fname%.*}

    sample_dir="$OUTDIR/$sample"
    mkdir -p "$sample_dir"/{input,logs,mobileelementfinder,integronfinder,isescan,tables}
    ln -sf "$fasta" "$sample_dir/input/$fname"

    echo "[INFO] Processing sample: $sample"

    # 1) MobileElementFinder
    me_prefix="$sample_dir/mobileelementfinder/$sample"
    echo "[INFO]   MobileElementFinder"
    conda run -n "$ME_ENV" mefinder find \
        --contig "$fasta" \
        --gff \
        -t "$THREADS" \
        "$me_prefix" \
        > "$sample_dir/logs/mobileelementfinder.stdout.log" \
        2> "$sample_dir/logs/mobileelementfinder.stderr.log"

    # 2) IntegronFinder
    # IntegronFinder docs note that increasing CPU too much (usually >4)
    # may decrease performance, so cap here at 4.
    if_cpu="$THREADS"
    if [[ "$if_cpu" -gt 4 ]]; then
        if_cpu=4
    fi
    echo "[INFO]   IntegronFinder"
    conda run -n "$IF_ENV" integron_finder \
        "$fasta" \
        --local-max \
        --linear \
        --cpu "$if_cpu" \
        --outdir "$sample_dir/integronfinder" \
        > "$sample_dir/logs/integronfinder.stdout.log" \
        2> "$sample_dir/logs/integronfinder.stderr.log"

    # 3) ISEScan
    echo "[INFO]   ISEScan"
    conda run -n "$IS_ENV" isescan.py \
        --seqfile "$fasta" \
        --output "$sample_dir/isescan" \
        --nthread "$THREADS" \
        > "$sample_dir/logs/isescan.stdout.log" \
        2> "$sample_dir/logs/isescan.stderr.log"
done

echo "[INFO] Normalizing outputs"
conda run -n "$IF_ENV" python "$PY_MERGER" --root "$OUTDIR"

echo "[INFO] Done"
echo "[INFO] Combined table: $OUTDIR/combined/mge_annotations.all.tsv"
echo "[INFO] BED-like table: $OUTDIR/combined/mge_annotations.all.bedlike.tsv"
