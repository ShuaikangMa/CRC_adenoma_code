#!/bin/bash

# ======= Parameter settings =======
SEQ_DIR="/root/autodl-tmp/seq"              # Root directory of SRA data, containing subdirectories starting with SRR
OUT_DIR="./fastq_output"                    # Output directory for decompressed FASTQ files
TMP_DIR="/root/autodl-tmp/middle"            # Temporary directory
THREADS=64                                  # Number of threads used for each decompression (do not set too high)

# ======= Prepare directories =======
mkdir -p "$OUT_DIR"
mkdir -p "$TMP_DIR"

# ======= Main processing workflow =======
for srr_dir in "$SEQ_DIR"/SRR*; do
  if [[ -d "$srr_dir" ]]; then
    sra_file=$(find "$srr_dir" -type f -name "*.sra" | head -n 1)

    if [[ -f "$sra_file" ]]; then
      sra_id=$(basename "$srr_dir")
      echo "============================"
      echo ">>> Processing $sra_id ..."
      echo ">>> SRA file path: $sra_file"

      fasterq-dump "$sra_file" \
        -e "$THREADS" \
        -O "$OUT_DIR" \
        --temp "$TMP_DIR" \
        --progress

      if [[ $? -eq 0 ]]; then
        echo "✓ Decompression finished, removing directory: $srr_dir"
        rm -rf "$srr_dir"
      else
        echo "✗ Decompression failed, keeping directory: $srr_dir for debugging"
      fi

    else
      echo "⚠ No .sra file found in directory $srr_dir, skipping"
    fi
  fi
done

echo "✓ All SRR datasets have been successfully decompressed"


INPUT_DIR="./fastq_output"
OUTPUT_DIR="/root/autodl-tmp/fastq_clean"
THREADS_PER_JOB=16
MAX_PARALLEL_JOBS=10

mkdir -p "$OUTPUT_DIR"

export INPUT_DIR OUTPUT_DIR THREADS_PER_JOB

process_sample() {
    fq="$1"
    fname=$(basename "$fq")
    sample=${fname%.fastq}

    echo ">>> Processing: $sample"
    
    fastp \
     -i "$fq" \
     --length_required 75 \
     --dedup \
     --qualified_quality_phred 20 \
     --thread "$THREADS_PER_JOB" \
     --out1 "$OUTPUT_DIR/clean_${sample}.fastq" \
     -h "$OUTPUT_DIR/report_${sample}.html" \
     -j "$OUTPUT_DIR/report_${sample}.json" \
     2>&1 | tee "$OUTPUT_DIR/log_${sample}.txt"

    # Check whether fastp executed successfully; if successful, delete original file
    if [ $? -eq 0 ]; then
        echo "✓ $sample filtering completed, removing original file $fq"
        rm -f "$fq"
    else
        echo "✗ $sample filtering failed, keeping original file $fq for inspection"
    fi
}

export -f process_sample

find "$INPUT_DIR" -name "*.fastq" | parallel -j "$MAX_PARALLEL_JOBS" process_sample {}

echo "✓ All samples have been processed by fastp in parallel. Results saved in: $OUTPUT_DIR"


# === Parameter settings ===
INPUT_DIR="/root/autodl-tmp/fastq_clean"               # Input FASTQ directory
OUTPUT_DIR="/root/autodl-tmp/dehosted_fastq"           # Output directory
BOWTIE2_INDEX="/root/autodl-tmp/index/genome_index"    # Bowtie2 index prefix (without .bt2 suffix)
THREADS_PER_JOB=16                                     # Threads per process
MAX_PARALLEL_JOBS=10                                   # Maximum number of parallel jobs (adjust based on CPU)

mkdir -p "$OUTPUT_DIR"

export INPUT_DIR OUTPUT_DIR BOWTIE2_INDEX THREADS_PER_JOB

# === Sample processing function ===
process_sample() {
    fq="$1"
    fname=$(basename "$fq")
    sample=${fname%.fastq}

    echo ">>> Removing host contamination for: $sample"

    # Step 1: Align reads to host reference genome
    bowtie2 -x "$BOWTIE2_INDEX" \
        -U "$fq" \
        -p "$THREADS_PER_JOB" \
        --very-sensitive-local \
        -S "$OUTPUT_DIR/${sample}.sam" \
        2> "$OUTPUT_DIR/${sample}_bowtie2.log"

    if [[ ! -s "$OUTPUT_DIR/${sample}.sam" ]]; then
        echo "✗ $sample failed: ${sample}.sam not generated. Check bowtie2 log: ${sample}_bowtie2.log"
        return 1
    fi

    # Step 2: Extract unmapped (non-host) reads
    samtools fastq -f 4 "$OUTPUT_DIR/${sample}.sam" > "$OUTPUT_DIR/clean_${sample}.fastq"

    # Optional: remove intermediate SAM file to save disk space
    rm "$OUTPUT_DIR/${sample}.sam"

    echo "✓ $sample processing completed. Output: clean_${sample}.fastq"
}

export -f process_sample

# === Parallel processing of all FASTQ files ===
find "$INPUT_DIR" -name "*.fastq" | parallel -j "$MAX_PARALLEL_JOBS" process_sample {}

echo "✓ All samples have completed host decontamination. Output directory: $OUTPUT_DIR"


INPUT_DIR="/root/autodl-tmp/dehosted_fastq"
OUTPUT_DIR="/root/autodl-tmp/end"
DB_DIR="/root/autodl-tmp/database"
THREADS=64  # Adjust according to server resources

mkdir -p "$OUTPUT_DIR"
cd "$INPUT_DIR" || exit 1

# Iterate over all FASTQ files matching naming pattern
for file in clean_clean_SRR*.fastq; do
  [ -e "$file" ] || continue  # Skip if no matching file exists

  # Extract sample name (without extension)
  sample_name=$(basename "$file" .fastq)

  echo ">>> Processing: $file"

  # Kraken2 taxonomic classification
  kraken2 \
    --db "$DB_DIR" \
    --use-names \
    --report-zero-counts \
    --report "$OUTPUT_DIR/${sample_name}.kraken.report" \
    --threads "$THREADS" \
    "$file"

  if [ $? -eq 0 ]; then
    echo "✓ Kraken2 classification completed, starting Bracken..."

    # Bracken abundance estimation
    bracken \
      -d "$DB_DIR" \
      -i "$OUTPUT_DIR/${sample_name}.kraken.report" \
      -o "$OUTPUT_DIR/${sample_name}.bracken.out" \
      -w "$OUTPUT_DIR/${sample_name}.bracken.kreport" \
      -l S

    if [ $? -eq 0 ]; then
      echo "✓ Bracken completed, removing original file $file ..."
      rm -f "$file"
    else
      echo "✗ Bracken failed, keeping original file $file for troubleshooting"
    fi

  else
    echo "✗ Kraken2 failed for $file"
  fi
done

echo "✓ All samples have completed taxonomic classification and abundance estimation. Results saved in: $OUTPUT_DIR"
