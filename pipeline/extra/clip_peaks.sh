#!/usr/bin/env bash
# PureCLIP binding regions of one reader in one species, from public CLIP libraries, with nf-core/clipseq 1.0.0
# Meant for a compute node of a cluster (UCR HPCC, SLURM; see clip_peaks.sbatch): FASTQ files and the Nextflow work
# directory live in node-local scratch and are deleted when the call ends; the reference and container image too,
# unless CACHE names a directory that keeps them between calls
# usage: clip_peaks.sh SELECTION ROW OUT   (ROW: row of data/raw/clip_standard/selection.tsv, counted without the
#        header)
#        INDEX_ONLY=1 clip_peaks.sh SELECTION ROW OUT   builds the STAR index of the species of the row with the
#        pipeline (its own index step, on 100000 reads of the library) into CACHE and stops: the step needs more memory
#        than the other steps, so it runs once per node, before the array
# writes OUT/{species}_{reader}_{library}.bed.gz (the PureCLIP binding regions of each IP library, as the pipeline
# writes them) and OUT/{species}_{reader}.provenance.tsv (read counts of the pipeline's QC tables, tools and versions)
set -euo pipefail

SELECTION=$1
ROW=$2
OUT=$3
THREADS=${THREADS:-${SLURM_CPUS_PER_TASK:-8}}
MEMORY_GB=${MEMORY_GB:-44}  # for Nextflow and the pipeline's resource limits
READ_LIMIT=${READ_LIMIT:-}  # reads kept per FASTQ file, for tests
INDEX_ONLY=${INDEX_ONLY:-}
[ -n "$INDEX_ONLY" ] && READ_LIMIT=100000
WORK=${WORK:-$(mktemp -d "${TMPDIR:-/tmp}/clip_peaks.XXXXXX")}
CACHE=${CACHE:-$WORK/cache}  # tools, container image and genomes; shared by the calls on a node when CACHE is set
mkdir -p "$WORK/fastq" "$CACHE" "$OUT"
trap 'rm -rf "$WORK"' EXIT

field() { awk -F'\t' -v row="$ROW" -v name="$1" \
  'NR == 1 {for (i = 1; i <= NF; i++) col[$i] = i; next} NR == row + 1 {print $col[name]}' "$SELECTION"; }
READER=$(field reader)
SPECIES=$(field species)
UMI=$(field umi)
echo "$SPECIES $READER: $(field accession), $(field method), UMI $UMI"
[ -n "$INDEX_ONLY" ] && [ -d "$CACHE/$SPECIES/index" ] && { echo "index exists"; exit 0; }
[ -s "$OUT/${SPECIES}_${READER}.provenance.tsv" ] && [ -z "$INDEX_ONLY" ] && { echo "done already"; exit 0; }

# Nextflow 22.10.6 (the pipeline is written in DSL1) and Apptainer 1.5.4 pinned; micromamba, aria2 and pigz unpinned
(
  flock 9
  if ! [ -x "$CACHE/env/bin/nextflow" ]; then
    curl -sL https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xj -C "$CACHE" bin/micromamba
    MAMBA_ROOT_PREFIX=$CACHE/mamba "$CACHE/bin/micromamba" create -y -q -p "$CACHE/env" -c conda-forge -c bioconda \
      nextflow=22.10.6 openjdk=17 apptainer=1.5.4 aria2 pigz
  fi
) 9> "$CACHE/env.lock"
export PATH=$CACHE/env/bin:$PATH
export APPTAINER_CACHEDIR=$CACHE/apptainer APPTAINER_TMPDIR=$WORK/tmp NXF_HOME=$WORK/nextflow \
  NXF_SINGULARITY_CACHEDIR=$CACHE/images
export NXF_OPTS="-Xms512m -Xmx2g"
mkdir -p "$APPTAINER_TMPDIR" "$NXF_SINGULARITY_CACHEDIR"

# release 1.0.0 of the pipeline (cloned: Nextflow's own download uses the GitHub API, whose limit the cluster
# addresses reach), its container image, and the GENCODE genome of the repository (human v44 on GRCh38, mouse vM25 on
# GRCm38)
case $SPECIES in
  human) BASE=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44
         FASTA=GRCh38.primary_assembly.genome.fa.gz; GTF=gencode.v44.basic.annotation.gtf.gz ;;
  mouse) BASE=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25
         FASTA=GRCm38.primary_assembly.genome.fa.gz; GTF=gencode.vM25.basic.annotation.gtf.gz ;;
esac
(
  flock 9
  [ -s "$NXF_SINGULARITY_CACHEDIR/nfcore-clipseq-1.0.0.img" ] ||
    apptainer pull -F "$NXF_SINGULARITY_CACHEDIR/nfcore-clipseq-1.0.0.img" docker://nfcore/clipseq:1.0.0 > /dev/null
  [ -d "$CACHE/clipseq" ] || git clone -q --depth 1 --branch 1.0.0 https://github.com/nf-core/clipseq "$CACHE/clipseq"
  mkdir -p "$CACHE/$SPECIES"
  for file in $FASTA $GTF; do
    [ -s "$CACHE/$SPECIES/$file" ] ||
      aria2c -q -x4 -s4 --max-tries=10 --retry-wait=15 -d "$CACHE/$SPECIES" -o "$file" "$BASE/$file"
  done
) 9> "$CACHE/$SPECIES.lock"

# IP libraries as single-end FASTQ files ({library}_reads.fastq.gz, apart from the {sample}.fastq.gz link of the
# pipeline): runs of a library joined; the UMI of the read name moved to the end of the read id (umi_tools reads it
# there)
name_umi() {
  pigz -dc | awk -v layout="$UMI" 'NR % 4 == 1 {
    umi = ""
    if (layout == "name_prefix" && match($1, /^@[ACGTN]+:/)) {
      umi = substr($1, 2, RLENGTH - 2); $1 = "@" substr($1, RLENGTH + 1) }
    else if (layout == "name_suffix" && match($NF, /_[ACGTN]+(\/[12])?$/)) {
      umi = $NF; sub(/\/[12]$/, "", umi); sub(/^.*_/, "", umi); sub(/_[ACGTN]+(\/[12])?$/, "", $NF) }
    else if (layout == "name_at" && match($2, /^[ACGTN]+@/)) {
      umi = substr($2, 1, RLENGTH - 1); $2 = substr($2, RLENGTH + 1) }
    if (umi == "") { print "no UMI in read name: " $0 > "/dev/stderr"; exit 1 }
    sub(/\/[12]$/, "", $1)
    print $1 "_" umi (NF > 1 ? " " $2 : "")
    next } { print }' | pigz
}
echo "sample,fastq" > "$WORK/samplesheet.csv"
LIBRARIES=$(field ip_libraries | tr ';' ' ')
[ -n "$INDEX_ONLY" ] && LIBRARIES=${LIBRARIES%% *}
for library in $LIBRARIES; do
  : > "$WORK/fastq/$library.raw"
  for url in $(field fastq | tr ';' '\n' | awk -F'|' -v l="$library" '$1 == l {print $2}'); do
    if [ -n "$READ_LIMIT" ]; then
      { curl -sL "$url" | pigz -dc | head -n $((4 * READ_LIMIT)) | pigz; } >> "$WORK/fastq/$library.raw" || true
    else
      aria2c -q -x4 -s4 --max-tries=10 --retry-wait=15 -d "$WORK" -o piece.gz "$url"
      cat "$WORK/piece.gz" >> "$WORK/fastq/$library.raw"
      rm "$WORK/piece.gz"
    fi
  done
  case $UMI in
    name_*) name_umi < "$WORK/fastq/$library.raw" > "$WORK/fastq/${library}_reads.fastq.gz"
            rm "$WORK/fastq/$library.raw" ;;
    *) mv "$WORK/fastq/$library.raw" "$WORK/fastq/${library}_reads.fastq.gz" ;;
  esac
  echo "$library,$WORK/fastq/${library}_reads.fastq.gz" >> "$WORK/samplesheet.csv"
done

# nf-core/clipseq with its defaults, except the genome (and its STAR index when built), the small RNA premapping of the
# species, PureCLIP as the peak caller, the UMI layout of the dataset, and the resource limits of the job
INDEX=$CACHE/$SPECIES/index
# Two settings the data force. The pipeline's default for the PureCLIP chromosomes (-iv 'all') is no chromosome name
# and stops PureCLIP with a NaN parameter, so the option is left empty: PureCLIP then learns on all contigs, its own
# default.
# preseq (library complexity, read only by MultiQC) stops on libraries without duplicates ("too many defects in the
# approximation") and with it the whole run, so its failure is ignored
printf "executor { name = 'local'; cpus = %s; memory = '%s GB' }\n"\
"params.pureclip_iv = ''\nprocess { withName: preseq { errorStrategy = 'ignore' } }\n" \
  "$THREADS" "$MEMORY_GB" > "$WORK/resources.config"
args=(--input "$WORK/samplesheet.csv" --fasta "$CACHE/$SPECIES/$FASTA" --gtf "$CACHE/$SPECIES/$GTF"
  --smrna_org "$SPECIES" --max_cpus "$THREADS" --max_memory "${MEMORY_GB}.GB" --max_time 48.h --outdir "$WORK/results")
if [ -n "$INDEX_ONLY" ]; then args+=(--save_index); else args+=(--peakcaller pureclip); fi
[ -d "$INDEX" ] && [ -z "$INDEX_ONLY" ] && args+=(--star_index "$INDEX")
PIPELINE=$CACHE/clipseq
case $UMI in
  none) args+=(--deduplicate false)
    # Reads without UMI cannot be deduplicated, and release 1.0.0 stops on its --deduplicate false branch: the branch
    # hands the aligned reads (one DSL1 channel, which is read once) to three processes and leaves the PureCLIP input
    # channel undefined, and the QC script needs the folder of deduplication logs. The branch gets the three copies
    # that the deduplication branch outputs and the QC step an empty folder, in a copy of the pipeline made for this
    # call
    PIPELINE=$WORK/clipseq
    cp -r "$CACHE/clipseq" "$PIPELINE"
    sed -i -e 's/^    ch_dedup = ch_aligned$/    ch_aligned.into { ch_dedup; ch_dedup_pureclip; ch_dedup_rseqc }/' \
      -e '/^    ch_dedup_rseqc = ch_aligned$/d' \
      -e 's/^    clip_qc.py \$clip_qc_args$/    mkdir -p dedup\n&/' "$PIPELINE/main.nf"
    grep -q 'ch_aligned.into { ch_dedup; ch_dedup_pureclip; ch_dedup_rseqc }' "$PIPELINE/main.nf"
    grep -q '^    mkdir -p dedup$' "$PIPELINE/main.nf" ;;
  name_*) args+=(--umi_separator _) ;;
  *) args+=(--move_umi "$UMI" --umi_separator _) ;;
esac
(cd "$WORK" && nextflow run "$PIPELINE" -profile singularity -c "$WORK/resources.config" -w "$WORK/work" "${args[@]}")
if [ -n "$INDEX_ONLY" ]; then
  rm -rf "$INDEX"
  mv "$WORK"/results/STAR_index/STAR_* "$INDEX"
  exit 0
fi

# binding regions of the libraries, and provenance from the pipeline's QC tables
QC=$WORK/results/clipqc
TOOLS=$(awk -F'\t' '{printf "%s%s %s", (NR > 1 ? "; " : ""), $1, $2}' \
  "$WORK/results/pipeline_info/software_versions.csv")
{
  printf 'reader\tspecies\tsource\taccession\tcell\tmethod\tlibrary\tinput_reads\tsmallrna_reads\tgenome_reads\t'\
'deduplicated_reads\tcrosslink_sites\tbinding_regions\ttools\n'
  for library in $LIBRARIES; do
    cp "$WORK/results/pureclip/$library.8nt.peaks.bed.gz" "$OUT/${SPECIES}_${READER}_${library}.bed.gz"
    reads=$(awk -F'\t' -v l="$library" '$1 == l {print $2 "\t" $3 "\t" $4}' "$QC/mapping_metrics.tsv")
    deduplicated=$(awk -F'\t' -v l="$library" '$1 == l {print $3}' "$QC/dedup_metrics.tsv" 2> /dev/null || echo NA)
    sites=$(awk -F'\t' -v l="$library" '$1 == l {print $4}' "$QC/peaks_metrics.tsv")
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$READER" "$SPECIES" "$(field source)" \
      "$(field accession)" "$(field cell)" "$(field method)" "$library" "$reads" "${deduplicated:-NA}" "$sites" \
      "$(pigz -dc "$OUT/${SPECIES}_${READER}_${library}.bed.gz" | wc -l)" "$TOOLS"
  done
} > "$OUT/${SPECIES}_${READER}.provenance.tsv"
