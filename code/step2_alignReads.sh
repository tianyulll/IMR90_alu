# Step 1 
# Remove reads do not have ITR
INDIR=fastq
OUTDIR=itrAlign
mkdir -p "$OUTDIR"

ITR="tctgcgcgctcgctcgctca"  
for R2 in "$INDIR"/*_R2_001.fastq.gz; do

    sample=$(basename "$R2" _R2_001.fastq.gz)
    R1="${INDIR}/${sample}_R1_001.fastq.gz"

    echo "Processing $sample"

    summary=$(cutadapt \
      -g "^${ITR15}" \
      -e 1 -O 15 \
      --no-indels \
      --discard-untrimmed \
      --pair-filter=any \
      -o "$OUTDIR/${sample}_R2.fastq.gz" \
      -p "$OUTDIR/${sample}_R1.fastq.gz" \
      "$R2" "$R1" 2>&1)

    echo "$summary" | awk '
        /Total read pairs processed:/ ||
        /Read 1 with adapter:/ ||
        /Pairs written \(passing filters\):/
    '
    echo
    echo
done

minimap2 -t 8 -d IMR90_hap1.mmi IMR90_hap1.fasta.gz

# Step2: Align to IMR90 haplotype 1 reference
mkdir -p alnIMR90
ref=genome/IMR90_hap1.mmi

for r2 in itrAlign/*_R2.fastq.gz; do
  base=$(basename "$r2" _R2.fastq.gz)
  r1="itrAlign/${base}_R1.fastq.gz"

  echo "[$(date '+%F %T')] Aligning R1 (single-end) to IMR90: $base"
  minimap2 -t 8 -ax sr $ref $r1 \
  | samtools sort -@ 8 -o "alnIMR90/${base}.R1.IMR90.bam" -
  samtools index "alnIMR90/${base}.R1.IMR90.bam"

  echo "[$(date '+%F %T')] Aligning R2 (single-end) to IMR90: $base"
  minimap2 -t 8 -ax sr $ref $r2 \
  | samtools sort -@ 8 -o "alnIMR90/${base}.R2.IMR90.bam" -
  samtools index "alnIMR90/${base}.R2.IMR90.bam"
done

