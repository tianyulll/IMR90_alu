# Step3: Filter for uniquely mapped read pairs
OUTDIR=alu_calls
BAMDIR=alnIMR90
MAPQ=50
mkdir -p "$OUTDIR"


for r1 in "$BAMDIR"/*.R1.IMR90.bam; do
  base=$(basename "$r1" .R1.IMR90.bam)
  r2="$BAMDIR/$base.R2.IMR90.bam"
  [[ -f "$r2" ]] || { echo "Missing $r2" >&2; continue; }

  echo "[$(date '+%F %T')] $base (R1 MAPQ >= $MAPQ → subset R1 & R2)"

  ids=$(mktemp)

  # Collect read IDs from R1: mapped, primary, non-supplementary, MAPQ>=threshold
  samtools view -q "$MAPQ" -F 4 -F 256 -F 2048 "$r1" \
    | cut -f1 | sort -u > "$ids"

  [[ -s "$ids" ]] || { echo "  No reads pass"; rm -f "$ids"; continue; }

  # Subset R1 and R2 BAMs using those read IDs
  samtools view -b -N "$ids" "$r1" | samtools sort -o "$OUTDIR/$base.R1.mapq$MAPQ.bam" -
  samtools index "$OUTDIR/$base.R1.mapq$MAPQ.bam"

  samtools view -b -N "$ids" "$r2" | samtools sort -o "$OUTDIR/$base.R2.fromR1.mapq$MAPQ.bam" -
  samtools index "$OUTDIR/$base.R2.fromR1.mapq$MAPQ.bam"

  rm -f "$ids"
done


# Step4
# Look for Alu in R2; Annotate with coordinates from R1
ALU_BED=genome/aluRef/IMR90_hap1.Alu.bed
OUT=alu_hap1_results.tsv

# output columns:
# sample, read_id,
# R1_chr, R1_start, R1_end, R1_strand, R1_mapq,
# R2_strand,
# alu_type, alu_chr, alu_start, alu_end, alu_strand
echo -e "sample\tread_id\tR1_chr\tR1_start\tR1_end\tR1_strand\tR1_mapq\tR2_strand\talu_type\talu_chr\talu_start\talu_end\talu_strand" > "$OUT"

for r2bam in alu_calls/*.R2.fromR1.mapq*.bam; do
  fname=$(basename "$r2bam")
  sample=${fname%%.*}  # everything before first "."
  echo "Processing $sample"

  # Derive matching R1 BAM name: base.R1.mapqX.bam
  r1bam="${r2bam/.R2.fromR1/.R1}"
  [[ -f "$r1bam" ]] || { echo "Missing R1 BAM for $r2bam -> expected $r1bam" >&2; continue; }

  # Build R1 lookup: read_id -> chr,start,end,strand,mapq (BED from R1 BAM)
  r1_lut=$(mktemp)
  samtools view -b -F 4 -F 256 -F 2048 "$r1bam" \
  | bedtools bamtobed -i - \
  | awk 'BEGIN{OFS="\t"}{print $4, $1, $2, $3, $6, $5}' \
  > "$r1_lut"
  # cols: read_id chr start end strand mapq

  # Intersect R2 with Alu, then attach R1 coords by read_id
  samtools view -b -F 4 -F 256 -F 2048 "$r2bam" \
  | bedtools bamtobed -i - \
  | bedtools intersect -wa -wb -a - -b "$ALU_BED" \
  | awk -v sample="$sample" 'BEGIN{OFS="\t"}
      FNR==NR { r1_chr[$1]=$2; r1_s[$1]=$3; r1_e[$1]=$4; r1_str[$1]=$5; r1_q[$1]=$6; next }
      {
        # read BED (1–6): chr start end read_id mapq strand
        # alu  BED (7–12): chr start end alu_type score strand
        read_id=$4
        r2_strand=$6
        alu_type=$10
        alu_chr=$7; alu_s=$8; alu_e=$9; alu_str=$12

        if (read_id in r1_chr) {
          print sample, read_id,
                r1_chr[read_id], r1_s[read_id], r1_e[read_id], r1_str[read_id], r1_q[read_id],
                r2_strand,
                alu_type, alu_chr, alu_s, alu_e, alu_str
        }
      }' "$r1_lut" - \
  >> "$OUT"

  rm -f "$r1_lut"
done
