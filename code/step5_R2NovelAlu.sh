
# ids is R2 with alu signature
# Overlap with R1 alignment to find uniquely mapped reads

BAMDIR=alu_calls
ids=data/r2_alu_sig.tsv
ALU_BED=genome/aluRef/IMR90_hap1.Alu.bed
OUT=data/r2_alu_res_v1.tsv
genome=genome/IMR90_hap1.genome

rm $OUT
echo -e "sample\tread_id\tR1_chr\tR1_start\tR1_end\tR1_strand\tR1_mapq\talu_type\talu_chr\talu_start\talu_end\talu_strand\tdist_bp" > "$OUT"

for r1 in $BAMDIR/*R1.mapq50.bam; do
	
	fname=$(basename "$r1")
	sample=${fname%%.*}   

	echo "[$(date '+%F %T')] $sample "
	
	samtools view -b -N $ids $r1 \
	| bedtools bamtobed -i - \
	| bedtools closest -g $genome -sorted \
	-D a -t first -a - -b "$ALU_BED" \
	| awk -v sample="$sample" -v maxd=50000 'BEGIN{OFS="\t"}
	  {
	    dist = $(NF)
	    if (dist < 0) dist = -dist   # abs(distance)
	    if (dist <= maxd) {
	      print sample,
	            $4, $1, $2, $3, $6, $5,
	            $10, $7, $8, $9, $12,
	            $(NF)
	    }
	  }' >> "$OUT"

done