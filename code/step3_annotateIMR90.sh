# Alu annotation
RepeatMasker \
  -pa 4 \
  -xsmall \
  -gff \
  IMR90_hap1.fasta.gz

# extract alu coordinates
awk 'BEGIN{OFS="\t"}
     $11=="SINE/Alu" {
       strand = ($9=="C" ? "-" : "+")
       print $5, $6-1, $7, $10, ".", strand
     }' IMR90_hap2.fasta.out > aluRef/IMR90_hap2.Alu.bed


awk 'BEGIN{OFS="\t"}
     $11=="SINE/Alu" {
       strand = ($9=="C" ? "-" : "+")
       print $5, $6-1, $7, $10, ".", strand
     }' IMR90_hap1.fasta.out > aluRef/IMR90_hap1.Alu.bed