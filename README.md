Identify novel alu from aav-IMR90

1. Remove reads that does not start with ITR using cutadapt.
2. Annotated ALU on IMR90 with RepeatMasker. ALU coordinates recorded in BED format.
3. Mapped reads to IMR90. Filtered for read pairs that are uniquely mapped (mapq>50, primary).
4. R2 were annotated with ALU signature using HMM
5. Intersected R2 from step4 with reads from step3 to get uniquely mapped read pairs that has an ALU signature. 
6. Annotated reads (from step 5) with distance to closest known ALU (from step 2).
