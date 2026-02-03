# Identify novel alu from aav-IMR90
We want to find integrations in unannotated ALU from AAV experiments. In this trial we used IMR90 cells. \

1. Remove reads that does not start with ITR using cutadapt. This step filters out irrelevant reads.
2. Annotated ALU on IMR90 with RepeatMasker. ALU coordinates recorded in BED format. IMR90 does not have a known annotation so so ALU was annotated from scratch.
3. Mapped reads to IMR90. Filtered for read pairs that are uniquely mapped (mapq>50, primary). This step filters out reads with low quality and unaligned reads.
4. R2 were annotated with ALU signature using HMM. This step looks for R2 that are likely to be in an ALU.
10. Intersected annotated R2 (from step 4) with filtered reads (from step 3) to get uniquely mapped read pairs that has an ALU signature. 
11. Annotated reads (from step 5) with distance to closest known ALU (from step 2). R1 was used to infer intergation site. Reads were then collapsed and summarized.
