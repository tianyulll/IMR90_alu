
# Transfer bcl files
rsync -av --no-inc-recursive --whole-file tianyu@microb120.med.upenn.edu:/media/sequencing/Illumina/250730_M03249_0041_000000000-LT9P4/ data/
  
# Convert bcl 2 fastq
bcl-convert \
  --bcl-input-directory data \
  --sample-sheet 250730_samplesheet_bcl.csv  \
  --output-directory fastq
