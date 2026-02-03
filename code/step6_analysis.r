
library(dplyr)

# John's result from HMM signature forALU
r1 <- readxl::read_excel("data/AAV_ALU_HMM_hunt_v1.xlsx")


# alu <- readr::read_tsv("genome/aluRef/IMR90_hap1.Alu.bed", 
#                        col_names = c("chr", "start", "end", "alu", "N", "strand"))
# r2 <- r %>% filter(!sample %in% c("250730_PositiveControl_1_R1_S49_L001", "Undetermined_S0_L001",
#                                  "250730_PositiveControl_2_R2_S50_L001", "250730_PositiveControl_4_R4_S52_L001"))
# r2 <- r2 %>% select(-R2_strand, -alu_chr, -R1_mapq, -alu_start, -alu_end, -alu_strand)
# 
# WriteXLS::WriteXLS(r2, ExcelFileName = "data/imr90_aav_alu_v1.xlsx")


# Alu signature with unique map and distance to known alu
r <- readr::read_tsv("data/r2_alu_res_v1.tsv")

tmp <- r %>% select(-read_id) %>% unique()

# Merge sites
gap_bp <- 3L

sites <- r %>%
  mutate(
    R1_pos = as.integer(R1_start),
    abs_dist = abs(dist_bp)
  ) %>%
  arrange(sample, R1_chr, R1_strand, R1_pos, abs_dist) %>%
  
  # 1) define sites (within 3 bp)
  group_by(sample, R1_chr, R1_strand) %>%
  mutate(
    site_id = cumsum(c(TRUE, diff(R1_pos) > gap_bp))
  ) %>%
  ungroup() %>%
  
  # 2) collapse per site
  group_by(sample, R1_chr, R1_strand, site_id) %>%
  summarise(
    ## site coordinates
    site_start = min(R1_pos),
    site_end   = max(R1_pos),
    
    n_reads    = n(),
    
    dist_bp    = dist_bp[which.min(abs_dist)],
    ## representative read (closest Alu)
    read_id_rep = read_id[which.min(abs_dist)],
    
    ## representative Alu annotation
    alu_type   = alu_type[which.min(abs_dist)],
    alu_chr    = alu_chr[which.min(abs_dist)],
    alu_start  = alu_start[which.min(abs_dist)],
    alu_end    = alu_end[which.min(abs_dist)],
    alu_strand = alu_strand[which.min(abs_dist)],
  
    .groups = "drop"
  ) %>%
  select(-site_id) %>%
  filter(abs(dist_bp) > 1000) 

res <- left_join(sites, r1, by = join_by(read_id_rep == targetName)) %>% select(-read_id_rep) 


WriteXLS::WriteXLS(res, ExcelFileName = "data/AAV_ALU_HMM_hunt_v1_filtered.xlsx")





