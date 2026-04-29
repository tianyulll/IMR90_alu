library(dplyr)
library(openxlsx)
library(ShortRead)

# John's result from HMM signature forALU
r1 <- readxl::read_excel("data/AAV_ALU_HMM_hunt_v1.xlsx")


# Alu signature with unique map and distance to known alu
r <- readr::read_tsv("data/r2_alu_res_v1.tsv")

# Merge sites
# Filter based on distance, orientation and signal strength
gap_bp <- 5L

sites <- r %>%
  mutate(R1_start_i = as.integer(R1_start), abs_dist = abs(dist_bp)) %>%
  arrange(sample, R1_chr, R1_strand, R1_start_i, abs_dist, desc(R1_mapq)) %>%
  group_by(sample, R1_chr, R1_strand) %>%
  mutate(site_id = cumsum(c(TRUE, diff(R1_start_i) > gap_bp))) %>%
  group_by(sample, R1_chr, R1_strand, site_id) %>%
  summarise(
    pick = 1L,  # already sorted: closest Alu, then highest MAPQ
    read_id   = read_id[pick],
    R1_start  = R1_start[pick],
    R1_end    = R1_end[pick],
    alu_type   = alu_type[pick],
    alu_chr    = alu_chr[pick],
    alu_start  = alu_start[pick],
    alu_end    = alu_end[pick],
    alu_strand = alu_strand[pick],
    dist_bp    = dist_bp[pick],
    R1_mapq   = R1_mapq[pick],
    .groups = "drop"
  ) %>%
  select(-pick, -site_id)

res <- left_join(sites, r1, by = join_by(read_id == targetName)) %>%
  filter(abs(dist_bp) > 2000, fullEval < 1) %>%
  arrange(abs(dist_bp))


# Add R1 seqs back to the df
add_fastq_seq_stream <- function(
    df,
    id_col = "read_id",
    pattern = "fastq/*R1_001.fastq.gz",
    chunk_size = 1e6
) {
  fastq_files <- Sys.glob(pattern)
  
  ids_wanted <- unique(df[[id_col]])
  ids_wanted <- sub("/1$", "", sub("\\s.*$", "", ids_wanted))
  
  found <- list()
  k <- 1
  
  for (fq in fastq_files) {
    message("Reading: ", fq)
    
    streamer <- FastqStreamer(fq, n = chunk_size)
    
    repeat {
      fq_chunk <- yield(streamer)
      if (length(fq_chunk) == 0) break
      
      fq_ids <- sub("/1$", "", sub("\\s.*$", "", as.character(id(fq_chunk))))
      hit <- fq_ids %in% ids_wanted
      
      if (any(hit)) {
        found[[k]] <- tibble(
          read_id = fq_ids[hit],
          R1_seq = as.character(sread(fq_chunk))[hit]
        )
        k <- k + 1
      }
    }
    
    close(streamer)
  }
  
  seq_df <- bind_rows(found) %>%
    distinct(read_id, .keep_all = TRUE)
  
  df %>%
    mutate(!!id_col := sub("/1$", "", sub("\\s.*$", "", .data[[id_col]]))) %>%
    left_join(seq_df, by = setNames("read_id", id_col))
}

df2 <- add_fastq_seq_stream(res)


# Add link to UCSC browser session with IMR90 tracks
export_ucsc_excel <- function(
    df,
    file = "df_ucsc_links.xlsx",
    hub_url = "https://microb191.med.upenn.edu/tracks/genome/IMR90Track/hub.txt",
    genome = "IMR90"
) {
  
  links <- paste0(
    "https://genome.ucsc.edu/cgi-bin/hgTracks?",
    "hubUrl=", URLencode(hub_url, reserved = TRUE),
    "&genome=", genome,
    "&position=", df$R1_chr, ":", df$R1_start, "-", df$R1_end,
    "&highlight=", df$R1_chr, ":", df$R1_start, "-", df$R1_end
  )
  wb <- createWorkbook()
  addWorksheet(wb, "Sheet1")
  writeData(wb, "Sheet1", df)
  
  for (i in seq_len(nrow(df))) {
    writeFormula(
      wb, "Sheet1",
      x = sprintf(
        '=HYPERLINK("%s","%s")',
        links[i],
        as.character(df[[1]][i])
      ),
      startCol = 1,
      startRow = i + 1
    )
  }
  
  saveWorkbook(wb, file, overwrite = TRUE)
}

export_ucsc_excel(df = df2, file = "data/IMR90_Alu_sites_v1.xlsx")


  

# alu <- readr::read_tsv("genome/aluRef/IMR90_hap1.Alu.bed", 
#                        col_names = c("chr", "start", "end", "alu", "N", "strand"))
# r2 <- r %>% filter(!sample %in% c("250730_PositiveControl_1_R1_S49_L001", "Undetermined_S0_L001",
#                                  "250730_PositiveControl_2_R2_S50_L001", "250730_PositiveControl_4_R4_S52_L001"))
# r2 <- r2 %>% select(-R2_strand, -alu_chr, -R1_mapq, -alu_start, -alu_end, -alu_strand)
# 
# WriteXLS::WriteXLS(r2, ExcelFileName = "data/imr90_aav_alu_v1.xlsx")




