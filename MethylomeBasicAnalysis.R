library(Biostrings)
library(ggplot2)
library(dplyr)
library(ggbreak)

#function to find GANTC motifs in genome (from a fasta sequence) and return list
find_gantc_motifs <- function(genome_file) 
{
  genome <- readDNAStringSet(genome_file)
  genome_seq <- as.character(genome[[1]])

  #change to GATC for E. coli, or to motif of interest for the organism of interest 
  motif_pattern <- "GA.TC"
  
  # Find positions on forward strand
  forward_matches <- gregexpr(motif_pattern, genome_seq)[[1]]
  # Position of A (relative to start of GANTC)
  forward_positions <- forward_matches + 1  
  
  motifs_list <- list()
  
  for (i in seq_along(forward_positions)) 
  {
    if (forward_positions[i] > 0) 
    {
      pos_fwd <- forward_positions[i]
      # Reverse strand A is 2 positions ahead for GANTC (change for different motifs)
      pos_rev <- pos_fwd + 2  
      motifs_list[[i]] <- data.frame(MotifID = i, Position = c(pos_fwd, pos_rev), MotifStrand = c("+", "-"))
    }
  }
  
  gantc_motifs <- do.call(rbind, motifs_list)
  return(gantc_motifs)
}

# function to filter modkit data to only include GANTC motif positions
filter_to_gantc <- function(data, gantc_motifs) 
{
  data_filtered <- data %>%
    inner_join(gantc_motifs, by = c("Position" = "Position", "Strand" = "MotifStrand"))
  return(data_filtered)
}

#function to pool data from + and - negative strands if desired. each unique location is defined using MotifID
#using coverage weighted average
pool_strands <- function(data) 
{
  pooled <- data %>%
    group_by(MotifID) %>%
    summarise(
      Frm = sum(Frm * Coverage, na.rm = TRUE) / sum(Coverage, na.rm = TRUE),
      Coverage = sum(Coverage, na.rm = TRUE), 
      n_strands = n(),
      Position = min(Position)  # Forward strand position
    ) %>%
    ungroup()
  
  return(pooled)
}

#set working directory to location of modkit output
setwd("~/Desktop/Methylome/")

#input files
sample1 <- read.table("sample1.tsv", header = TRUE, sep = "\t")
sample2 <- read.table("sample2.tsv", header = TRUE, sep = "\t")

#Find GANTC motifs and filter
genome_file <- "fasta/GANTC/Methylobacterium_PA1.fasta"  #genome sequence of reference organism
gantc_motifs <- find_gantc_motifs(genome_file)           #making a list of all motifs of interest from that reference

cat("Found", length(unique(gantc_motifs$MotifID)), "GANTC motifs\n")
cat("Total positions (2 per motif):", nrow(gantc_motifs), "\n")

#Simplifying modkit output to a few columns of interest
sample1 <- sample1[, -c(1, 7, 8, 9, 10, 14, 15, 17, 18)]
colnames(sample1) = c("Start","Position","Mod","Coverage","Strand","Frm","Num mod","Num og","Fail")
sample2 <- sample2[, -c(1, 7, 8, 9, 10, 14, 15, 17, 18)]
colnames(sample2) = c("Start","Position","Mod","Coverage","Strand","Frm","Num mod","Num og","Fail")

#removing data not corresponding to the GANTC / whatever motif is chosen
sample1_filtered <- filter_to_gantc(sample1, gantc_motifs)
sample2_filtered <- filter_to_gantc(sample2, gantc_motifs)
cat("Sample 1: Filtered to", nrow(sample1_filtered), "rows (strand-matched)\n")
cat("Sample 2: Filtered to", nrow(sample2_filtered), "rows (strand-matched)\n")

#Pool data across strands (optional) - we will always save strand specific data anyway
sample1_pooled <- pool_strands(sample1_filtered)
sample2_pooled <- pool_strands(sample2_filtered)
cat("Sample 1: Pooled to", nrow(sample1_pooled), "positions\n")
cat("Sample 2: Pooled to", nrow(sample2_pooled), "positions\n")

sample1_pooled$Sample <- "Sample1"
sample2_pooled$Sample <- "Sample2"

# Combine both samples - in case of replicates - just one sample if not
combined_data <- rbind(sample1_filtered, sample2_filtered)

# Calculate strand-specific Frm for each MotifID
# Average across both samples for each strand separately
strand_specific_df <- combined_data %>%
  group_by(MotifID, Strand) %>%
  summarise(
    # Coverage-weighted average Frm for this strand across both samples
    Frm_strand = sum(Frm * Coverage, na.rm = TRUE) / sum(Coverage, na.rm = TRUE),
    Coverage_strand = sum(Coverage, na.rm = TRUE),
    n_obs = n(),
    .groups = "drop"
  ) %>%
  # Reshape to wide format with Frm+ and Frm- columns
  tidyr::pivot_wider(
    id_cols = MotifID,
    names_from = Strand,
    values_from = c(Frm_strand, Coverage_strand),
    names_glue = "{.value}_{Strand}"
  ) %>%
  # Rename columns to FrmP and FrmN
  rename(
    FrmP = `Frm_strand_+`,
    FrmN = `Frm_strand_-`,
    CoverageP = `Coverage_strand_+`,
    CoverageN = `Coverage_strand_-`
  )

#add the overall coverage-weighted Frm (pooled across strands)
overall_frm <- combined_data %>%
  group_by(MotifID) %>%
  summarise(
    Frm = sum(Frm * Coverage, na.rm = TRUE) / sum(Coverage, na.rm = TRUE),
    Total_Coverage = sum(Coverage, na.rm = TRUE),
    Position = min(Position)
  )

#Merge strand-specific with overall dataframe
final_df <- strand_specific_df %>%
  left_join(overall_frm, by = "MotifID") %>%
  select(MotifID, Position, Frm, FrmP, FrmN, Total_Coverage, CoverageP, CoverageN)

#View results
head(final_df)

#save as a csv file
write.csv(final_df,"csv_files/MethanolRep2.csv",row.names = FALSE)
