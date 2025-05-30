here::i_am("GraftM-analysis/update_annotations.R")
library(here)
library(tidyverse)

dram_annotations_directory <- here("data", "DRAM_annotations_v4")
output_annotations_directory <- here("data", "DRAM_annotations_v4")

read_dram_annotations <- function(filename = "all_annotations.tsv") {
  annotations <- read_tsv(here(dram_annotations_directory, filename)) %>%
    rename(gene_id = 1)
}

dram_annotations <- read_dram_annotations(filename = "all_annotations_raw.tsv")

graftm_assignments <- tribble(
  ~graftm_call, ~gene_name,
  # dsrA
  "Root; dsrA; Acidobacteriota", "dsrA",
  "Root; dsrA; Acidobacteriota_2", "dsrA",
  "Root; dsrA; Acidobacteriota_3", "dsrA",
  "Root; dsrA; Bacteroidota", "dsrA",
  "Root; dsrA; Chloroflexota", "dsrA",
  "Root; dsrA; Desulfovibrionales", "dsrA",
  "Root; dsrA; Myxococcota", "dsrA",
  "Root; dsrA; Nitrospirota", "dsrA",
  "Root; dsrA; Syntrophia", "dsrA",
  "Root; dsrA; Syntrophia_2", "dsrA",
  "Root; dsrA; Syntrophobacteria", "dsrA",
  "Root; dsrA; Syntrophorhabdia", "dsrA",
  "Root; dsrA; Verrucomicrobiota", "dsrA",
  "Root; rdsrA; Alphaproteobacteria", "rdsrA",
  "Root; rdsrA; Gammaproteobacteria", "rdsrA",
  "Root; rdsrA; Myxococcota_2", "rdsrA",
  # dsrB
  "Root; dsrB", "dsrB",
  "Root; dsrB; Acidobacteriota", "dsrB",
  "Root; dsrB; Acidobacteriota_2", "dsrB",
  "Root; dsrB; Bacteroidota", "dsrB",
  "Root; dsrB; Chloroflexota", "dsrB",
  "Root; dsrB; Desulfovibrionales", "dsrB",
  "Root; dsrB; Myxococcota", "dsrB",
  "Root; dsrB; Nitrospirota", "dsrB",
  "Root; dsrB; Syntrophia", "dsrB",
  "Root; dsrB; Syntrophia_2", "dsrB",
  "Root; dsrB; Syntrophobacteria", "dsrB",
  "Root; dsrB; Syntrophorhabdia", "dsrB",
  "Root; dsrB; Verrucomicrobiota", "dsrB",
  "Root; rdsrB", "rdsrB",
  "Root; rdsrB; Alphaproteobacteria", "rdsrB",
  "Root; rdsrB; Gammaproteobacteria_2", "rdsrB",
  # narG_nxrA
  "Root; NarG_Proteobacteria", "narG",
  "Root; NarG_Actinobacteriota", "narG",
  "Root; NarG_Acidobacteriota", "narG",
  "Root; NarG_Cyanobacteria", "narG",
  "Root; NarG_Mixed", "narG",
  "Root; periplasmic_NxrA", "nxrA",
  "Root; cytoplasmic_NxrA_Mixed", "nxrA",
  # narH_nxrB
  "Root; NarH_Acidobacteriota", "narH",
  "Root; NarH_Actinobacteriota", "narH",
  "Root; NarH_Desulfobacterota", "narH",
  "Root; NarH_Proteobacteria", "narH",
  "Root; NxrB", "nxrB",
  "Root; NxrB_2", "nxrB",
  # amoA_pmoA
  "Root; amoA_Nitrospirae", "amoA",
  "Root; amoA_Mixed", "amoA",
  "Root; pmoA_Mixed", "pmoA",
  "Root; pmoA_Mixed_2", "pmoA",
  "Root; pmoA_Mixed_3", "pmoA",
  "Root; pmoA_Rhizobiales", "pmoA",
  # amoB_pmoB
  "Root; amoB_Nitrospirae", "amoB",
  "Root; amoB_Mixed", "amoB",
  "Root; pmoB_Mixed", "pmoB",
  "Root; pmoB_Methylococcales", "pmoB",
  # amoC_pmoC
  "Root; amoC_Nitrospirae", "amoC",
  "Root; amoC_Mixed", "amoC",
  "Root; pmoC_Mixed", "pmoC",
  "Root; pmoC_Methylococcales", "pmoC",
  # norZ_norB
  "Root", "norB",
  "Root; cNor_Alphaproteobacteria", "norB",
  "Root; cNor_Bacteroidota", "norB",
  "Root; cNor_Gammaproteobacteria", "norB",
  "Root; cNor_Verrucomicrobiota", "norB",
  "Root; cNor_Mixed", "norB",
  "Root; qNor_Alphaproteobacteria", "norB",
  "Root; qNor_Mixed", "norB",
  "Root; qNor_Mixed_2", "norB",
  "Root; qNor_Planctomycetota", "norB"
)

new_ko_ids <- tribble(
  ~gene_name, ~new_ko_id,
  # dsrA
  "dsrA", "K111800",
  "rdsrA", "K111801",
  # dsrB
  "dsrB", "K111810",
  "rdsrB", "K111811",
  # narG_nxrA
  "narG", "K003700",
  "nxrA", "K003701",
  # narH_nxrB
  "narH", "K003710",
  "nxrB", "K003711",
  # amoA_pmoA
  "pmoA", "K109440",
  "amoA", "K109441",
  # amoB_pmoB
  "pmoB", "K109450",
  "amoB", "K109451",
  # amoC_pmoC
  "pmoC", "K109460",
  "amoC", "K109461",
  # norZ_norB
  "norB", "K045610",
  "norZ", "K045611"
)

graftm_inputs <- tribble(
  ~gene, ~ko_id,
  "dsrA", "K11180",
  "dsrB", "K11181",
  "narG_nxrA", "K00370",
  "narH_nxrB", "K00371",
  "amoA_pmoA", "K10944",
  "amoB_pmoB", "K10945",
  "amoC_pmoC", "K10946",
  "norZ_norB", "K04561"
  ) %>%
  mutate(
    data = map(gene, ~ read_tsv(here("GraftM-analysis", ., "emerge_all_read_tax.tsv"), col_names = c("gene_id", "graftm_call")))
  ) %>%
  unnest(data) %>%
  left_join(graftm_assignments) %>%
  left_join(new_ko_ids)

# Annotations ko_id: change to KO from graftm_inputs new_ko_id
new_annotations <- dram_annotations %>%
  left_join(graftm_inputs %>% select(gene_id, new_ko_id)) %>%
  mutate(
    ko_id = map2_chr(ko_id, new_ko_id, ~ ifelse(is.na(.y), .x, .y))
  )

write_tsv(
  new_annotations %>% select(gene_id, ko_id),
  here(output_annotations_directory, "updated_ko_ids.tsv"),
  na = ""
  )

# Replace annotations ko_id column with new column
# DRAM_DIR=data/DRAM_annotations_v4
# awk 'BEGIN{OFS=FS="\t"}FNR==NR{a[NR]=$2;next}{$10=a[FNR]}1' $DRAM_DIR/updated_ko_ids.tsv $DRAM_DIR/all_annotations_raw.tsv > $DRAM_DIR/all_annotations.tsv
