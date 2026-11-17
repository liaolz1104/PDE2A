#!/usr/bin/env Rscript

# Freeze the official human PROGENy and CollecTRI resources used at Checkpoint 10.
# CollecTRI is reconstructed from the OmniPath static table using the exact
# source/target/weight rules implemented by decoupleR 2.16.0 get_collectri().

suppressPackageStartupMessages({
  library(decoupleR)
  library(OmnipathR)
  library(dplyr)
  library(stringr)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 10_freeze_activity_resources.R <project_root>")
root <- normalizePath(args[[1]], mustWork = TRUE)
out <- file.path(root, "results", "10_regulatory_activity", "resources")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

progeny_file <- file.path(out, "PROGENy_human_top500_resource.tsv")
collectri_file <- file.path(out, "CollecTRI_human_resource.tsv")
if (file.exists(progeny_file) || file.exists(collectri_file)) {
  stop("Frozen resource TSV already exists; refusing to overwrite.")
}

progeny_net <- decoupleR::get_progeny(organism = "human", top = 500) %>%
  arrange(source, p_value, target)

collectri_raw <- OmnipathR::static_table(
  query = "interactions",
  resource = "collectri",
  organism = 9606L,
  strict_evidences = FALSE,
  wide = FALSE
)

collectri_interactions <- collectri_raw %>%
  filter(!str_detect(source, "COMPLEX")) %>%
  select(source_genesymbol, target_genesymbol, is_stimulation, is_inhibition)

collectri_complex <- collectri_raw %>%
  filter(str_detect(source, "COMPLEX")) %>%
  select(source_genesymbol, target_genesymbol, is_stimulation, is_inhibition) %>%
  mutate(source_genesymbol = case_when(
    str_detect(source_genesymbol, "JUN") | str_detect(source_genesymbol, "FOS") ~ "AP1",
    str_detect(source_genesymbol, "REL") | str_detect(source_genesymbol, "NFKB") ~ "NFKB",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(source_genesymbol))

collectri_net <- bind_rows(collectri_interactions, collectri_complex) %>%
  distinct(source_genesymbol, target_genesymbol, .keep_all = TRUE) %>%
  mutate(weight = case_when(is_stimulation == 1 ~ 1, is_stimulation == 0 ~ -1)) %>%
  transmute(source = source_genesymbol, target = target_genesymbol, mor = weight) %>%
  filter(!is.na(source), !is.na(target), !is.na(mor)) %>%
  arrange(source, target)

if (nrow(progeny_net) == 0L || nrow(collectri_net) == 0L) stop("A frozen resource is empty.")
if (anyDuplicated(progeny_net[c("source", "target")])) stop("Duplicate PROGENy source-target pairs.")
if (anyDuplicated(collectri_net[c("source", "target")])) stop("Duplicate CollecTRI source-target pairs.")

write_tsv(progeny_net, progeny_file)
write_tsv(collectri_net, collectri_file)
saveRDS(progeny_net, file.path(out, "PROGENy_human_top500_resource.rds"))
saveRDS(collectri_net, file.path(out, "CollecTRI_human_resource.rds"))

progeny_counts <- progeny_net %>% count(source, name = "target_gene_n")
write_tsv(progeny_counts, file.path(out, "PROGENy_pathway_target_counts.tsv"))

metadata <- tibble(
  field = c(
    "freeze_date", "decoupleR_version", "OmnipathR_version",
    "PROGENy_source", "PROGENy_top_targets_per_pathway", "PROGENy_interaction_n", "PROGENy_pathway_n",
    "CollecTRI_source", "CollecTRI_retrieval", "CollecTRI_raw_interaction_n",
    "CollecTRI_processed_interaction_n", "CollecTRI_TF_n", "complex_handling"
  ),
  value = c(
    as.character(Sys.Date()), as.character(packageVersion("decoupleR")), as.character(packageVersion("OmnipathR")),
    "decoupleR::get_progeny via OmniPath PROGENy", "500", nrow(progeny_net), n_distinct(progeny_net$source),
    "OmniPath CollecTRI", "OmnipathR static table fallback; official live endpoint unavailable",
    nrow(collectri_raw), nrow(collectri_net), n_distinct(collectri_net$source),
    "split_complexes=FALSE equivalent: AP1/NFKB complexes collapsed as in decoupleR 2.16.0"
  )
)
write_tsv(metadata, file.path(out, "ACTIVITY_RESOURCE_METADATA.tsv"))

cat("Frozen PROGENy:", nrow(progeny_net), "interactions across", n_distinct(progeny_net$source), "pathways\n")
cat("Frozen CollecTRI:", nrow(collectri_net), "interactions across", n_distinct(collectri_net$source), "TFs\n")
