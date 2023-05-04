library(dplyr)
library(tidyr)
library(tibble)
library(speedyseq)


if(!interactive()){
  args <- commandArgs(trailingOnly = TRUE) 
  if (length(args) != 4) {
    stop(paste("USAGE: Rscript make_project_phylo.R path_to_manifest_file.txt project_name <metaphlan|kraken|humann> output.rds\n", args))
  }
  if (!file.exists(args[1])) {
    stop(paste("file ", args[1], " does not exist!"))
  }
  outfile = args[4]
  if (file.exists(outfile)) {
    stop(paste("output file ", outfile, " already exists!"))
  }
  
  file_paths <- readLines(args[1])
  proj_name <- args[2]
  print("paths to aggregate")
  print(file_paths)
  res_type = tolower(args[3])
 } else {
  proj_name = "test"

    file_paths = c(
    "~/GitHub/Shouval-CART-2023/pipeline_results_1m/brinkvd/isabl/data_lake/analyses/38/80/3880/metaphlan/08788_1511N_HW7L2BBXX_metaphlan3_profile.txt",
    "~/GitHub/Shouval-CART-2023/pipeline_results_1m/brinkvd/isabl/data_lake/analyses/40/10/4010/metaphlan/13481_2512B_HWNLKDSX3_metaphlan3_profile.txt"
  )
  res_type = "metaphlan"
  outfile = "tmp_metaphlan.rds"
  
  file_paths = c(
    "~/GitHub/Shouval-CART-2023/pipeline_results_1m/brinkvd/isabl/data_lake/analyses/38/80/3880/humann/08788_1511N_HW7L2BBXX_humann3_pathabundance_cpm.tsv",
    "~/GitHub/Shouval-CART-2023/pipeline_results_1m/brinkvd/isabl/data_lake/analyses/40/10/4010/humann/13481_2512B_HWNLKDSX3_humann3_pathabundance_cpm.tsv"
  )
  res_type = "humann"
  outfile = "tmp_hummann.rds"

  file_paths = c(
    "/Users/nwaters/GitHub/Shouval-CART-2023/pipeline_results_full/brinkvd/isabl/data_lake/analyses/28/93/2893/kraken2/08788_1715A_HW7L2BBXX_kraken2.bracken.S.report",
    "/Users/nwaters/GitHub/Shouval-CART-2023/pipeline_results_full/brinkvd/isabl/data_lake/analyses/28/90/2890/kraken2/10194_2108A_H7WG3BBXY_kraken2.bracken.S.report"
  )
  res_type = "braken"
  outfile = "tmp_hummann.rds"
  
}

if (!res_type %in% c("metaphlan", "humann", "bracken")) {
  stop(paste("can only aggregate metaphlan, humann, or braken results, not ", res_type))
}



aggregate_metaphlan <- function(file_paths, outfile){
  print("  parsing input files")
  bigdf_at_species <- lapply(file_paths, function(x){
    read.csv(x, col.names = c("clade_name", "taxids","perc", "perc_at", "reads_from_clade")) %>% 
      mutate(experiment_id = gsub("_metaphlan3_profile.txt", "", basename(x))) %>% 
      filter(grepl("\\|s__", clade_name) & !grepl("\\|t__", clade_name)) 
  }) %>% bind_rows() %>% 
    select(-perc, -perc_at) %>% 
    pivot_wider(names_from = experiment_id, values_from = reads_from_clade, values_fill = 0)

  

  print("  building taxonomy table")
  # we will populate the actual sample_data later, for now make an empty df entry
  tax_tab <- build_tax_table(
    lapply(gsub("\\|", ";", bigdf_at_species$clade_name), phyloseq::parse_taxonomy_qiime))  %>% as_tibble() %>% 
    mutate(clade_name = bigdf_at_species$clade_name, 
           taxids = bigdf_at_species$taxids) %>% column_to_rownames("clade_name") %>%
    mutate(Species_taxid = gsub(".*?\\|.*?\\|.*?\\|.*?\\|.*?\\|.*?\\|(.*)", "\\1", taxids)) %>% as.matrix() %>% 
    tax_table()
  
  print("  creating phyloseq object")
  raw_mpa_phy <- phyloseq(
    otu_table(bigdf_at_species %>%
                select(clade_name, 3:ncol(.)), taxa_are_rows = TRUE),
    tax_tab,
    sample_data(data.frame(dummymetadata="dummy", fid=colnames(bigdf_at_species)[3:ncol(bigdf_at_species)]) %>% 
                  column_to_rownames("fid"))
  ) 
  print("  saving results")
  
  saveRDS(raw_mpa_phy, file = outfile)
}

aggregate_humann <- function(file_paths, outfile){
  print("  parsing input files")
  
  bigdf <- lapply(file_paths, function(x){
    read.csv(x, col.names = c("funct_id_raw", "count"), sep="\t") %>% 
      mutate(experiment_id = gsub("_humann*$", "", basename(x))) #%>% filter(grepl("\\|g__",funct_id_raw ))
  }) %>% 
    bind_rows() %>% 
    pivot_wider(names_from = experiment_id, values_from = count, values_fill = 0)
  
  
  print("  building functional 'taxonomy' table")
  funct_taxonomy <- data.frame(funct_id_raw = bigdf$funct_id_raw) %>% 
    mutate(
      funct_id_clean = ifelse(!grepl("\\|", funct_id_raw), paste(funct_id_raw, "|all"), funct_id_raw),
      funct_id = gsub("(.*)\\|.*", "\\1", funct_id_clean),
      sub_funct = gsub("(.*)\\|(.*)", "\\2", funct_id_clean)
    ) %>% 
    column_to_rownames("funct_id_raw") %>% 
    select(-funct_id_clean) %>% 
    as.matrix()
  
  print("  creating phyloseq object")
  raw_humann_phy <- phyloseq(
    otu_table(bigdf, taxa_are_rows = TRUE),
    tax_table(funct_taxonomy),
    sample_data(data.frame(dummymetadata="dummy", fid=colnames(bigdf)[2:ncol(bigdf)]) %>% 
                  column_to_rownames("fid"))
  ) 
  print("  saving results")
  saveRDS(raw_humann_phy, file = outfile)
}

aggregate_braken <- function(file_paths, outfile){
  print("  parsing input files")
  
  bigbrakendf <- lapply(file_paths, function(x){
    read.csv(x, col.names = c("perc", "reads", "reads_at_species", "taxlev", "taxid", "taxon"), sep="\t", header = FALSE) %>%
      mutate(
        taxon = trimws(taxon),
        orig_tax_level = taxlev,
        experiment_id = gsub("_kraken2.bracken.S.report$", "", basename(x)),
      ) %>%
      filter(taxlev %in% c("D", "P",  "C", "O", "F", "G", "S")) %>%
      pivot_wider(names_from = taxlev, values_from = taxon) %>% 
      fill(S, .direction = "up") %>%
      fill(D:G, .direction = "down") %>% 
      filter(orig_tax_level == "S") %>% 
      rename(K=D) %>% 
      mutate(clade_name=paste0("k__", K, "|p__", P, "|c__", C, "|o__", O, "|f__", F, "|g__", G, "|s__", S)) %>%
      select(clade_name, K:S,  taxid, reads_at_species, experiment_id)
  }) %>% 
    bind_rows() %>% 
    pivot_wider(names_from = experiment_id, values_from = reads_at_species, values_fill = 0)
  
  
  print("  building taxonomy table")
  brack_tax <- bigbrakendf %>% select(clade_name,  K:taxid) %>% 
    rename("Kingdom"=K, "Phylum"=P, "Class"=C,  "Order"=O, "Family"=F, "Genus"=G, "Species"="S") %>% 
    column_to_rownames("clade_name") %>% 
    as.matrix()
  otu_df <- bigbrakendf %>% select(-c(K:taxid))
  print("  creating phyloseq object")
  raw_brack_phy <- phyloseq(
    otu_table(otu_df %>% column_to_rownames("clade_name"), taxa_are_rows = TRUE),
    tax_table(brack_tax ),
    sample_data(data.frame(dummymetadata="dummy", fid=colnames(otu_df)[2:ncol(otu_df)]) %>% 
                  column_to_rownames("fid"))
  ) 
  print("  saving results")
  saveRDS(raw_brack_phy, file = outfile)
}
print(paste("Aggregating", res_type, "results"))
if (res_type =="metaphlan"){
  aggregate_metaphlan(file_paths, outfile)
} else if (res_type == "humann"){
  aggregate_humann(file_paths, outfile)
} else if (res_type == "bracken"){
  aggregate_braken(file_paths, outfile)
} else {
  stop("unsupported results type!")
}
print("Done. Exiting")


