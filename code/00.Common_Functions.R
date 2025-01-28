
#######################################################################################################
## Function to convert top N-gene logic-gated CARs identified by LogiCAR Designer to standard format ##
#######################################################################################################

convert_logic <- function(logic_gate, genes) {
  genes <- strsplit(genes, ", ")[[1]]
  # Map based on the logic gate
  if (logic_gate == "_") { # # Single gene
    return(sprintf("%s", genes[[1]]))
  } else if (logic_gate == "O") { # # Doublets
    return(sprintf("%s | %s", genes[1], genes[2]))
  } else if (logic_gate == "A") {
    return(sprintf("%s & %s", genes[1], genes[2]))
  } else if (logic_gate == "OO") { # # Triplets
    return(sprintf("(%s | %s) | %s", genes[1], genes[2], genes[3]))
  } else if (logic_gate == "OA") {
    return(sprintf("(%s | %s) & %s", genes[1], genes[2], genes[3]))
  } else if (logic_gate == "AO") {
    return(sprintf("(%s & %s) | %s", genes[1], genes[2], genes[3]))
  } else if (logic_gate == "AA") {
    return(sprintf("(%s & %s) & %s", genes[1], genes[2], genes[3]))
  } else if (logic_gate == "AA_A") { # # Quadruples
    return(sprintf("((%s & %s) & %s) & %s", genes[1], genes[2], genes[3], genes[4]))
  } else if (logic_gate == "AO_A") {
    return(sprintf("((%s & %s) | %s) & %s", genes[1], genes[2], genes[3], genes[4]))
  } else if (logic_gate == "OA_A") {
    return(sprintf("((%s | %s) & %s) & %s", genes[1], genes[2], genes[3], genes[4]))
  } else if (logic_gate == "OO_A") {
    return(sprintf("((%s | %s) | %s) & %s", genes[1], genes[2], genes[3], genes[4]))
  } else if (logic_gate == "AA_O") {
    return(sprintf("((%s & %s) & %s) | %s", genes[1], genes[2], genes[3], genes[4]))
  } else if (logic_gate == "AO_O") {
    return(sprintf("((%s & %s) | %s) | %s", genes[1], genes[2], genes[3], genes[4]))
  } else if (logic_gate == "OA_O") {
    return(sprintf("((%s | %s) & %s) | %s", genes[1], genes[2], genes[3], genes[4]))
  } else if (logic_gate == "OO_O") {
    return(sprintf("((%s | %s) | %s) | %s", genes[1], genes[2], genes[3], genes[4]))
  } else if (logic_gate == "A_A_O") { # # Triplets
    return(sprintf("(%s & %s) & (%s | %s)", genes[1], genes[2], genes[3], genes[4]))
  } else if (logic_gate == "A_O_A") {
    return(sprintf("(%s & %s) | (%s & %s)", genes[1], genes[2], genes[3], genes[4]))
  } else if (logic_gate == "O_A_O") {
    return(sprintf("(%s | %s) & (%s | %s)", genes[1], genes[2], genes[3], genes[4]))
  } else if (logic_gate == "O_O_A") {
    return(sprintf("(%s | %s) | (%s & %s)", genes[1], genes[2], genes[3], genes[4]))
  } else if (logic_gate == "AA_A_A") { # # Quintuplets
    return(sprintf("(((%s & %s) & %s) & %s) & %s", genes[1], genes[2], genes[3], genes[4], genes[5]))
  } else if (logic_gate == "AO_A_A") {
    return(sprintf("(((%s & %s) | %s) & %s) & %s", genes[1], genes[2], genes[3], genes[4], genes[5]))
  } else if (logic_gate == "OA_A_A") {
    return(sprintf("(((%s | %s) & %s) & %s) & %s", genes[1], genes[2], genes[3], genes[4], genes[5]))
  } else if (logic_gate == "OO_A_A") {
    return(sprintf("(((%s | %s) | %s) & %s) & %s", genes[1], genes[2], genes[3], genes[4], genes[5]))
  } else if (logic_gate == "AA_O_A") {
    return(sprintf("(((%s & %s) & %s) | %s) & %s", genes[1], genes[2], genes[3], genes[4], genes[5]))
  } else if (logic_gate == "AO_O_A") {
    return(sprintf("(((%s & %s) | %s) | %s) & %s", genes[1], genes[2], genes[3], genes[4], genes[5]))
  } else if (logic_gate == "OA_O_A") {
    return(sprintf("(((%s | %s) & %s) | %s) & %s", genes[1], genes[2], genes[3], genes[4], genes[5]))
  } else if (logic_gate == "OO_O_A") {
    return(sprintf("(((%s | %s) | %s) | %s) & %s", genes[1], genes[2], genes[3], genes[4], genes[5]))
  } else if (logic_gate == "A_A_O_A") {
    return(sprintf("((%s & %s) & (%s | %s)) & %s", genes[1], genes[2], genes[3], genes[4], genes[5]))
  } else if (logic_gate == "A_O_A_A") {
    return(sprintf("((%s & %s) | (%s & %s)) & %s", genes[1], genes[2], genes[3], genes[4], genes[5]))
  } else if (logic_gate == "O_A_O_A") {
    return(sprintf("((%s | %s) & (%s | %s)) & %s", genes[1], genes[2], genes[3], genes[4], genes[5]))
  } else if (logic_gate == "O_O_A_A") {
    return(sprintf("((%s | %s) | (%s & %s)) & %s", genes[1], genes[2], genes[3], genes[4], genes[5]))
  } else if (logic_gate == "AA_A_O") {
    return(sprintf("(((%s & %s) & %s) & %s) | %s", genes[1], genes[2], genes[3], genes[4], genes[5]))
  } else if (logic_gate == "AO_A_O") {
    return(sprintf("(((%s & %s) | %s) & %s) | %s", genes[1], genes[2], genes[3], genes[4], genes[5]))
  } else if (logic_gate == "OA_A_O") {
    return(sprintf("(((%s | %s) & %s) & %s) | %s", genes[1], genes[2], genes[3], genes[4], genes[5]))
  } else if (logic_gate == "OO_A_O") {
    return(sprintf("(((%s | %s) | %s) & %s) | %s", genes[1], genes[2], genes[3], genes[4], genes[5]))
  } else if (logic_gate == "AA_O_O") {
    return(sprintf("(((%s & %s) & %s) | %s) | %s", genes[1], genes[2], genes[3], genes[4], genes[5]))
  } else if (logic_gate == "AO_O_O") {
    return(sprintf("(((%s & %s) | %s) | %s) | %s", genes[1], genes[2], genes[3], genes[4], genes[5]))
  } else if (logic_gate == "OA_O_O") {
    return(sprintf("(((%s | %s) & %s) | %s) | %s", genes[1], genes[2], genes[3], genes[4], genes[5]))
  } else if (logic_gate == "OO_O_O") {
    return(sprintf("(((%s | %s) | %s) | %s) | %s", genes[1], genes[2], genes[3], genes[4], genes[5]))
  } else if (logic_gate == "A_A_O_O") {
    return(sprintf("((%s & %s) & (%s | %s)) | %s", genes[1], genes[2], genes[3], genes[4], genes[5]))
  } else if (logic_gate == "A_O_A_O") {
    return(sprintf("((%s & %s) | (%s & %s)) | %s", genes[1], genes[2], genes[3], genes[4], genes[5]))
  } else if (logic_gate == "O_A_O_O") {
    return(sprintf("((%s | %s) & (%s | %s)) | %s", genes[1], genes[2], genes[3], genes[4], genes[5]))
  } else if (logic_gate == "O_O_A_O") {
    return(sprintf("((%s | %s) | (%s & %s)) | %s", genes[1], genes[2], genes[3], genes[4], genes[5]))
  } else if (logic_gate == "AA_a_A") {
    return(sprintf("((%s & %s) & %s) & (%s & %s)", genes[1], genes[2], genes[3], genes[4], genes[5]))
  } else if (logic_gate == "AO_a_A") {
    return(sprintf("((%s & %s) | %s) & (%s & %s)", genes[1], genes[2], genes[3], genes[4], genes[5]))
  } else if (logic_gate == "OA_a_A") {
    return(sprintf("((%s | %s) & %s) & (%s & %s)", genes[1], genes[2], genes[3], genes[4], genes[5]))
  } else if (logic_gate == "OO_a_A") {
    return(sprintf("((%s | %s) | %s) & (%s & %s)", genes[1], genes[2], genes[3], genes[4], genes[5]))
  } else if (logic_gate == "AA_o_A") {
    return(sprintf("((%s & %s) & %s) | (%s & %s)", genes[1], genes[2], genes[3], genes[4], genes[5]))
  } else if (logic_gate == "AO_o_A") {
    return(sprintf("((%s & %s) | %s) | (%s & %s)", genes[1], genes[2], genes[3], genes[4], genes[5]))
  } else if (logic_gate == "OA_o_A") {
    return(sprintf("((%s | %s) & %s) | (%s & %s)", genes[1], genes[2], genes[3], genes[4], genes[5]))
  } else if (logic_gate == "OO_o_A") {
    return(sprintf("((%s | %s) | %s) | (%s & %s)", genes[1], genes[2], genes[3], genes[4], genes[5]))
  } else if (logic_gate == "AA_a_O") {
    return(sprintf("((%s & %s) & %s) & (%s | %s)", genes[1], genes[2], genes[3], genes[4], genes[5]))
  } else if (logic_gate == "AO_a_O") {
    return(sprintf("((%s & %s) | %s) & (%s | %s)", genes[1], genes[2], genes[3], genes[4], genes[5]))
  } else if (logic_gate == "OA_a_O") {
    return(sprintf("((%s | %s) & %s) & (%s | %s)", genes[1], genes[2], genes[3], genes[4], genes[5]))
  } else if (logic_gate == "OO_a_O") {
    return(sprintf("((%s | %s) | %s) & (%s | %s)", genes[1], genes[2], genes[3], genes[4], genes[5]))
  } else if (logic_gate == "AA_o_O") {
    return(sprintf("((%s & %s) & %s) | (%s | %s)", genes[1], genes[2], genes[3], genes[4], genes[5]))
  } else if (logic_gate == "AO_o_O") {
    return(sprintf("((%s & %s) | %s) | (%s | %s)", genes[1], genes[2], genes[3], genes[4], genes[5]))
  } else if (logic_gate == "OA_o_O") {
    return(sprintf("((%s | %s) & %s) | (%s | %s)", genes[1], genes[2], genes[3], genes[4], genes[5]))
  } else if (logic_gate == "OO_o_O") {
    return(sprintf("((%s | %s) | %s) | (%s | %s)", genes[1], genes[2], genes[3], genes[4], genes[5]))
  } else {
    stop(paste("Unknown LogicGate:", logic_gate))
  }
}


####################################################################################################
## Function (calcEquivalentExpression2) to calculate "equivalent expressions" of logic-gated CARs ##
####################################################################################################

calcEquivalentExpression2 <- function(logic_expr_vec, data_input) {
  flag_dataframe <- if (is.data.frame(data_input)) 1 else 0
  # 1) Define custom binary operators for AND and OR
  `%AND%` <- function(a, b) pmin(a, b)
  `%OR%`  <- function(a, b) pmax(a, b)
  # 2) Create a single environment for vectorized evaluation (parent=baseenv to find '(' etc.)
  env <- new.env(parent = baseenv())
  assign("%AND%", `%AND%`, envir = env)
  assign("%OR%",  `%OR%`,  envir = env)
  # 3) Extract all unique gene names from the logic expressions
  #    Use a regex to match gene names, including those with "_not", and remove spaces
  gene_names <- unique(unlist(regmatches(logic_expr_vec, gregexpr("\\b[A-Za-z0-9_]+\\b", logic_expr_vec))))
  # 4) Assign only the necessary columns of data_df as vectors in env
  for (gene in gene_names) {
    if (gene %in% colnames(data_input)) {
      if (flag_dataframe ==1){
        assign(gene, data_input[[gene]], envir = env) ## for data frame 
      }else{
        assign(gene, data_input[, gene, drop = FALSE], envir = env) ## for matrix
      }
    }
  }
  # Parse expressions ONCE outside the loop
  parsed_exprs <- lapply(logic_expr_vec, function(expr) {
    expr_mod <- gsub("&", "%AND%", expr, fixed = TRUE)
    expr_mod <- gsub("\\|", "%OR%", expr_mod)
    parse(text = expr_mod)[[1]]
  })
  result_list <- vector("list", length(parsed_exprs))
  for (i in seq_along(parsed_exprs)) {
    result_list[[i]] <- eval(parsed_exprs[[i]], envir = env)
  }
  # 7) Return results as a matrix or list
  names(result_list) = logic_expr_vec
  return(result_list)  # Returns a list
}




##############################################
## Function to get all protein coding genes ##
##############################################

library(AnnotationHub)
library(org.Hs.eg.db)
library(biomaRt)

# Function to get protein-coding genes from multiple sources
get_protein_coding_genes <- function() {
  # Install required packages if needed
  if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  required_packages <- c("org.Hs.eg.db", "biomaRt", "AnnotationHub")
  for (pkg in required_packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      BiocManager::install(pkg)
      library(pkg, character.only = TRUE)
    }
  }
  
  # 1. Get genes from org.Hs.eg.db
  genes_orgdb <- AnnotationDbi::select(org.Hs.eg.db,
                                       keys = keys(org.Hs.eg.db, "GENETYPE"),
                                       columns = c("SYMBOL", "GENETYPE"),
                                       keytype = "GENETYPE")
  protein_coding_orgdb <- unique(genes_orgdb$SYMBOL[genes_orgdb$GENETYPE == "protein-coding"])
  
  # 2. Get genes from biomaRt
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  protein_coding_biomart <- getBM(
    attributes = c("hgnc_symbol"),
    filters = c("biotype"),
    values = list(biotype = "protein_coding"),
    mart = ensembl
  )$hgnc_symbol
  
  # 3. Get genes from latest Ensembl via AnnotationHub
  ah <- AnnotationHub()
  # Get latest EnsDb for human
  edb_query <- query(ah, c("EnsDb", "Homo sapiens"))
  latest_edb <- edb_query[[length(edb_query)]]
  
  genes_ensdb <- genes(latest_edb)
  protein_coding_ensdb <- unique(genes_ensdb$symbol[genes_ensdb$gene_biotype == "protein_coding"])
  
  # Combine and compare results
  all_sources <- list(
    org.Hs.eg.db = protein_coding_orgdb,
    biomaRt = protein_coding_biomart,
    EnsDb = protein_coding_ensdb
  )
  
  # Get consensus (genes present in all sources)
  consensus_genes <- Reduce(intersect, all_sources)
  
  # Get union (genes present in any source)
  union_genes <- Reduce(union, all_sources)
  
  # Create summary
  summary <- list(
    consensus_genes = consensus_genes,
    union_genes = union_genes,
    by_source = all_sources,
    counts = list(
      consensus = length(consensus_genes),
      union = length(union_genes),
      by_source = sapply(all_sources, length)
    )
  )
  
  return(summary)
}

# Function to print summary
print_protein_coding_summary <- function(summary) {
  cat("Summary of Human Protein-Coding Genes:\n\n")
  cat("Number of genes by source:\n")
  for (source in names(summary$counts$by_source)) {
    cat(sprintf("%s: %d genes\n", source, summary$counts$by_source[[source]]))
  }
  cat(sprintf("\nConsensus (present in all sources): %d genes\n", summary$counts$consensus))
  cat(sprintf("Union (present in any source): %d genes\n", summary$counts$union))
  
  # Calculate source-specific genes
  cat("\nGenes unique to each source:\n")
  for (source in names(summary$by_source)) {
    unique_genes <- setdiff(summary$by_source[[source]], 
                            unlist(summary$by_source[names(summary$by_source) != source]))
    cat(sprintf("%s: %d genes\n", source, length(unique_genes)))
  }
}


#####################################################################
## Function to find and update gene names using the mygene package ##
#####################################################################

library(mygene)
library(dplyr)

# Function to map gene symbols to the latest HGNC symbols
map_genes_to_hgnc <- function(gene_list) {
  # Query mygene.info to get updated symbols
  result <- queryMany(gene_list, scopes = c("symbol", "alias"), fields = "symbol", species = "human")
  
  # Extract mapping results
  mapping <- data.frame(
    input_gene = result$query,
    updated_gene = result$symbol,
    stringsAsFactors = FALSE
  )
  
  # Handle multiple results for a single query
  mapping <- mapping %>%
    group_by(input_gene) %>%
    summarize(updated_gene = ifelse(any(updated_gene == input_gene), 
                                    input_gene, 
                                    dplyr::first(updated_gene)))
  
  # Handle unmapped genes (NA in updated_gene)
  mapping$updated_gene[is.na(mapping$updated_gene)] <- mapping$input_gene[is.na(mapping$updated_gene)]
  
  return(mapping)
}


