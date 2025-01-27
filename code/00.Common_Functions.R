
#######################################################################################################
## Function to convert top N-gene logic-gated CARs identified by LogiCAR Designer to standard format ##
#######################################################################################################

convert_logic <- function(logic_gate, genes) {
  genes <- strsplit(genes, ", ")[[1]]
  # Map based on the logic gate
  if (logic_gate == "_") {
    # Single gene
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
  } else if (logic_gate == "AA_O") { # # Triplets
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
  # return(do.call(cbind, result_list))  # Returns a matrix
  names(result_list) = logic_expr_vec
  return(result_list)  # Returns a list
}


