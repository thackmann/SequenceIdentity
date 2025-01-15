# Calculate Relative Evolutionary Divergence of a Set of Tree Tips
#
# This script calculates relative evolutionary divergence (RED) over
# a set of tips in a phylogenetic tree.  The calculation is done following
# Parks et al. (2018. Nat Biotechnol. 36:996-1004).
# For a pair of tips, it finds the most recent common ancestor (MCRA) and gets
# the value of RED with the castor::get_reds() function.  It can also use
# RED values from the phylogenetic tree, if provided.  It can do the calculation
# over a set of tips, and it can do it in parallel by breaking the set into chunks.
#
# Author:  Timothy Hackmann
# Date: 15 January 2025

#=== Define functions ====
  #' Install Missing CRAN Packages
  #'
  #' This function checks for missing CRAN packages and installs them if they are not already installed.
  #'
  #' @param packages A character vector of CRAN package names to check and install if missing.
  #' @return None. The function installs missing packages and provides a message if installation occurs.
  #' @examples
  #' cran_packages <- c("dplyr", "ggplot2")
  #' install_missing_cran_packages(cran_packages)
  install_missing_cran_packages <- function(packages) {
    # Identify missing packages
    missing_cran <- packages[!(packages %in% installed.packages()[, "Package"])]
    
    # Check if there are any missing packages
    if (length(missing_cran) > 0) {
      message("Installing missing CRAN packages: ", paste(missing_cran, collapse = ", "))
      install.packages(missing_cran)
    } else {
      message("All packages are already installed.")
    }
  }
  
  #' Install Missing Bioconductor Packages
  #'
  #' This function checks for missing Bioconductor packages and installs them using BiocManager if they are not already installed.
  #'
  #' @param packages A character vector of Bioconductor package names to check and install if missing.
  #' @return None. The function installs missing packages and provides a message if installation occurs.
  #' @examples
  #' bioc_packages <- c("Biostrings", "GenomicRanges")
  #' install_missing_bioc_packages(bioc_packages)
  install_missing_bioc_packages <- function(packages) {
    # Ensure BiocManager is installed
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    
    # Identify missing packages
    missing_bioc <- packages[!(packages %in% installed.packages()[, "Package"])]
    
    # Check if there are any missing packages
    if (length(missing_bioc) > 0) {
      message("Installing missing Bioconductor packages: ", paste(missing_bioc, collapse = ", "))
      BiocManager::install(missing_bioc)
    } else {
      message("All Bioconductor packages are already installed.")
    }
  }
  
  #' Get Most Recent Common Ancestor (MRCA) for a Set of Tips
  #'
  #' Identifies the most recent common ancestor (MRCA) for a given set of tips in a phylogenetic tree.
  #' This function is based on \code{ape::getMRCA} but uses a C++ implementation for faster execution.
  #'
  #' @param phy A phylogenetic tree of class \code{"phylo"}.
  #' @param tip A character vector of tip labels or an integer vector of tip indices for which the MRCA is to be determined.
  #'
  #' @return An integer representing the node index of the MRCA if the tips are valid and \code{phy} is of class \code{"phylo"}.
  #' If less than two tips are provided, the function returns \code{NULL}.
  #'
  #' @examples
  #' library(ape)
  #' tree <- read.tree(text = "(((A:1,B:1):1,(C:1,D:2):2):2,E:3);")
  #' get_mrca(tree, c("A", "B"))  # Should return the node index for the MRCA of tips A and B
  #'
  #' @export
  get_mrca <- function(phy, tip) {
    if (!inherits(phy, "phylo")) 
      stop("object \"phy\" is not of class \"phylo\"")
    if (length(tip) < 2) 
      return(NULL)
    
    Ntip <- length(phy$tip.label)
    tip_indices <- if (is.character(tip)) match(tip, phy$tip.label) else tip
    
    # Call the C++ function
    mrca <- get_mrca_cpp(as.matrix(phy$edge), Ntip, tip_indices)
    return(mrca)
  }

  #' Get RED of the Common Ancestor of Pair of Tips
  #'
  #' This function calculates or extracts the Relative Evolutionary Divergence (RED)
  #' of the most recent common ancestor (MRCA) of a given pair of tips in a phylogenetic tree.
  #'
  #' @param tree A phylogenetic tree of class "phylo".
  #' @param tip1 A character string representing the first tip.
  #' @param tip2 A character string representing the second tip.
  #' @param reds Optional. A numeric vector of REDs for nodes in the tree.
  #'             If not provided, REDs will be calculated using the `castor` package.
  #' @param use_cpp Logical indicating whether to use the C++ implementation for faster computation.
  #'
  #' @return A numeric value representing the RED of the most recent common ancestor (MRCA) of the tips.
  #'
  #' @examples
  #' tree <- ape::read.tree(text = "(((A:1,B:1):1,(C:1,D:2):2):2,E:3);")
  #' reds <- castor::get_reds(tree)
  #' calculate_red_pair(tree, "A", "B", reds)
  calculate_red_pair <- function(tree, tip1, tip2, reds = NULL, use_cpp=TRUE) {
    if (!inherits(tree, "phylo")) 
      stop("object 'tree' is not of class 'phylo'")
    
    # Get RED values if not provided
    if (is.null(reds)) {
      reds <- castor::get_reds(tree)
    }
    
    if(use_cpp==TRUE)
    {
      # Get tip indices
      tip_indices <- match(c(tip1, tip2), tree$tip.label)
      if (any(is.na(tip_indices))) stop("One or both tips not found in tree.")
      
      # Get RED value for tips
      red_value <- calculate_red_pair_cpp(
                      edge = as.matrix(tree$edge), 
                      reds = reds, 
                      Ntip = length(tree$tip.label), 
                      tip_indices = tip_indices - 1  # Adjust for 0-based indexing
                      )
    }else if(use_cpp==FALSE)
    { 
         # Get the MRCA of the two tips
         mrca <- ape::getMRCA(tree, c(tip1, tip2))
         
         # Get RED value for MRCA of tips
         red_value <- reds[mrca - length(tree$tip.label)] # Adjust for tip indices}
    }
    
    return(red_value)
  }
  
  #' Calculate Pairwise RED for a Set of Tips
  #'
  #' Calculates RED values for all pairs of tips in a phylogenetic tree.
  #'
  #' @param tree A phylogenetic tree of class "phylo".
  #' @param tips A character vector of tip labels to calculate pairwise RED values for.
  #' @param reds Optional. A numeric vector of REDs for nodes in the tree.
  #'             If not provided, REDs will be calculated using the `castor` package.
  #'
  #' @return A dataframe containing:
  #'   - `Tip1`: First tip.
  #'   - `Tip2`: Second tip.
  #'   - `RED`: RED value of their most recent common ancestor.
  #'
  #' @examples
  #' tree <- ape::read.tree(text = "(((A:1,B:1):1,(C:1,D:2):2):2,E:3);")
  #' reds <- castor::get_reds(tree)
  #' calculate_red_set(tree, c("A", "B", "C", "D"), reds)
  calculate_red_set <- function(tree, tips, reds=NULL) {
    if (!inherits(tree, "phylo")) 
      stop("object 'tree' is not of class 'phylo'")
    
    if (is.null(reds)) {
      reds <- castor::get_reds(tree)
    }
    
    # Initialize progress bar
    n <- length(tips)
    pb <- progress::progress_bar$new(
      format = "  Calculating REDs [:bar] :percent eta: :eta",
      total = n * (n - 1) / 2, # Total number of comparisons
      clear = FALSE, width = 60
    )
    
    # Get tip indices
    tip_indices <- match(tips, tree$tip.label)
    if (any(is.na(tip_indices))) stop("One or more tips not found in tree.")
    
    # Call C++ function, passing the progress bar environment
    pairwise_red <- calculate_red_set_cpp(
      edge = as.matrix(tree$edge), 
      reds = reds, 
      Ntip = length(tree$tip.label),
      tip_indices = tip_indices - 1,  # Adjust for 0-based indexing
      pb_env = pb
    )
    
    # Replace tip indices with labels
    pairwise_red$Tip1 <- pairwise_red$Tip1 + 1 # Adjust for 0-based indexing
    pairwise_red$Tip1 <- tree$tip.label[as.integer(pairwise_red$Tip1)]
    pairwise_red$Tip2 <- pairwise_red$Tip2 + 1 # Adjust for 0-based indexing
    pairwise_red$Tip2 <- tree$tip.label[as.integer(pairwise_red$Tip2)]
    
    return(pairwise_red)
  }
  
  #' Generate Chunks of Pairwise Indices for RED Calculation
  #'
  #' Splits all pairwise tip combinations into chunks for parallel processing.
  #'
  #' @param tip_indices An integer vector of tip indices.
  #' @param chunk_size An integer specifying the size of each chunk.
  #'
  #' @return A list of dataframes, where each dataframe contains two columns:
  #'   - `Var1`: Index of the first tip in a pair.
  #'   - `Var2`: Index of the second tip in a pair.
  #'
  #' @examples
  #' tips <- 1:5
  #' chunks <- generate_chunks(tips, 3)
  #' print(chunks)
  generate_chunks_RED <- function(tip_indices, chunk_size) {
    pairs <- expand.grid(Var1 = tip_indices, Var2 = tip_indices) %>%
      dplyr::filter(Var1 < Var2)  # Exclude diagonal and redundant pairs
    split(pairs, ceiling(seq_len(nrow(pairs)) / chunk_size))
  }
  
  #' Calculate Pairwise RED for a Chunk of Data
  #'
  #' Calculates pairwise Relative Evolutionary Divergence (RED) values for a set 
  #' of tips in a phylogenetic tree, dividing it into chunks for parallel processing, 
  #' and saving the results of each chunk to a specified output directory.
  #'
  #' @param tree A phylogenetic tree of class "phylo".
  #' @param tips A character vector of tip labels.
  #' @param reds A numeric vector of RED values for nodes in the tree. If NULL, they will be calculated.
  #' @param n_chunks Number of chunks for processing (default: 4).
  #' @param cpp_file A string specifying the path to the C++ file for calculating RED.
  #' @param chunk_dir A string specifying the directory to save chunk results as CSV files.
  #' @param workers Number of parallel workers to use (default: number of available cores).
  #' @param overwrite Logical indicating whether to overwrite existing chunk files (default: FALSE).
  #'
  #' @return This function does not return a value. It saves results as CSV files in the specified directory.
  #'
  #' @examples
  #' tree <- ape::read.tree(text = "(((A:1,B:1):1,(C:1,D:2):2):2,E:3);")
  #' reds <- castor::get_reds(tree)
  #' calculate_red_chunks(tree, tips = c("A", "B", "C", "D"), reds, n_chunks = 2, chunk_dir = "red_chunks")
  calculate_red_chunks <- function(tree, tips, reds = NULL, n_chunks = 4, cpp_file, chunk_dir, workers = parallelly::availableCores(), overwrite = FALSE) {
    if (!inherits(tree, "phylo")) stop("Tree must be a 'phylo' object.")
    if (is.null(chunk_dir)) stop("Output directory must be specified.")
    
    # Ensure RED values are available
    if (is.null(reds)) reds <- castor::get_reds(tree)
    
    # Create output directory if it does not exist
    if (!dir.exists(chunk_dir)) dir.create(chunk_dir, recursive = TRUE)
    
    # Get tip indices and generate chunks
    tip_indices <- match(tips, tree$tip.label)
    if (any(is.na(tip_indices))) stop("One or more tips not found in the tree.")
    chunk_size <- ceiling((length(tip_indices) * (length(tip_indices) - 1) / 2) / n_chunks)
    chunks <- generate_chunks_RED(tip_indices, chunk_size)
    
    # Identify already processed chunks
    if (!overwrite) {
      existing_files <- list.files(chunk_dir, pattern = "^chunk_\\d{3}\\.csv$", full.names = TRUE)
      processed_chunks <- as.numeric(gsub("chunk_(\\d{3})\\.csv", "\\1", basename(existing_files)))
      chunks_to_process <- setdiff(seq_along(chunks), processed_chunks)
      
      if (length(chunks_to_process) == 0) {
        message("All chunks have already been processed.")
        return(invisible(NULL))
      }
      
      message("Resuming calculation for the following chunks: ", paste(chunks_to_process, collapse = ", "))
    } else {
      # Process all chunks if overwrite is TRUE
      chunks_to_process <- seq_along(chunks)
      message("Overwrite mode enabled. Processing all chunks.")
    }
    
    # Set up parallel processing
    future::plan(future::multisession, workers = workers)
    
    # Calculate RED values for each chunk in parallel
    furrr::future_map(chunks_to_process, ~ {
      chunk_index <- .x
      chunk <- chunks[[.x]]
      
      # Source Rcpp file in each worker
      Rcpp::sourceCpp(cpp_file)
      
      # Adjust indices for 0-based indexing in C++
      chunk$Var1 <- chunk$Var1 - 1
      chunk$Var2 <- chunk$Var2 - 1
      
      # Process the chunk
      pairwise_df <- calculate_red_chunk_cpp(
        edge = as.matrix(tree$edge),
        reds = reds,
        Ntip = length(tree$tip.label),
        chunk = chunk
      )
      
      # Replace tip indices with labels
      pairwise_df$Tip1 <- pairwise_df$Tip1 + 1 # Adjust for 0-based indexing
      pairwise_df$Tip1 <- tree$tip.label[as.integer(pairwise_df$Tip1)]
      pairwise_df$Tip2 <- pairwise_df$Tip2 + 1 # Adjust for 0-based indexing
      pairwise_df$Tip2 <- tree$tip.label[as.integer(pairwise_df$Tip2)]
      
      # Save results to a CSV file
      chunk_file <- file.path(chunk_dir, sprintf("chunk_%03d.csv", chunk_index))
      write.csv(pairwise_df, chunk_file, row.names = FALSE)
      
      invisible(NULL)  # Suppress return value
    }, .options = furrr::furrr_options(seed = TRUE))
    
    # Stop parallel processing
    future::plan(future::sequential)
    
    # No return value
    invisible(NULL)
  }

  #' Combine Processed Chunks into a Single Object
  #'
  #' Combines the results of processed chunks from a specified output directory into a single dataframe or an Arrow dataset.
  #'
  #' @param chunk_dir A string specifying the directory containing processed chunk files as CSVs.
  #' @param use_arrow A logical indicating whether to create an Arrow dataset instead of a dataframe. Default is FALSE.
  #'
  #' @return Either a dataframe or an Arrow dataset object containing:
  #'   - `Seq1`: Name of the first sequence.
  #'   - `Seq2`: Name of the second sequence.
  #'   - `Identity`: Percent identity between the sequences.
  #'
  #' @examples
  #' # Assuming chunk files are saved in the "chunk_directory":
  #' df <- combine_processed_chunks("chunk_directory")
  #' print(df)
  #'
  #' # For large datasets using Arrow:
  #' dataset <- combine_processed_chunks("chunk_directory", use_arrow = TRUE)
  #' print(dataset)
  combine_processed_chunks <- function(chunk_dir, use_arrow = FALSE) {
    # Ensure the output directory exists
    if (!dir.exists(chunk_dir)) {
      stop("The specified output directory does not exist.")
    }
    
    # List all CSV files in the output directory
    chunk_files <- list.files(chunk_dir, pattern = "^chunk_\\d{3}\\.csv$", full.names = TRUE)
    
    # Check if there are any chunk files
    if (length(chunk_files) == 0) {
      stop("No chunk files found in the specified output directory.")
    }
    
    if (use_arrow) {
      # Ensure Arrow is installed
      if (!requireNamespace("arrow", quietly = TRUE)) {
        stop("The 'arrow' package is required for this functionality. Please install it using install.packages('arrow').")
      }
      
      # Create an Arrow dataset from the chunk files
      dataset <- arrow::open_dataset(chunk_files, format = "csv")
      return(dataset)
    } else {
      # Read and combine all chunk files into a single dataframe
      df <- chunk_files %>%
        lapply(read.csv) %>%
        dplyr::bind_rows()
      
      # Return the combined dataframe
      return(df)
    }
  }

  #' Test and Validate RED Calculations
  #'
  #' This function tests and validates calculations using a a phylogenetic tree with known
  #' values of RED. It evaluates calculate_RED_pair(), 
  #' calculate_RED_set(), calculate_RED_chunks(), and dependent functions. 
  #' The test tree is from https://www.nature.com/articles/nbt.4229/figures/1.
  #' A temporary test directory is created for intermediate files and removed 
  #' automatically after the tests.
  #'
  #' @param cpp_file A string specifying the path to the compiled C++ file containing the implementation for RED calculation.
  #' @param test_dir Optional. A string specifying the path to a test directory for intermediate files. If not provided, 
  #'   a directory named "test_[random_string]" is created in the current working directory.
  #'
  #' @return A list containing the following elements:
  #'   \describe{
  #'     \item{\code{pair_result}}{The RED value for a pair of tips.}
  #'     \item{\code{set_result}}{A dataframe of pairwise RED values for all tips in the tree.}
  #'     \item{\code{chunk_result}}{The combined results of chunked RED calculations.}
  #'   }
  #'
  #' @details
  #' The function creates a temporary test directory to store intermediate results for chunked RED calculations.
  #' This directory is removed at the end of the function execution, even if an error occurs.
  #'
  #' The following tests are performed:
  #' \itemize{
  #'   \item Pairwise RED calculation for two tips in the tree.
  #'   \item Pairwise RED calculation for all tips in the tree.
  #'   \item Chunked RED calculation, with results combined for validation.
  #' }
  #'
  #' @examples
  #' cpp_file <- "path/to/your/cpp_file.cpp"
  #' test_results <- test_red_calculation(cpp_file)
  #' print(test_results$pair_result)  # Pairwise RED result
  #' print(test_results$set_result)  # All pairwise RED results
  #' print(test_results$chunk_result)  # Combined chunked RED results
  #'
  #' @export
  #' Test and Validate RED Calculations
  #'
  #' This function tests and validates calculations using a a phylogenetic tree with known
  #' values of RED. It evaluates calculate_RED_pair(), 
  #' calculate_RED_set(), calculate_RED_chunks(), and dependent functions. 
  #' The test tree is from https://www.nature.com/articles/nbt.4229/figures/1.
  #' A temporary test directory is created for intermediate files and removed 
  #' automatically after the tests.
  #'
  #' @param cpp_file A string specifying the path to the compiled C++ file containing the implementation for RED calculation.
  #' @param test_dir Optional. A string specifying the path to a test directory for intermediate files. If not provided, 
  #'   a directory named "test_[random_string]" is created in the current working directory.
  #'
  #' @return A list containing the following elements:
  #'   \describe{
  #'     \item{\code{pair_result}}{The RED value for a pair of tips.}
  #'     \item{\code{set_result}}{A dataframe of pairwise RED values for all tips in the tree.}
  #'     \item{\code{chunk_result}}{The combined results of chunked RED calculations.}
  #'   }
  #'
  #' @details
  #' The function creates a temporary test directory to store intermediate results for chunked RED calculations.
  #' This directory is removed at the end of the function execution, even if an error occurs.
  #'
  #' The following tests are performed:
  #' \itemize{
  #'   \item Pairwise RED calculation for two tips in the tree.
  #'   \item Pairwise RED calculation for all tips in the tree.
  #'   \item Chunked RED calculation, with results combined for validation.
  #' }
  #'
  #' @examples
  #' cpp_file <- "path/to/your/cpp_file.cpp"
  #' test_results <- test_red_calculation(cpp_file)
  #' print(test_results$pair_result)  # Pairwise RED result
  #' print(test_results$set_result)  # All pairwise RED results
  #' print(test_results$chunk_result)  # Combined chunked RED results
  #'
  #' @export
  test_red_calculation <- function(cpp_file, test_dir = NULL) {
    test_results <- list()  # Track test results
    
    # Generate a test directory if not provided
    if (is.null(test_dir)) {
      random_string <- stringr::str_sub(stringr::str_c(sample(c(letters, LETTERS, 0:9), 8, replace = TRUE), collapse = ""), 1, 8)
      test_dir <- file.path(getwd(), paste0("test_", random_string))
    }
    
    # Create the test directory
    if (!dir.exists(test_dir)) dir.create(test_dir, recursive = TRUE)
    
    # Ensure the directory is removed at the end
    on.exit({
      if (dir.exists(test_dir)) {
        unlink(test_dir, recursive = TRUE)
        cat("Test directory removed:", test_dir, "\n")
      }
    }, add = TRUE)
    
    cat("Test directory created:", test_dir, "\n")
    
    # Get test data
    newick_string <- "(((A:1,B:1):1,(C:1,D:2):2):2,E:3);"
    tree <- ape::read.tree(text = newick_string)
    
    # Get cpp file
    Rcpp::sourceCpp(cpp_file)
    
    # Test RED pair calculation
    cat("\n=== Testing calculate_red_pair ===\n")
    result1 <- calculate_red_pair(tree, tip1 = "A", tip2 = "B")
    expected1 <- 0.7105263
    print(result1)
    if (all.equal(result1, expected1, tolerance = 1e-6)) {
      test_results$calculate_RED_pair <- TRUE
      cat("Output matches expected value for calculate_red_pair.\n\n")
    } else {
      test_results$calculate_RED_pair <- FALSE
      cat("Output does NOT match expected value for calculate_red_pair!\n")
      cat("Expected:\n")
      print(expected1)
      cat("Actual:\n")
      print(result1)
    }
    
    # Test RED set calculation
    cat("\n=== Testing calculate_red_set ===\n")
    result2 <- calculate_red_set(tree = tree, tips = tree$tip.label)
    expected_result2 <- data.frame(
      Tip1 = c("A", "A", "A", "A", "B", "B", "B", "C", "C", "D"),
      Tip2 = c("B", "C", "D", "E", "C", "D", "E", "D", "E", "E"),
      RED = c(0.7105263, 0.4210526, 0.4210526, 0.0000000, 0.4210526,
              0.4210526, 0.0000000, 0.7518797, 0.0000000, 0.0000000)
    )
    print(result2)
    if (all.equal(result2, expected_result2, tolerance = 1e-6)) {
      test_results$calculate_RED_set <- TRUE
      cat("Output matches expected values for calculate_red_set.\n\n")
    } else {
      test_results$calculate_RED_set <- FALSE
      cat("Output does NOT match expected values for calculate_red_set!\n")
      cat("Expected:\n")
      print(expected_result2)
      cat("Actual:\n")
      print(result2)
    }
    
    # Test chunked RED calculation
    cat("\n=== Testing calculate_red_chunks ===\n")
    calculate_red_chunks(tree, tips = tree$tip.label, n_chunks = 2, cpp_file = cpp_file, chunk_dir = test_dir)
    
    # Combine chunks
    cat("\n=== Combining processed chunks ===\n")
    result3 <- combine_processed_chunks(chunk_dir = test_dir)
    print(result3 %>% dplyr::arrange(Tip1))
    if (all.equal(result3 %>% dplyr::arrange(Tip1), expected_result2, tolerance = 1e-6)) {
      test_results$calculate_RED_chunks <- TRUE
      cat("Output matches expected values for combined chunks.\n\n")
    } else {
      test_results$calculate_RED_chunks <- FALSE
      cat("Output does NOT match expected values for combined chunks!\n")
      cat("Expected:\n")
      print(expected_result2)
      cat("Actual:\n")
      print(result3 %>% dplyr::arrange(Tip1))
    }
    
    # Final summary
    if (all(unlist(test_results))) {
      cat("\nAll function outputs match expected values.\n")
    } else {
      failed_tests <- names(test_results)[!unlist(test_results)]
      cat("\nOutput does not match expected values for the following functions: ", paste(failed_tests, collapse = ", "), "\n")
    }
    
    invisible(list(pair_result = result1, set_result = result2, chunk_result = result3))
  }
  
#=== Install packages ===
  # Define required packages
  cran_packages <- c("arrow", "ape", "castor", "import", "readr", "Rcpp", "dplyr", "import", "future", "furrr", "stringr", "magrittr")
  
  # Install missing packages
  install_missing_cran_packages(cran_packages)

  # Set package options
  options(future.globals.maxSize = 8 * 1024^3)  # Set limit to 8 GB  
  
# === Import functions ===
  import::from(magrittr, "%>%")
 
# === Define file paths ===
  base_path <- "C:\\My Directory" # Set to actual directory
  RED_cpp_file <- file.path(base_path, "calculateRelativeEvolutionaryDivergence.cpp") 
  
#=== Test functions === 
  test_red_calculation(cpp_file = RED_cpp_file)
  