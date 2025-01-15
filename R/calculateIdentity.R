# Calculate Identity Over A Set Of Sequences
#
# This script calculates pairwise identity among a set of sequences
# It can calculate this in parallel by breaking the set into chunks
# Sequences are aligned using the Needleman-Wunsch algorithm with scores of 
# 2 for match, -1 for mismatch, -10 for gap opening, and -0.5 for gap extension
# Identity is calculated as the number of letters that match after alignment, 
# divided by the number of letters in the alignment.  By default, gaps are 
# ignored and do not contribute to the number of matching letters or length of alignment.
# The script calls a cpp file that in turn uses SeqAn (v 2.4) to do calculations.
#
# Author:  Timothy Hackmann
# Date: 15 January 2025

# === Define functions ===
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
  
  #' Calculate Pairwise Identity of Two Sequences
  #'
  #' Aligns two sequences and calculates their percent identity, alignment score, and other alignment metrics.
  #'
  #' @param seq1 A string representing the first sequence.
  #' @param seq2 A string representing the second sequence.
  #' @param match An integer score for matching characters (default: 2).
  #' @param mismatch An integer penalty for mismatching characters (default: -1).
  #' @param gap_open An integer penalty for opening a gap (default: -10).
  #' @param gap_extend An integer penalty for extending a gap (default: -0.5).
  #' @param ignore_gap A logical indicating whether to ignore gaps when calculating percent identity (default: TRUE).
  #'
  #' @return A list containing:
  #'   - `alignment_score`: Alignment score.
  #'   - `alignment_string`: String representation of the alignment.
  #'   - `aligned_sequence1`: Aligned version of `seq1`.
  #'   - `aligned_sequence2`: Aligned version of `seq2`.
  #'   - `percent_identity`: Percent identity between the sequences.
  #'   - `match_count`: Number of matches.
  #'   - `mismatch_count`: Number of mismatches.
  #'   - `total_aligned`: Total number of aligned positions.
  #'
  #' @examples
  #' seq1 <- "ATGCGT"
  #' seq2 <- "ATGCGC"
  #' result <- calculate_identity_pair(seq1, seq2)
  #' print(result$percent_identity)
  calculate_identity_pair <- function(seq1, seq2, match = 2, mismatch = -1, gap_open = -10, gap_extend = -0.5, ignore_gap = TRUE) {
    # Perform alignment using the C++ function
    alignment_result <- perform_alignment(seq1, seq2, match, mismatch, gap_open, gap_extend)
    
    # Extract aligned sequences
    aligned_row1 <- alignment_result$aligned_row1
    aligned_row2 <- alignment_result$aligned_row2
    
    # Calculate percent identity
    identity_result <- calculate_identity_from_alignment(aligned_row1, aligned_row2, ignore_gap)
    
    # Return results
    list(
      alignment_score = alignment_result$alignment_score,
      alignment_string = alignment_result$alignment,
      aligned_sequence1 = aligned_row1,
      aligned_sequence2 = aligned_row2,
      percent_identity = identity_result$percent_identity,
      match_count = identity_result$match_count,
      mismatch_count = identity_result$mismatch_count,
      total_aligned = identity_result$totalAligned
    )
  }
  
  #' Calculate Pairwise Identities for Set of Sequences
  #'
  #' Calculates pairwise identities for all sequences in a dataset
  #'
  #' @param dna_set A named vector or DNAStringSet of sequences.
  #'
  #' @return A dataframe containing:
  #'   - `Seq1`: Name of the first sequence.
  #'   - `Seq2`: Name of the second sequence.
  #'   - `Identity`: Percent identity between the sequences.
  #'
  #' @examples
  #' dna_set <- c("ATGCGT" = "ATGCGT", "ATGCGC" = "ATGCGC", "TTGCAA" = "TTGCAA")
  #' result <- calculate_identity_set(dna_set)
  #' print(result)
  calculate_identity_set <- function(dna_set) {
    # Initialize progress bar in R
    n <- length(dna_set)
    pb <- progress::progress_bar$new(
      format = "  Aligning [:bar] :percent eta: :eta",
      total = n * (n - 1) / 2, # Total number of comparisons
      clear = FALSE, width = 60
    )
    
    # Call C++ function, passing the progress bar environment
    pairwise_identities <- calculate_identity_set_cpp(as.character(dna_set), pb)
    
    # Replace Var1 and Var2 with sequence names and rename columns
    pairwise_identities <- pairwise_identities %>%
      dplyr::mutate(
        Seq1 = names(dna_set)[Var1],
        Seq2 = names(dna_set)[Var2]
      ) %>%
      dplyr::select(Seq1, Seq2, Identity = Freq)
    
    # Return results as a dataframe
    return(pairwise_identities)
  }
  
  #' Generate Chunks of Pairwise Indices for Identity Calculation
  #'
  #' Splits all pairwise sequence combinations into chunks for parallel processing.
  #'
  #' @param n An integer specifying the number of sequences.
  #' @param chunk_size An integer specifying the size of each chunk.
  #'
  #' @return A list of dataframes, where each dataframe contains two columns:
  #'   - `Var1`: Index of the first sequence in a pair.
  #'   - `Var2`: Index of the second sequence in a pair.
  #'
  #' @examples
  #' chunks <- generate_chunks(4, 2)
  #' print(chunks)
  generate_chunks_identity <- function(n, chunk_size) {
    pairs <- expand.grid(Var1 = seq_len(n), Var2 = seq_len(n)) %>%
      dplyr::filter(Var1 < Var2)  # Exclude diagonal and redundant pairs
    split(pairs, ceiling(seq_len(nrow(pairs)) / chunk_size))
  }
  
  #' Calculate Pairwise Identities for a Chunk of Data
  #'
  #' Processes pairwise identity calculations for a given dataset, dividing it into chunks
  #' for parallel processing, and saving the results of each chunk to a specified output directory.
  #'
  #' @param dna_set A named vector or DNAStringSet of sequences.
  #' @param n_chunks An integer specifying the number of chunks to divide the dataset into (default: 4).
  #' @param cpp_file A string specifying the path to the C++ file for alignment calculations.
  #' @param chunk_dir A string specifying the directory to save chunk results as CSV files.
  #' @param workers An integer specifying the number of parallel workers to use (default: the number of available cores).
  #' @param overwrite A logical value indicating whether to overwrite previously processed chunks (default: FALSE).
  #'
  #' @return This function does not return a value. It saves results as CSV files in the specified directory.
  #'
  #' @examples
  #' dna_set <- c("ATGCGT" = "ATGCGT", "ATGCGC" = "ATGCGC", "TTGCAA" = "TTGCAA")
  #' calculate_identity_chunks(
  #'   dna_set = dna_set,
  #'   n_chunks = 4,
  #'   cpp_file = "path/to/cpp/file.cpp",
  #'   chunk_dir = "chunk_directory",
  #'   overwrite = TRUE
  #' )
  calculate_identity_chunks <- function(dna_set, n_chunks = 4, cpp_file, chunk_dir, workers = parallelly::availableCores(), overwrite = FALSE) {
    # Ensure output directory is specified
    if (is.null(chunk_dir)) {
      stop("An output directory must be specified to save the results.")
    }
    
    # Ensure output directory exists
    if (!dir.exists(chunk_dir)) dir.create(chunk_dir, recursive = TRUE)
    
    # Split sequences into chunks for processing
    n <- length(dna_set)
    chunk_size <- ceiling((n * (n - 1) / 2) / n_chunks)
    chunks <- generate_chunks_identity(n, chunk_size)
    
    # Determine already processed chunks
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
    
    # Calculate pairwise identities for sequences in each chunk
    furrr::future_map(chunks_to_process, ~ {
      chunk_index <- .x  
      
      # Source Rcpp file in each worker
      Rcpp::sourceCpp(cpp_file)
      
      # Process the chunk
      df <- calculate_identity_chunk_cpp(as.character(dna_set), chunks[[chunk_index]])
      
      # Replace Var1 and Var2 with sequence names and rename columns
      df <- df %>%
        dplyr::mutate(
          Seq1 = names(dna_set)[Var1],
          Seq2 = names(dna_set)[Var2]
        ) %>%
        dplyr::select(Seq1, Seq2, Identity = Freq)
      
      # Save the result to a file using the original chunk index
      chunk_file <- file.path(chunk_dir, sprintf("chunk_%03d.csv", chunk_index))
      write.csv(df, chunk_file, row.names = FALSE)
      
      # Suppress return value
      invisible(NULL)
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

  #' Test and Validate Identity Calculations
  #'
  #' This function tests and validates calculations using a set of DNA sequences with
  #' known values of identity. It evaluates calculate_identity_pair(), 
  #' calculate_identity_set(), calculate_identity_chunks(), and dependent functions. 
  #' The test sequences are short strings of nucleotides, and the 
  #' expected identity values have been verified by hand.  
  #' A temporary test directory is created for intermediate files and removed 
  #' automatically after the tests.
  #'
  #' @param cpp_file A string specifying the path to the compiled C++ file containing the implementation of identity calculation.
  #' @param test_dir Optional. A string specifying the path to a test directory for intermediate files. If not provided, 
  #'   a directory named "test_[random_string]" is created in the current working directory.
  #'
  #' @return A list containing the following elements:
  #'   \describe{
  #'     \item{\code{pair_result}}{The result of pairwise identity calculation for two sequences.}
  #'     \item{\code{set_result}}{The result of pairwise identity calculation for all sequences.}
  #'     \item{\code{chunk_result}}{The combined results of chunked identity calculations.}
  #'   }
  #'
  #' @details
  #' The function creates a temporary test directory to store intermediate results for chunked identity calculations.
  #' This directory is removed at the end of the function execution, even if an error occurs. 
  #'
  #' The following tests are performed:
  #' \itemize{
  #'   \item Pairwise identity calculation of two sequences.
  #'   \item Pairwise identity calculation of all sequences.
  #'   \item Chunked pairwise identity calculation, with results combined for validation.
  #' }
  #'
  #' @examples
  #' cpp_file <- "path/to/your/cpp_file.cpp"
  #' test_results <- test_identity_calculation(cpp_file)
  #' print(test_results$pair_result)  # Pairwise identity result
  #' print(test_results$set_result)  # All sequence pairwise identity result
  #' print(test_results$chunk_result)  # Chunked identity calculation result
  #'
  #' @export
  test_identity_calculation <- function(cpp_file, test_dir = NULL) {
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
    dna_set <- Biostrings::DNAStringSet(
      c(
        "1" = "ATGCTGCGC",
        "2" = "ATGCGCGC",
        "3" = "TTGCAA",
        "4" = "TTGCCG"
      )
    )
    
    # Get cpp file
    Rcpp::sourceCpp(cpp_file)
    
    # Test pairwise identity of two sequences
    cat("\n=== Testing calculate_identity_pair ===\n")
    seq1 <- as.character(dna_set[1])
    seq2 <- as.character(dna_set[2])
    result1 <- calculate_identity_pair(seq1 = seq1, seq2 = seq2, ignore_gap = TRUE)
    
    # Expected output for pairwise identity
    expected_result1 <- list(
      alignment_score = 6,
      alignment_string = "      0     .     \n        ATGCTGCGC\n        |||| ||||\n        ATGC-GCGC\n\n\n",
      aligned_sequence1 = "ATGCTGCGC",
      aligned_sequence2 = "ATGC-GCGC",
      percent_identity = 100,
      match_count = 8,
      mismatch_count = 0,
      total_aligned = 8
    )
    
    print(result1)
    if (all.equal(result1, expected_result1, tolerance = 1e-6)) {
      test_results$calculate_identity_pair <- TRUE
      cat("Output matches expected value for calculate_identity_pair.\n\n")
    } else {
      test_results$calculate_identity_pair <- FALSE
      cat("Output does NOT match expected value for calculate_identity_pair!\n")
      cat("Expected:\n")
      print(expected_result1)
      cat("Actual:\n")
      print(result1)
    }
    
    # Test pairwise identity of all sequences
    cat("\n=== Testing calculate_identity_set ===\n")
    result2 <- calculate_identity_set(dna_set)
        result2 <- result2 %>%
      dplyr::mutate(Identity = as.character(round(as.numeric(Identity), 5)))
    result2 <- dplyr::mutate(result2, dplyr::across(everything(), as.character))

    # Expected output for identity set
    expected_result2 <- data.frame(
      Seq1 = c(1, 1, 1, 2, 2, 3),
      Seq2 = c(2, 3, 4, 3, 4, 4),
      Identity = c(100.00000, 50.00000, 66.66667, 50.00000, 50.00000, 66.66667)
    )
    expected_result2 <- expected_result2 %>%
      dplyr::mutate(Identity = as.character(round(as.numeric(Identity), 5)))
    expected_result2 <- dplyr::mutate(expected_result2, dplyr::across(everything(), as.character))
    
    print(result2 %>% dplyr::arrange(Seq1, Seq2))
    if (isTRUE(all.equal(result2 %>% dplyr::arrange(Seq1, Seq2), expected_result2, tolerance = 1e-6))) {
      test_results$calculate_identity_set <- TRUE
      cat("Output matches expected values for calculate_identity_set.\n\n")
    } else {
      test_results$calculate_identity_set <- FALSE
      cat("Output does NOT match expected values for calculate_identity_set!\n")
      cat("Expected:\n")
      print(expected_result2)
      cat("Actual:\n")
      print(result2 %>% dplyr::arrange(Seq1, Seq2))
    }
    
    # Test chunked identity calculations
    cat("\n=== Testing calculate_identity_chunks ===\n")
    calculate_identity_chunks(dna_set, n_chunks = 2, cpp_file = cpp_file, chunk_dir = test_dir)
    
    # Combine chunks
    cat("\n=== Combining processed chunks ===\n")
    result3 <- combine_processed_chunks(chunk_dir = test_dir)
    
    # Ensure consistent column types
    result3 <- result3 %>%
      dplyr::mutate(
        Seq1 = as.numeric(Seq1),
        Seq2 = as.numeric(Seq2),
        Identity = as.numeric(Identity)
      )
    
    expected_result2 <- expected_result2 %>%
      dplyr::mutate(
        Seq1 = as.numeric(Seq1),
        Seq2 = as.numeric(Seq2),
        Identity = as.numeric(Identity)
      )
    
    print(result3 %>% dplyr::arrange(Seq1, Seq2))
    if (isTRUE(all.equal(result3 %>% dplyr::arrange(Seq1, Seq2), expected_result2, tolerance = 1e-6))) {
      test_results$calculate_identity_chunks <- TRUE
      cat("Output matches expected values for combined chunks.\n\n")
    } else {
      test_results$calculate_identity_chunks <- TRUE
      cat("Output does NOT match expected values for combined chunks!\n")
      cat("Expected:\n")
      print(expected_result2)
      cat("Actual:\n")
      print(result3 %>% dplyr::arrange(Seq1, Seq2))
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
  cran_packages <- c("arrow", "import", "Rcpp", "dplyr", "future", "furrr")
  bioc_packages <- c("RSeqAn", "Biostrings")
  
  # Install missing packages
  install_missing_cran_packages(cran_packages)
  install_missing_bioc_packages(bioc_packages)
  
  # Set package options
  options(future.globals.maxSize = 8 * 1024^3)  # Set limit to 8 GB
  
# === Import functions ===
  import::from(magrittr, "%>%")
  
# === Set file paths ===
  base_path <- "C:\\My Directory" # Set to actual directory
  identity_cpp_file <- file.path(base_path, "calculateIdentity.cpp")
  identity_chunk_dir <- file.path(base_path, "chunks\\")
  fasta_file <- file.path(base_path, "data\\16S_ribosomal_sequences.fasta")

#=== Test functions === 
  test_identity_calculation(cpp_file = identity_cpp_file)
  
#===  Calculate identities ===
  # Load DNA sequences from a FASTA file
  dna_set <- Biostrings::readDNAStringSet(fasta_file)
  
  # Perform calculation
  calculate_identity_chunks(dna_set, n_chunks = 1000, cpp_file = identity_cpp_file, chunk_dir = identity_chunk_dir)
  
  # Combine chunks
  result <- combine_processed_chunks(chunk_dir = identity_chunk_dir, use_arrow = TRUE)