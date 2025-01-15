# Get 16S Ribosomal Sequences from LPSN
#
# This script downloads the 16S rDNA sequences for type strains of bacterial species from the LPSN database.
#
# Author: Timothy Hackmann
# Date: 15 January 2025

# === Define functions ===
  #' Retrieve the Body of a Web Page
  #'
  #' This function retrieves the body content of a web page using the `polite` package.
  #'
  #' @param url A character string representing the URL of the web page to scrape.
  #' @param user_agent A character string specifying the user agent to use for the request. 
  #'        Default is `"bingbot"`.
  #'
  #' @return An object containing the HTML body of the web page.
  #' @export
  #'
  #' @examples
  #' # Example usage:
  #' body <- get_web_page_body("https://example.com", user_agent = "bingbot")
  #'
  get_web_page_body <- function(url, user_agent = "bingbot") {
    # Initiate a session with the specified URL
    session <- polite::bow(url, user_agent = user_agent)
    
    # Scrape the body of the page
    body <- polite::scrape(bow = session)
    
    return(body)
  }
  
  #' Extract Links from Web Page Body
  #'
  #' This function extracts all hyperlinks (URLs) from the body of a web page.
  #'
  #' @param body An HTML document object, typically returned by functions like `polite::scrape()` or `rvest::read_html()`.
  #' @param tag A character string specifying the HTML tag to search for. Default is `"a"`.
  #' @param attribute A character string specifying the attribute to extract. Default is `"href"`.
  #'
  #' @return A character vector of extracted links.
  #' @export
  #'
  #' @examples
  #' # Example usage:
  #' body <- polite::scrape(polite::bow("https://example.com"))
  #' links <- extract_links(body)
  #'
  extract_links <- function(body, tag = "a", attribute = "href") {
    # Extract specified attributes from the specified HTML tags
    links <- body %>%
      rvest::html_nodes(tag) %>%
      rvest::html_attr(attribute)
    
    return(links)
  }
  
  #' Download a FASTA file from a given link
  #'
  #' This function downloads a FASTA file from the provided URL and saves it to the specified file path.
  #' If the download fails, an error message is returned.
  #'
  #' @param fasta_link A character string containing the URL of the FASTA file to download.
  #' @param fp A character string specifying the file path where the downloaded FASTA file will be saved.
  #'
  #' @return A list with the following elements:
  #' \describe{
  #'   \item{success}{A logical value indicating whether the download was successful.}
  #'   \item{filepath}{The file path where the FASTA file was saved (NULL if unsuccessful).}
  #'   \item{message}{A message indicating the status of the operation (e.g., success or error).}
  #' }
  #'
  #' @examples
  #' # Example usage
  #' fasta_link <- "https://example.com/sample.fasta"
  #' fp <- "path/to/save/sample.fasta"
  #' result <- download_fasta_file(fasta_link, fp)
  #' print(result$message)
  #'
  #' @export
  download_fasta_file <- function(fasta_link, fp) {
    tryCatch({
      # Download and save the file
      fasta_file <- httr::GET(fasta_link)
      writeBin(httr::content(fasta_file, "raw"), fp)
      
      list(success = TRUE, filepath = fp, message = paste("Downloaded:", fp))
    }, error = function(e) {
      list(success = FALSE, filepath = NULL, message = paste("Error:", e$message))
    })
  }
  
  #' Combine DNAStringSet Objects
  #'
  #' This function combines a list of `DNAStringSet` objects into a single `DNAStringSet`.
  #'
  #' @param fasta_list A list of `DNAStringSet` objects to be combined.
  #'
  #' @return A combined `DNAStringSet` object containing all sequences and their corresponding IDs.
  #' @import Biostrings ShortRead
  #' @export
  #'
  #' @examples
  #' \dontrun{
  #' combined_fasta <- combine_fasta(fasta_list)
  #' }
  combine_fasta <- function(fasta_list) {
    # Combine sequences
    combined_reads <- do.call(c, lapply(fasta_list, ShortRead::sread))
    
    # Combine IDs
    combined_ids <- do.call(c, lapply(fasta_list, ShortRead::id))
    
    # Create a new ShortRead object
    ShortRead::ShortRead(sread = combined_reads, id = combined_ids)
  }

  
# === Import functions ===
  import::from(magrittr, "%>%")
  
# === Set file paths ===
  base_path <- "C:\\My Directory" # Set to actual directory
  lpsn_file <- file.path(base_path, "data\\lpsn_gss_2024-10-17.csv") # From https://lpsn.dsmz.de/downloads
  
# === Read in data ===
  lpsn_data <- utils::read.csv(lpsn_file)

# === Get metadata for species ===
  species_data = lpsn_data
  species_data = species_data %>% dplyr::filter(grepl(pattern = "correct name", x = status, ignore.case = TRUE))
  species_data = species_data %>% dplyr::filter(sp_epithet!="")
  
  # For subspecies, keep only entries that have subspecies specified (e.g., keep Selenomonas ruminantium lactilytica but not Selenomonas ruminantium)
  species_data <- species_data %>% dplyr::group_by(genus_name, sp_epithet) %>% dplyr::filter(!(dplyr::n_distinct(subsp_epithet) > 1 & subsp_epithet == "")) %>% dplyr::ungroup()
  
# === Download ribosomal sequences ===
  # Handle most strains (those that are not type subspecies)
    most_strains <- species_data %>% dplyr::filter(sp_epithet!=subsp_epithet)
    addresses <- most_strains$address
    LPSN_ID <- most_strains$record_no
    base_url <- "https://lpsn.dsmz.de"
    download_status <- vector("character", length(addresses))
    
    for (i in seq_along(addresses)) {
      # Get body of the page
      body <- get_web_page_body(url = addresses[i], user_agent = "me")
      
      # Extract link to FASTA file
      links <- extract_links(body)
      fasta_link <- links[grepl("\\.fasta$", links)]
      fasta_link <- ifelse(grepl("^http", fasta_link), fasta_link, paste0(base_url, fasta_link))
      
      if (length(fasta_link) > 0) {
        # Construct the file path
        fp <- paste0(base_path, "data\\", LPSN_ID[i], ".fasta")
        
        # Download the FASTA file
        download_result <- download_fasta_file(fasta_link[1], fp)
        download_status[i] <- download_result$message} else {
        download_status[i] <- "No FASTA file found"
      }
      
      # Show progress
      svMisc::progress(value = i, max.value = length(addresses))
    }
  
  # Handle strains of type subspecies
    ## These strains do not have sequences on their own page--instead they are on the page of the parent taxon
    type_subspecies <- species_data %>% dplyr::filter(sp_epithet==subsp_epithet)
    addresses <- type_subspecies$address
    LPSN_ID <- type_subspecies$record_no
    base_url <- "https://lpsn.dsmz.de"
    download_status <- vector("character", length(addresses))
    
    for (i in seq_along(addresses)) {
      # Get body of the page
      body <- get_web_page_body(url = addresses[i], user_agent = "me")
      
      # Extract parent taxon link
      parent_taxon_link <- body %>%
        rvest::html_nodes(xpath = "//b[contains(text(), 'Parent taxon:')]/following-sibling::a[1]") %>%
        rvest::html_attr("href")
      
      if (!is.na(parent_taxon_link) && length(parent_taxon_link) > 0) {
        # Get URL of parent taxon
        parent_url <- paste0(base_url, parent_taxon_link)
        body <- get_web_page_body(url = parent_url, user_agent = "me")
        
        # Extract link to FASTA file
        links <- extract_links(body)
        fasta_link <- links[grepl("\\.fasta$", links)]
        fasta_link <- ifelse(grepl("^http", fasta_link), fasta_link, paste0(base_url, fasta_link))
        
        if (length(fasta_link) > 0) {
          # Construct the file path
          fp <- paste0(base_path, "data\\", LPSN_ID[i], ".fasta")
          
          # Download the FASTA file
          download_result <- download_fasta_file(fasta_link[1], fp)
          download_status[i] <- download_result$message
        } else {
          download_status[i] <- "No FASTA file found"
        }
      } else {
        download_status[i] <- "No parent taxon link found"
      }
      
      # Show progress
      svMisc::progress(value = i, max.value = length(addresses))
    }

# === Remove extra sequences ===    
    # Type subspecies
    # Extra sequences are those that contain "subsp."
    LPSN_ID <- type_subspecies$record_no
    
    for (i in 1:length(LPSN_ID)) {
      fp <- paste0(base_path, "data\\", LPSN_ID[i], ".fasta")
      
      # Check if the file exists
      if (!file.exists(fp)) {
        message(paste("File not found:", fp, "- Exiting loop."))
        next
      }
      
      tryCatch({
        fasta_content <- Biostrings::readDNAStringSet(fp, format = "fasta")
        
        if (length(fasta_content) > 0) {
          # Filter out sequences with "subsp." in their IDs
          filtered_sequences <- fasta_content[!grepl("subsp\\.", names(fasta_content))]
          
          if (length(filtered_sequences) > 0) {
            # Save the filtered sequences back to the file
            Biostrings::writeXStringSet(filtered_sequences, filepath = fp, format = "fasta")
            message(paste("Processed and saved:", fp))
          } else {
            # Remove the file if no valid sequences remain
            file.remove(fp)
            message(paste("No valid sequences found in:", fp, "- File removed."))
          }
        }
      }, error = function(e) {
        message(paste("Error processing file:", fp, "-", e$message))
      })
    }
    
    # Most strains
    # Extra sequences are those that contain "subsp.", if there are two or more sequences
    LPSN_ID <- most_strains$LPSN_ID
    
    for (i in 1:length(LPSN_ID)) {
      fp <- paste0(base_path, "data\\", LPSN_ID[i], ".fasta")
      
      # Check if the file exists
      if (!file.exists(fp)) {
        message(paste("File not found:", fp, "- Exiting loop."))
        next
      }
      
      tryCatch({
        fasta_content <- Biostrings::readDNAStringSet(fp, format = "fasta")
        
        if (length(fasta_content) > 1) {
          # Filter out sequences with "subsp." in their IDs
          filtered_sequences <- fasta_content[!grepl("subsp\\.", names(fasta_content))]
          
          if (length(filtered_sequences) > 0) {
            # Save the filtered sequences back to the file
            Biostrings::writeXStringSet(filtered_sequences, filepath = fp, format = "fasta")
            message(paste("Processed and saved:", fp))
          } else {
            # Remove the file if no valid sequences remain
            file.remove(fp)
            message(paste("No valid sequences found in:", fp, "- File removed."))
          }
        }
      }, error = function(e) {
        message(paste("Error processing file:", fp, "-", e$message))
      })
    }
    
# === Put ribosomal sequences in a single object ===
  # Get names of FASTA files downloaded above
  fasta_files <- list.files(paste0(base_path, "data\\"), 
                            pattern = "\\.fasta$", full.names = TRUE)
  
  # Read each FASTA file and return as a list
  fasta_list <- lapply(seq_along(fasta_files), function(i) {
    svMisc::progress(i, max.value = length(fasta_files))
    return(ShortRead::readFasta(fasta_files[i]))
  })
  
  # Rename FASTA files according to file name (record)
  fasta_names <- list.files(paste0(base_path, "data\\"), 
                            pattern = "\\.fasta$", full.names = FALSE)
  fasta_names <- gsub(pattern="\\.fasta$", replacement="", x = fasta_names)
  
  fasta_list_renamed <-  lapply(seq_along(fasta_list), function(i) {
    fasta <- fasta_list[[i]]
    new_ids  <- Biostrings::BStringSet(rep(fasta_names[i], length(ShortRead::sread(fasta))))
    ShortRead::ShortRead(sread = ShortRead::sread(fasta), id = new_ids)
  })
  
  # Combine files
  combined_fasta <- combine_fasta(fasta_list_renamed)
  
  # Put in DNAStringSet
  seq = ShortRead::sread(combined_fasta)
  names(seq) = ShortRead::id(combined_fasta)

  
# === Add ribosomal sequences to organism data ===  
  idx <- match(species_data$record_no, names(seq))
  species_data$`16S_ribosomal_sequence` <- NA
  species_data$`16S_ribosomal_sequence`[!is.na(idx)] <- as.character(seq[idx[!is.na(idx)]])
  
# === Export ===
  Biostrings::writeXStringSet(seq, filepath = paste0(base_path, "data\\16S_ribosomal_sequences.fasta"), format = "fasta")
  write.csv(species_data, paste0(base_path, "data\\lpsn_ribosomal_sequences.csv"))
  