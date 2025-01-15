# Get Phylogeny from LPSN
#
# This script retrieves phylogeny of type strains from the LPSN database.
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
  
  #' Extract Taxonomic Ranks from Links
  #'
  #' This function extracts Genus, Family, and other taxonomic ranks from a vector of links.
  #' The links are from the LPSN and follow the format `"/rank/name"`
  #'
  #' @param phylogeny A character vector of links
  #' @param ranks A character vector of  ranks to extract
  #' @return A named character vector with the extracted taxonomic names for each rank. If a rank is not 
  #'         present in the phylogeny list, its value will be `NA`.
  #'
  #' @examples
  #' # Example phylogeny input
  #' phylogeny_example <- c(
  #'   "/genus/abditibacterium",
  #'   "/family/abditibacteriaceae",
  #'   "/order/abditibacteriales",
  #'   "/class/abditibacteriia",
  #'   "/phylum/abditibacteriota",
  #'   "/domain/bacteria"
  #' )
  #'
  #' # Desired ranks
  #' ranks <- c("Genus", "Family", "Order", "Class", "Phylum", "Domain")
  #'
  #' # Extract taxonomy
  #' extract_phylogeny(phylogeny_example, ranks)
  #'
  #' @export
  extract_phylogeny <- function(phylogeny, ranks) {
    # Initialize named vector with NA for all ranks
    named_phylogeny <- setNames(rep(NA, length(ranks)), ranks)
    
    # Loop through the phylogeny elements and match to the ranks
    for (item in phylogeny) {
      # Extract the rank and name (assuming format "/rank/name")
      rank <- stringr::str_extract(item, "(?<=/)[^/]+")  # Get the rank
      name <- stringr::str_extract(item, "[^/]+$")      # Get the name
      
      # If the rank matches one in our list, assign the name to it
      if (!is.na(rank) && rank %in% tolower(ranks)) {
        named_phylogeny[which(tolower(ranks) == rank)] <- stringr::str_to_title(name)  # Convert name to title case
      }
    }
    
    return(named_phylogeny)
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

# === Retrieve phylogeny of type strains ===
  addresses <- species_data$address
  base_url <- "https://lpsn.dsmz.de"
  phylogeny = vector("list", length(addresses))
  
  for (i in seq_along(addresses)) {
    # Visit page of type strain
      # Get body of the page
        url = addresses[i]
        body <- get_web_page_body(url = url, user_agent = "me")
           
      # Extract link after "Parent taxon:"
          parent_taxon_link <- body %>%
            rvest::html_nodes(xpath = "//b[contains(text(), 'Parent taxon:')]/following-sibling::a[1]") %>%
            rvest::html_attr("href")
      
      # Get first parent taxon
      parent_taxonomy = c(parent_taxon_link)
          
    # Visit pages of parent taxa
    # Repeat until there are no parent taxa left
    repeat {
        # Stop if no parent taxa are found
        if (is.na(parent_taxon_link) || length(parent_taxon_link) == 0) {
          break
        }
      
        # Get url of parent taxon
        url <- paste0(base_url, parent_taxon_link)
        
        # Get body of the page
        body <- get_web_page_body(url = url, user_agent = "me")
        
        # Extract the parent taxon link
        parent_taxon_link <- body %>%
          rvest::html_nodes(xpath = "//b[contains(text(), 'Parent taxon:')]/following-sibling::a[1]") %>%
          rvest::html_attr("href")
        
        # Stop if no parent taxa are found
        if (is.na(parent_taxon_link) || length(parent_taxon_link) == 0) {
          break
        }
        
        # Add the parent taxon link to the taxonomy vector
        parent_taxonomy <- c(parent_taxonomy, parent_taxon_link)
      }
  
      # Store phylogeny
      phylogeny[[i]] <- parent_taxonomy
  
      # Show progress of loop
      svMisc::progress(value = i, max.value = length(addresses))
  }
  
  # Extract phylogeny
  phylogeny_df <- purrr::map_dfr(
    phylogeny,
    ~ extract_phylogeny(.x, c("Domain", "Phylum", "Class", "Order", "Family", "Genus"))
  )
  
  # Combine phylogeny with species data
  species_data <- cbind(species_data, phylogeny_df)
  
# === Export ===
  fp <- paste0(base_path, "data\\LPSN_phylogeny.csv")
  write.csv(species_data, file = fp, row.names = FALSE)
