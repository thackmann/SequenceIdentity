#include <Rcpp.h>
using namespace Rcpp;

// Function to get most recent common ancestor (MRCA) for two tips in phylogenetic tree
// [[Rcpp::export]]
int get_mrca_cpp(IntegerMatrix edge, int Ntip, IntegerVector tip_indices) {
  int Nnode = edge.nrow() - Ntip;  // Number of internal nodes
  IntegerVector parents(Ntip + Nnode + 1, -1);  // Node-to-parent mapping
  
  // Convert tip_indices to 1-based indices
  for (int i = 0; i < tip_indices.size(); ++i) {
    tip_indices[i] += 1;
  }
  
  // Build the parent mapping
  for (int i = 0; i < edge.nrow(); i++) {
    parents[edge(i, 1)] = edge(i, 0);
  }
  
  // Find ancestors of the first tip
  std::unordered_set<int> ancestors;
  int current = tip_indices[0];
  while (current != -1) {
    ancestors.insert(current);
    current = parents[current];
  }
  
  // Traverse the second tip to find the first common ancestor
  current = tip_indices[1];
  while (ancestors.find(current) == ancestors.end()) {
    current = parents[current];
    if (current == -1) {
      stop("Failed to find MRCA: Reached root without finding common ancestor.");
    }
  }
  
  return current;
}

// Function to calculate RED for common ancestor of two tips
// [[Rcpp::export]]
double calculate_red_pair_cpp(IntegerMatrix edge, NumericVector reds, int Ntip, IntegerVector tip_indices) {
  try {
    // Get the MRCA
    int mrca = get_mrca_cpp(edge, Ntip, tip_indices);
    
    // Validate that the MRCA is an internal node
    if (mrca <= Ntip) {
      Rcpp::Rcerr << "Error: MRCA is a tip. Problematic pair: "
                  << tip_indices[0] << " (tip 1), "
                  << tip_indices[1] << " (tip 2), "
                  << "MRCA: " << mrca << std::endl;
      stop("Invalid MRCA: Node is a tip, not an internal node.");
    }
    
    // Fetch the RED value for the MRCA
    int red_index = mrca - Ntip - 1;  // Adjust MRCA index for 0-based REDs
    if (red_index < 0 || red_index >= reds.size()) {
      Rcpp::Rcerr << "Error: Invalid RED index for MRCA. Problematic pair: "
                  << tip_indices[0] << " (tip 1), "
                  << tip_indices[1] << " (tip 2), "
                  << "MRCA: " << mrca << ", RED index: " << red_index << std::endl;
      stop("Invalid RED index: Out of bounds.");
    }
    
    return reds[red_index];
  } catch (std::exception &e) {
    Rcpp::Rcerr << "Exception caught for tip pair: "
                << tip_indices[0] << " (tip 1), "
                << tip_indices[1] << " (tip 2). "
                << "Message: " << e.what() << std::endl;
    throw;
  }
}


// Function to calculate RED over a set of tips
// [[Rcpp::export]]
DataFrame calculate_red_set_cpp(IntegerMatrix edge, NumericVector reds, int Ntip, IntegerVector tip_indices, Rcpp::Environment pb_env) {
  // Ensure RED values are provided
  if (is_true(all(is_na(reds)))) {
    stop("RED values must be precomputed before calling this function.");
  }
  
  // Precompute all combinations of Tip1 and Tip2
  int total_combinations = tip_indices.size() * (tip_indices.size() - 1) / 2;
  std::vector<int> tip1_indices;
  std::vector<int> tip2_indices;
  NumericVector red_vector(total_combinations);
  
  // Access progress bar tick function
  Rcpp::Function pb_tick = pb_env["tick"];
  
  // Iterate through all pair combinations
  int index = 0;
  for (int i = 0; i < tip_indices.size() - 1; ++i) {
    for (int j = i + 1; j < tip_indices.size(); ++j) {
      // Extract tip indices for the pair
      IntegerVector pair = IntegerVector::create(tip_indices[i], tip_indices[j]); 
      
      // Calculate RED for the common ancestor
      double red_value = calculate_red_pair_cpp(edge, reds, Ntip, pair);
      
      // Store results
      tip1_indices.push_back(tip_indices[i]);
      tip2_indices.push_back(tip_indices[j]);
      red_vector[index++] = red_value;
      
      // Update progress bar
      pb_tick();
    }
  }
  
  // Return the resulting dataframe
  return DataFrame::create(
    Named("Tip1") = wrap(tip1_indices),
    Named("Tip2") = wrap(tip2_indices),
    Named("RED") = wrap(red_vector)
  );
}

// Function to calculate RED in chunks
// [[Rcpp::export]]
DataFrame calculate_red_chunk_cpp(IntegerMatrix edge, NumericVector reds, int Ntip, DataFrame chunk) {
  // Extract necessary information from the chunk
  IntegerVector var1_indices = chunk["Var1"];
  IntegerVector var2_indices = chunk["Var2"];
  
  // Ensure RED values are provided
  if (is_true(all(is_na(reds)))) {
    stop("RED values must be precomputed before calling this function.");
  }
  
  // Vectors to store results
  std::vector<int> tip1_indices;
  std::vector<int> tip2_indices;
  std::vector<double> red_vector;
  
  // Iterate through each pair in the chunk
  for (int k = 0; k < var1_indices.size(); ++k) {
    int i = var1_indices[k];  // Tip index for the first tip in the pair
    int j = var2_indices[k];  // Tip index for the second tip in the pair
    
    // Create the pair
    IntegerVector pair = IntegerVector::create(i, j);
    
    // Calculate RED for the pair
    double red_value = calculate_red_pair_cpp(edge, reds, Ntip, pair);
    
    // Store results
    tip1_indices.push_back(i);
    tip2_indices.push_back(j);
    red_vector.push_back(red_value);
  }
  
  // Return the resulting dataframe
  return DataFrame::create(
    Named("Tip1") = wrap(tip1_indices),
    Named("Tip2") = wrap(tip2_indices),
    Named("RED") = wrap(red_vector)
  );
}