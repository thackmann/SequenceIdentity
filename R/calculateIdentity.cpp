// [[Rcpp::depends(RSeqAn)]]
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/align.h>
#include <seqan/score.h> // Include for scoring schemes
#include <Rcpp.h>

using namespace seqan;
using namespace Rcpp;

// Function to convert an alignment row to a std::string
std::string extractAlignedRow(seqan::Row<seqan::Align<seqan::String<char>, seqan::ArrayGaps>>::Type &row) {
  std::string alignedSequence;
  for (size_t i = 0; i < length(row); ++i) {
    alignedSequence += isGap(row, i) ? '-' : row[i];
  }
  return alignedSequence;
}

// Function to perform pairwise alignment of two sequences
// [[Rcpp::export]]
List perform_alignment(std::string seq1, std::string seq2, 
                       int match = 2, int mismatch = -1, 
                       int gap_open = -10, int gap_extend = -0.5) {
  // Convert input sequences to SeqAn string type
  seqan::String<char> sequence1 = seq1.c_str();
  seqan::String<char> sequence2 = seq2.c_str();
  
  // Create the alignment object and initialize rows
  seqan::Align<seqan::String<char>, seqan::ArrayGaps> align;
  resize(rows(align), 2);
  assignSource(row(align, 0), sequence1);
  assignSource(row(align, 1), sequence2);
  
  // Define the scoring scheme with AffineGaps() tag
  seqan::Score<int> scoringScheme(match, mismatch, gap_extend, gap_open);
  
  // Perform the global alignment using the AffineGaps() tag
  int alignmentScore = globalAlignment(align, scoringScheme, AffineGaps());
  
  // Extract the aligned sequences
  std::string alignedRow1 = extractAlignedRow(row(align, 0));
  std::string alignedRow2 = extractAlignedRow(row(align, 1));
  
  // Convert the alignment object to a string for output
  std::ostringstream alignmentOutput;
  alignmentOutput << align;
  
  return List::create(
    _["alignment_score"] = alignmentScore,
    _["alignment"] = alignmentOutput.str(),
    _["aligned_row1"] = alignedRow1,
    _["aligned_row2"] = alignedRow2
  );
}

// Function to calculate percent identity from aligned sequences
// [[Rcpp::export]]
List calculate_identity_from_alignment(std::string aligned_row1, std::string aligned_row2, bool ignore_gap = true) {
  int matchCount = 0, mismatchCount = 0, totalAligned = 0;
  
  // Trim leading and trailing gaps
  size_t start = 0, end = aligned_row1.size();
  while (start < end && (aligned_row1[start] == '-' || aligned_row2[start] == '-')) ++start;
  while (end > start && (aligned_row1[end - 1] == '-' || aligned_row2[end - 1] == '-')) --end;
  
  // Calculate matches, mismatches, and total aligned positions
  for (size_t i = start; i < end; ++i) {
    if (ignore_gap) {
      // Ignore positions where either sequence has a gap
      if (aligned_row1[i] != '-' && aligned_row2[i] != '-') {
        ++totalAligned;  // Count only positions where neither sequence has a gap
        if (aligned_row1[i] == aligned_row2[i]) {
          ++matchCount;
        } else {
          ++mismatchCount;
        }
      }
    } else {
      // Count all positions, including gaps
      if (aligned_row1[i] == aligned_row2[i]) {
        ++matchCount;  // Count as a match if characters are the same
      } else {
        ++mismatchCount;  // Count as a mismatch if characters differ or one is a gap
      }
      ++totalAligned;  // Count every position
    }
  }
  
  // Calculate percent identity
  double percentIdentity = (totalAligned > 0) ? (100.0 * matchCount / totalAligned) : 0.0;
  
  return List::create(
    _["percent_identity"] = percentIdentity,
    _["match_count"] = matchCount,
    _["mismatch_count"] = mismatchCount,
    _["totalAligned"] = totalAligned
  );
}

// Function to calculate identity over set of sequences
// [[Rcpp::export]]
DataFrame calculate_identity_set_cpp(CharacterVector sequences, Rcpp::Environment pb_env) {
  int n = sequences.size();
  
  // Vectors to store results
  std::vector<int> var1_indices;
  std::vector<int> var2_indices;
  std::vector<double> freq;
  
  // Access progress bar tick function
  Rcpp::Function pb_tick = pb_env["tick"];
  
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {  // Exclude diagonal
      // Perform alignment
      std::string seq1 = Rcpp::as<std::string>(sequences[i]);
      std::string seq2 = Rcpp::as<std::string>(sequences[j]);
      List alignmentResult = perform_alignment(seq1, seq2, 2, -1, -10, -0.5);
      
      // Calculate identity
      std::string alignedRow1 = Rcpp::as<std::string>(alignmentResult["aligned_row1"]);
      std::string alignedRow2 = Rcpp::as<std::string>(alignmentResult["aligned_row2"]);
      List identityResult = calculate_identity_from_alignment(alignedRow1, alignedRow2);
      double percentIdentity = Rcpp::as<double>(identityResult["percent_identity"]);
      
      // Append indices to results
      var1_indices.push_back(i + 1);  // Use 1-based indexing for R compatibility
      var2_indices.push_back(j + 1);  // Use 1-based indexing for R compatibility
      freq.push_back(percentIdentity);
      
      // Update progress bar
      pb_tick();
    }
  }
  
  // Return results as a DataFrame
  return DataFrame::create(
    _["Var1"] = var1_indices,
    _["Var2"] = var2_indices,
    _["Freq"] = freq
  );
}

// Function to calculate identity in chunks
// [[Rcpp::export]]
DataFrame calculate_identity_chunk_cpp(CharacterVector sequences, DataFrame chunk) {
  // Extract Var1 and Var2 indices from the chunk
  IntegerVector var1_indices = chunk["Var1"];
  IntegerVector var2_indices = chunk["Var2"];
  
  // Vectors to store results
  std::vector<int> var1_results;
  std::vector<int> var2_results;
  std::vector<double> freq;

  for (int k = 0; k < var1_indices.size(); ++k) {
    int i = var1_indices[k] - 1;  // Convert to 0-based indexing
    int j = var2_indices[k] - 1;  // Convert to 0-based indexing
    
    // Perform alignment
    std::string seq1 = Rcpp::as<std::string>(sequences[i]);
    std::string seq2 = Rcpp::as<std::string>(sequences[j]);
    List alignmentResult = perform_alignment(seq1, seq2, 2, -1, -10, -0.5);
    
    // Calculate identity
    std::string alignedRow1 = Rcpp::as<std::string>(alignmentResult["aligned_row1"]);
    std::string alignedRow2 = Rcpp::as<std::string>(alignmentResult["aligned_row2"]);
    List identityResult = calculate_identity_from_alignment(alignedRow1, alignedRow2);
    double percentIdentity = Rcpp::as<double>(identityResult["percent_identity"]);
    
    // Append results
    var1_results.push_back(i + 1);  // Store 1-based indices
    var2_results.push_back(j + 1);  // Store 1-based indices
    freq.push_back(percentIdentity);
  }

  // Return results as a DataFrame
  return DataFrame::create(
    _["Var1"] = var1_results,
    _["Var2"] = var2_results,
    _["Freq"] = freq
  );
}