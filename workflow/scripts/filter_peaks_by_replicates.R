# setup -----------------------------------------------------------------------
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(GenomicRanges))

# define functions -------------------------------------------------------------
# function to merge multiple sets of ChIP peaks into a non-overlapping master set
# supports multiple methods
# peaks should be a GRangesList
# method should be either "overlap" or "merge": 
# overlap method will determine overlaps among all peak sets and combine non-overlapping peaks
# merge method will combine overlapping peaks into a single peak
combine_peaks <- function(peaks, method = "overlap", min_overlap = 1L) {
  require(GenomicRanges)
  
  # check arguments
  if (length(peaks) < 2) {
    stop("Only one set of peaks provided. For merging peaks, provide multiple peak sets.")
  }
  
  if (!inherits(peaks, "GRangesList") ) {
    stop("one or more peak set is not a Granges object. peaks must be a GRangesList")
  }
  
  if (!any(method == c("overlap", "merge"))) {
    stop("Invalid method provided. method should be 'overlap' or 'merge' ")
  }
  
  if (!is_integer(min_overlap)) {
    stop("min_overlap must be an integer")
  }
  
  
  # merge peaks using "overlap" method
  if (method == "overlap") {
    combined_peaks <- peaks[[1]]
    for (i in 2:length(peaks)) {
      i_set <- peaks[[i]]
      combined_peaks <- c(combined_peaks,
                          subsetByOverlaps(i_set, combined_peaks, minoverlap = min_overlap, invert = TRUE))
    }
    
  }
  # merge peaks using "merge" method
  if (method == "merge") {
    combined_peaks <- GenomicRanges::reduce(unlist(peaks))
    
  }
  mcols(combined_peaks) <- NULL
  return(combined_peaks)
  
}

# function to build an overlap table
# takes a Granges list object
# combines all nonoverlapping peaks from all samples, then builds a logical table indicating which peaks from each sample overlap each peak on the master list
peak_overlap_table <- function(peaks, method = "overlap", min_overlap = 1L, combine_peaks = TRUE) {
  
  # check arguments
  if (length(peaks) < 2) {
    stop("Only one set of peaks provided. For merging peaks, provide multiple peak sets.")
  }
  
  if (!inherits(peaks, "GRangesList") ) {
    stop("one or more peak set is not a Granges object. peaks must be a GRangesList")
  }
  
  if (!any(method == c("overlap", "merge"))) {
    stop("Invalid method provided. method should be 'overlap' or 'merge' ")
  }
  
  if (!is_integer(min_overlap)) {
    stop("min_overlap must be an integer")
  }
  
  
  # set names
  if(is.null(names(peaks))) {
    names(peaks) <- paste0("sample_", seq((peaks)))
  }
  
  
  if (combine_peaks) {
    # get a master set of all nonoverlapping peaks from all peak sets
    all_peaks_gr <- combine_peaks(peaks, method = method, min_overlap = min_overlap)
    
  } else {
    all_peaks_gr <- peaks[[1]]
  }
  all_peaks_df <- as.data.frame(all_peaks_gr) %>%
    dplyr::select(1:5)
  
  # build a logical overlap table indicating which peaks were detected in which samples
  for (i in 1:length(peaks)) {
    sample_name <- names(peaks)[i]
    i_peaks <- peaks[[i]]
    
    all_peaks_df[,sample_name] <- FALSE
    overlaps <- findOverlaps(all_peaks_gr, i_peaks, minoverlap = min_overlap)@from
    all_peaks_df[overlaps,sample_name] <- TRUE
  }
  return(all_peaks_df)
}

# read peaks into overlap table ------------------------------------------------
merged_peak_file <- snakemake@input[["merged_peaks"]]
names(merged_peak_file) <- "merged_peaks"

replicate_peak_files <- snakemake@input[["ind_peaks"]]
names(replicate_peak_files) <- paste0("rep", seq(replicate_peak_files), "_peaks")

peak_files <- c(merged_peak_file, replicate_peak_files)

overlap_table <- peak_files %>% 
  map(rtracklayer::import) %>% 
  GRangesList() %>% 
  peak_overlap_table(min_overlap = as.integer(snakemake@params[["min_overlap"]]))

# filter to include only peaks detected in multiple replicates -----------------
filtered_peaks <- overlap_table %>% 
  mutate(reps_w_peak = rowSums(.[7:ncol(.)])) %>% 
  filter(merged_peaks, reps_w_peak >= snakemake@params[["replicate_threshold"]]) %>% 
  makeGRangesFromDataFrame()

# export filtered peaks to file  -----------------------------------------------
if (snakemake@wildcards[["peak_type"]] == "narrow") {
  keep_cols <- c(1:3, 6:7, 5, 8:11)
} else if (snakemake@wildcards[["peak_type"]] == "broad") {
  keep_cols <- c(1:3, 6:7, 5, 8:10)
} else {
  stop("Peak format should be 'narrow' or 'broad'. Check input file format")
}

output <- merged_peak_file %>% 
  rtracklayer::import() %>% 
  subsetByOverlaps(filtered_peaks) %>% 
  as.data.frame() %>% 
  select(all_of(keep_cols)) %>% 
  mutate(strand = ".")

write_tsv(output, snakemake@output[[1]], col_names = FALSE)
