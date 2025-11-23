homology_between_group <- function(dna, group_df) {
  
  # ================================================================
  # 0. Recursive list processing
  # ================================================================
  if (is.list(dna) && all(sapply(dna, inherits, "DNAbin"))) {
    
    if (is.null(names(dna))) stop("List input must be named.")
    
    out <- lapply(names(dna), function(gene) {
      res <- homology_between_group(dna[[gene]], group_df)
      res$gene <- gene
      return(res)
    })
    
    names(out) <- names(dna)
    return(out)
  }
  
  # ================================================================
  # 1. Validate DNAbin
  # ================================================================
  if (!inherits(dna, "DNAbin"))
    stop("homology_between_group(): 'dna' must be DNAbin.")
  if (is.null(rownames(dna)))
    stop("homology_between_group(): DNAbin requires rownames.")
  
  # ================================================================
  # 2. Validate group_df
  # ================================================================
  if (!is.data.frame(group_df)) stop("group_df must be a data.frame.")
  if (ncol(group_df) < 2) stop("group_df must contain ≥2 columns.")
  
  colnames(group_df)[1:2] <- c("sample_id", "group")
  group_df$sample_id <- as.character(group_df$sample_id)
  group_df$group     <- as.character(group_df$group)
  
  if (length(unique(group_df$group)) < 2)
    stop("homology_between_group(): at least two groups required.")
  
  # ================================================================
  # 3. Duplicate checks
  # ================================================================
  if (any(duplicated(rownames(dna))))
    stop("Duplicate sample names in DNAbin.")
  if (any(duplicated(group_df$sample_id)))
    stop("Duplicate sample IDs in group_df.")
  
  # ================================================================
  # 4. Matching
  # ================================================================
  group_df <- group_df[group_df$sample_id %in% rownames(dna), ]
  if (nrow(group_df) == 0) stop("No matching sample IDs.")
  
  # ================================================================
  # 5. Compute homology matrix
  # ================================================================
  dist_mat <- as.matrix(ape::dist.dna(dna, model="raw", pairwise.deletion=TRUE))
  homology_mat <- 100 * (1 - dist_mat)
  
  # ================================================================
  # 6. Compute between-group homology
  # ================================================================
  groups <- unique(group_df$group)
  result_list <- list()
  pair_list   <- list()
  
  for (i in seq_along(groups)) {
    for (j in seq_along(groups)) if (j > i) {
      
      g1 <- groups[i]
      g2 <- groups[j]
      
      s1 <- group_df$sample_id[group_df$group == g1]
      s2 <- group_df$sample_id[group_df$group == g2]
      
      comb <- expand.grid(sample1=s1, sample2=s2, stringsAsFactors=FALSE)
      vals <- round(homology_mat[cbind(comb$sample1, comb$sample2)], 4)
      
      pair_list[[paste(g1, g2, sep="_")]] <- data.frame(
        group1=g1, group2=g2,
        sample1=comb$sample1, sample2=comb$sample2,
        homologylogy=vals
      )
      
      result_list[[paste(g1, g2, sep="_")]] <- data.frame(
        group1=g1, group2=g2,
        min=min(vals), mean=mean(vals), max=max(vals)
      )
    }
  }
  
  return(list(
    summary  = dplyr::bind_rows(result_list),
    pairwise = dplyr::bind_rows(pair_list)
  ))
}
