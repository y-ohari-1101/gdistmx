distance_mx <- function(dna, group_df, model = "raw") {
  
  # ================================================================
  # 0. Recursive handling for list input
  # ================================================================
  if (is.list(dna) && all(sapply(dna, inherits, "DNAbin"))) {
    
    if (is.null(names(dna))) {
      stop("distance_mx(): list input must have names for each gene.")
    }
    
    out <- lapply(names(dna), function(gene) {
      res <- distance_mx(dna[[gene]], group_df, model=model)
      res$gene <- gene
      return(res)
    })
    names(out) <- names(dna)
    return(out)
  }
  
  # ================================================================
  # 1. Validate DNAbin input
  # ================================================================
  if (!inherits(dna, "DNAbin"))
    stop("distance_mx(): dna must be a DNAbin object.")
  
  if (is.null(rownames(dna)))
    stop("distance_mx(): DNAbin requires rownames (sample IDs).")
  
  # ================================================================
  # 2. Compute within- and between-group distance
  # ================================================================
  within  <- dist_within_group(dna, group_df, model=model)
  between <- dist_between_group(dna, group_df, model=model)
  
  # ================================================================
  # 3. Construct distance matrix
  #    - diagonal: within (mean + range)
  #    - lower: between mean
  #    - upper: between min–max
  # ================================================================
  groups <- unique(group_df$group)
  G <- length(groups)
  mat <- matrix("", nrow=G, ncol=G, dimnames=list(groups, groups))
  
  # diagonal
  for (g in groups) {
    row <- within$summary[within$summary$group == g, ]
    if (nrow(row) == 0) {
      mat[g, g] <- ""
    } else {
      mm <- sprintf("%.4f (%.4f - %.4f)", row$mean, row$min, row$max)
      mat[g, g] <- mm
    }
  }
  
  # between
  if (!is.null(between$summary) && nrow(between$summary) > 0) {
    for (i in seq_along(groups)) {
      for (j in seq_along(groups)) if (j > i) {
        
        g1 <- groups[i]
        g2 <- groups[j]
        
        row <- between$summary[
          between$summary$group1 == g1 &
            between$summary$group2 == g2, ]
        
        if (nrow(row) == 0) next
        
        mean_str <- sprintf("%.4f", row$mean)
        range_str <- sprintf("%.4f - %.4f", row$min, row$max)
        
        mat[g2, g1] <- mean_str
        mat[g1, g2] <- range_str
      }
    }
  }
  
  # ================================================================
  # 4. Return
  # ================================================================
  return(list(
    within_summary   = within$summary,
    within_pairwise  = within$pairwise,
    between_summary  = between$summary,
    between_pairwise = between$pairwise,
    matrix           = as.data.frame(mat)
  ))
}
