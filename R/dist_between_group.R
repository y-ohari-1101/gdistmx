dist_between_group <- function(dna, group_df, model = "raw") {
  
  # ================================================================
  # 0. Recursive processing for list input
  # ================================================================
  if (is.list(dna) && all(sapply(dna, inherits, "DNAbin"))) {
    
    if (is.null(names(dna))) {
      stop("dist_between_group(): list input must have names for each gene.")
    }
    
    out <- lapply(names(dna), function(gene) {
      res <- dist_between_group(dna[[gene]], group_df, model=model)
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
    stop("dist_between_group(): 'dna' must be a DNAbin object.")
  if (is.null(rownames(dna)))
    stop("dist_between_group(): DNAbin must contain rownames (sample IDs).")
  
  # ================================================================
  # 2. Validate group_df
  # ================================================================
  if (!is.data.frame(group_df)) stop("group_df must be a data.frame.")
  if (ncol(group_df) < 2) stop("group_df must contain ≥2 columns.")
  
  colnames(group_df)[1:2] <- c("sample_id", "group")
  group_df$sample_id <- as.character(group_df$sample_id)
  group_df$group     <- as.character(group_df$group)
  
  if (length(unique(group_df$group)) < 2) {
    stop("dist_between_group(): at least two groups are required.")
  }
  
  # ================================================================
  # 3. Duplicate checks
  # ================================================================
  if (any(duplicated(rownames(dna))))
    stop("Duplicate sample names found in DNAbin.")
  if (any(duplicated(group_df$sample_id)))
    stop("Duplicate sample IDs found in group_df.")
  
  # ================================================================
  # 4. Match samples
  # ================================================================
  group_df <- group_df[group_df$sample_id %in% rownames(dna), ]
  if (nrow(group_df) == 0)
    stop("dist_between_group(): no sample IDs match.")
  
  # ================================================================
  # 5. Distance matrix
  # ================================================================
  dist_mat <- as.matrix(ape::dist.dna(dna, model=model, pairwise.deletion=TRUE))
  
  # ================================================================
  # 6. Between-group distances
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
      vals <- round(dist_mat[cbind(comb$sample1, comb$sample2)], 4)
      
      pair_list[[paste(g1, g2, sep="_")]] <- data.frame(
        group1=g1, group2=g2,
        sample1=comb$sample1, sample2=comb$sample2,
        distance=vals
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
