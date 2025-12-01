#' @title dist_within_group
#' @description \code{dist_within_group} Calculating genetic distances within each group
#' @param dna DNAbin inported by ape::read.dna()
#' @param group_df group definition file
#' @param model genetic distance model available in ape::dist.dna()
#' @return List consist of genetic distance among individuals and groups
#' @export
#' @import utils
#' @seealso ape


dist_within_group <- function(dna, group_df, model = "raw") {

  # ================================================================
  # 0. Recursive processing for list input
  # ================================================================
  if (is.list(dna) && all(sapply(dna, inherits, "DNAbin"))) {

    if (is.null(names(dna))) {
      stop("dist_within_group(): list input must have names for each gene.")
    }

    out <- lapply(names(dna), function(gene) {
      res <- dist_within_group(dna[[gene]], group_df, model = model)
      res$gene <- gene
      return(res)
    })

    names(out) <- names(dna)
    return(out)
  }

  # ================================================================
  # 1. Validate DNAbin input
  # ================================================================
  if (!inherits(dna, "DNAbin")) {
    stop("dist_within_group(): 'dna' must be a DNAbin object.")
  }
  if (is.null(rownames(dna))) {
    stop("dist_within_group(): DNAbin object must contain rownames (sample IDs).")
  }

  # ================================================================
  # 2. Validate group_df structure
  # ================================================================
  if (!is.data.frame(group_df))
    stop("dist_within_group(): group_df must be a data.frame.")
  if (ncol(group_df) < 2)
    stop("dist_within_group(): group_df must contain >=2 columns (sample_id, group).")

  colnames(group_df)[1:2] <- c("sample_id", "group")
  group_df$sample_id <- as.character(group_df$sample_id)
  group_df$group     <- as.character(group_df$group)

  # ================================================================
  # 3. Duplicate checks
  # ================================================================
  if (any(duplicated(rownames(dna)))) {
    stop("dist_within_group(): duplicate sample names found in DNAbin rownames.")
  }
  if (any(duplicated(group_df$sample_id))) {
    stop("dist_within_group(): duplicate sample IDs found in group_df.")
  }

  # ================================================================
  # 4. Match group_df samples to DNAbin
  # ================================================================
  group_df <- group_df[group_df$sample_id %in% rownames(dna), ]
  if (nrow(group_df) == 0) {
    warning("dist_within_group(): no sample IDs in group_df match DNAbin rownames.")
    return(list(
      summary  = data.frame(group=character(), min=numeric(), mean=numeric(), max=numeric()),
      pairwise = data.frame(group=character(), sample1=character(), sample2=character(), distance=numeric())
    ))
  }

  # ================================================================
  # 5. Compute distance matrix
  # ================================================================
  dist_mat <- as.matrix(ape::dist.dna(dna, model=model, pairwise.deletion=TRUE))

  # ================================================================
  # 6. Compute within-group distances
  # ================================================================
  groups <- unique(group_df$group)
  result_list <- list()
  pair_list   <- list()

  for (grp in groups) {

    samples <- group_df$sample_id[group_df$group == grp]

    if (length(samples) < 2) {
      result_list[[grp]] <- data.frame(group=grp, min=0, mean=0, max=0)
      next
    }

    dmat <- dist_mat[samples, samples, drop=FALSE]
    comb <- t(combn(samples, 2))
    vals <- round(dmat[comb], 4)

    pair_list[[grp]] <- data.frame(
      group=grp,
      sample1=comb[,1],
      sample2=comb[,2],
      distance=vals
    )

    result_list[[grp]] <- data.frame(
      group=grp,
      min=min(vals),
      mean=mean(vals),
      max=max(vals)
    )
  }

  return(list(
    summary  = dplyr::bind_rows(result_list),
    pairwise = dplyr::bind_rows(pair_list)
  ))
}
