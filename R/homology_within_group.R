#' @title homology_within_group
#' @description \code{dist_within_group} Calculating genetic distances within each group
#' @param dna DNAbin inported by ape::read.dna()
#' @param group_df group definition file
#' @return List consist of genetic distance among individuals and groups
#' @export
#' @import utils
#' @seealso ape

homology_within_group <- function(dna, group_df) {

  # ================================================================
  # 0. Recursive processing for list input
  # ================================================================
  if (is.list(dna) && all(sapply(dna, inherits, "DNAbin"))) {

    if (is.null(names(dna))) stop("homology_within_group(): list must have names.")

    out <- lapply(names(dna), function(gene) {
      res <- homology_within_group(dna[[gene]], group_df)
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
    stop("homology_within_group(): 'dna' must be a DNAbin object.")
  if (is.null(rownames(dna)))
    stop("homology_within_group(): DNAbin must have rownames.")

  # ================================================================
  # 2. Validate group_df
  # ================================================================
  if (!is.data.frame(group_df)) stop("group_df must be a data.frame.")
  if (ncol(group_df) < 2) stop("group_df must contain >=2 columns.")

  colnames(group_df)[1:2] <- c("sample_id", "group")
  group_df$sample_id <- as.character(group_df$sample_id)
  group_df$group     <- as.character(group_df$group)

  # ================================================================
  # 3. Duplicate checks
  # ================================================================
  if (any(duplicated(rownames(dna)))) stop("Duplicate sample names in DNAbin.")
  if (any(duplicated(group_df$sample_id))) stop("Duplicate sample IDs in group_df.")

  # ================================================================
  # 4. Match samples
  # ================================================================
  group_df <- group_df[group_df$sample_id %in% rownames(dna), ]
  if (nrow(group_df) == 0) stop("No matched samples found.")

  # ================================================================
  # 5. Compute homology matrix
  # ================================================================
  dist_mat <- as.matrix(ape::dist.dna(dna, model="raw", pairwise.deletion=TRUE))
  homology_mat <- 100 * (1 - dist_mat)

  # ================================================================
  # 6. Compute within-group homology
  # ================================================================
  groups <- unique(group_df$group)
  result_list <- list()
  pair_list   <- list()

  for (grp in groups) {
    samples <- group_df$sample_id[group_df$group == grp]

    if (length(samples) < 2) {
      result_list[[grp]] <- data.frame(group=grp, min=100, mean=100, max=100)
      next
    }

    hmat <- homology_mat[samples, samples, drop=FALSE]
    comb <- t(combn(samples, 2))
    vals <- round(hmat[comb], 4)

    pair_list[[grp]] <- data.frame(
      group=grp,
      sample1=comb[,1], sample2=comb[,2],
      homology=vals
    )

    result_list[[grp]] <- data.frame(
      group=grp,
      min=min(vals), mean=mean(vals), max=max(vals)
    )
  }

  return(list(
    summary  = dplyr::bind_rows(result_list),
    pairwise = dplyr::bind_rows(pair_list)
  ))
}
