#' @title export_matrix_csv
#' @description \code{dist_within_group} Calculating genetic distances within each group
#' @return List consist of genetic distance among individuals and groups
#' @param result test
#' @param out_dir test
#' @export
#' @seealso ape


export_matrix_csv <- function(result, out_dir = "./res") {

  # ================================================================
  # 1. Output directory
  # ================================================================
  if (is.null(out_dir) || out_dir == "") out_dir <- "."
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)


  # ================================================================
  # 2. Case 1: multiple-gene list (named list)
  # ================================================================
  if (is.list(result) && !("matrix" %in% names(result))) {

    invisible(lapply(names(result), function(gene_name) {

      df <- result[[gene_name]]$matrix
      if (is.null(df)) return(NULL)

      # safety: fallback name
      if (is.null(gene_name) || gene_name == "") gene_name <- "Gene"

      out_path <- file.path(out_dir, paste0(gene_name, "_matrix.csv"))
      write.csv(df, out_path, row.names = TRUE)

      message("Saved: ", out_path)
      NULL
    }))
    return(invisible(NULL))
  }


  # ================================================================
  # 3. Case 2: single-gene structure (one distance_mx() / homology_mx())
  # ================================================================
  if ("matrix" %in% names(result)) {

    gene <- if (!is.null(result$gene)) result$gene else "Gene"

    out_path <- file.path(out_dir, paste0(gene, "_matrix.csv"))
    write.csv(result$matrix, out_path, row.names = TRUE)

    message("Saved: ", out_path)
    return(invisible(NULL))
  }


  # ================================================================
  # 4. Otherwise → error
  # ================================================================
  stop("export_matrix_csv(): input is not a valid result from distance_mx() or homology_mx().")
}
