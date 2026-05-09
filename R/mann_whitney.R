#' Perform Mann-Whitney U Test for Gene Expression Data
#'
#' Compares treatment versus control groups for each gene using the
#' two-sided Mann-Whitney U test, which does not assume normality
#' and is suitable for small sample sizes.
#'
#' @param data A data frame containing expression values and group labels.
#' @param gene_list Character vector of gene column names to test.
#' @param group_col Character. Name of the column containing group labels
#'   (must include "control group" and "experimental group").
#' @param model_name Character. Name of the stress model for output labeling.
#'
#' @return A data frame with columns: Gene, P_value, Significant.
#' @export
#' @importFrom stats na.omit wilcox.test
#'
#' @examples
#' \dontrun{
#' data <- read.csv("expression_data.csv")
#' genes <- c("TLR4", "MyD88", "FADD")
#' result <- run_mann_whitney(data, genes, "Categories", "Reserpine")
#' }
run_mann_whitney <- function(data, gene_list, group_col, model_name) {

  results <- data.frame(
    Gene = character(),
    P_value = numeric(),
    Significant = character(),
    stringsAsFactors = FALSE
  )

  # Split data into control and treatment groups
  control_data <- data[data[[group_col]] == "control group", ]
  treat_data   <- data[data[[group_col]] == "experimental group", ]

  for (gene in gene_list) {

    # Extract and clean expression values
    ctrl_vals  <- as.numeric(stats::na.omit(control_data[[gene]]))
    treat_vals <- as.numeric(stats::na.omit(treat_data[[gene]]))

    if (length(ctrl_vals) < 3 | length(treat_vals) < 3) {
      results <- rbind(results, data.frame(
        Gene = gene, P_value = NA, Significant = "Insufficient data"
      ))
      next
    }

    # Run Mann-Whitney U test
    test_res <- stats::wilcox.test(treat_vals, ctrl_vals, exact = FALSE)
    p_val    <- round(test_res$p.value, 4)
    sig_lab  <- ifelse(p_val < 0.05, "Yes (*)", "No")

    results <- rbind(results, data.frame(
      Gene = gene, P_value = p_val, Significant = sig_lab
    ))
  }

  cat("\n==========", model_name, "==========\n")
  print(results)
  invisible(results)
}
