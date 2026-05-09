#' @importFrom utils read.delim
# Internal helper function to read BGI GO enrichment result files
#' @importFrom readxl read_excel
read_bgi_go <- function(file_path) {
  raw <- tryCatch(
    suppressMessages(readxl::read_excel(file_path)),
    error = function(e) {
      tryCatch(
        read.delim(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                   quote = "", comment.char = ""),
        error = function(e2) {
          stop("Unable to read file: ", file_path)
        }
      )
    }
  )
  raw <- as.data.frame(raw)

  # Map BGI column names to standard names used for plotting
  if ("go_p_term_id" %in% colnames(raw)) {
    raw$GO.ID <- raw$go_p_term_id
  }
  if ("go_p_term_desc" %in% colnames(raw)) {
    raw$Term <- raw$go_p_term_desc
  }
  if ("go_p_pvalue" %in% colnames(raw)) {
    raw$P.value <- as.numeric(raw$go_p_pvalue)
  }
  if ("gene_symbols" %in% colnames(raw)) {
    raw$Genes <- raw$gene_symbols
  }

  keep_cols <- c("GO.ID", "Term", "P.value", "Genes")
  available <- intersect(keep_cols, colnames(raw))
  raw[, available, drop = FALSE]
}

#' Plot a comparison of BGI GO enrichment results between two stress models
#'
#' Creates a comparative bubble plot showing the top 10 enriched GO terms
#' from the reserpine and chronic stress models, using BGI-generated GO
#' enrichment results.
#'
#' @param ontology Character. GO sub-ontology: "BP", "CC", or "MF". Default "BP".
#' @param reserpine_file Character. Path to the reserpine model GO result file.
#' @param cs_file Character. Path to the chronic stress model GO result file.
#' @param title Character. Title for the plot.
#'
#' @return A ggplot2 object (printed to the active graphics device).
#' @export
#' @importFrom ggplot2 ggplot aes geom_point scale_size scale_color_gradient labs theme_minimal theme element_text
#'
#' @examples
#' \dontrun{
#' plot_go_comparison("BP", "BGI/GO/reserpine_BP.txt", "BGI/GO/cs_BP.txt")
#' }
plot_go_comparison <- function(ontology = "BP",
                               reserpine_file = NULL,
                               cs_file = NULL,
                               title = "GO Enrichment Comparison: Reserpine vs Chronic Stress") {

  # Use default paths to bundled BGI data if not specified
  if (is.null(reserpine_file)) {
    reserpine_file <- system.file("extdata/BGI/Reserpine_GO/GO_enrichment/Control-vs-Reserpine/Control-vs-Reserpine_gene_go_P_rich_term.xls",
                                  package = "TLR.zebrafish.stress")
  }
  if (is.null(cs_file)) {
    cs_file <- system.file("extdata/BGI/CS_GO/GO_enrichment/Control-vs-Stress/Control-vs-Stress_gene_go_P_rich_term.xls",
                           package = "TLR.zebrafish.stress")
  }

  res_go <- read_bgi_go(reserpine_file)
  cs_go  <- read_bgi_go(cs_file)

  # Label each dataset with its stress model
  res_go$Model <- "Reserpine"
  cs_go$Model  <- "Chronic Stress"

  # Merge the two GO datasets
  combined <- rbind(res_go, cs_go)
  # Select top 10 most significant terms for display
  combined <- combined[order(combined$P.value), ]
  top_terms <- unique(combined$GO.ID)[1:10]
  plot_data <- combined[combined$GO.ID %in% top_terms, ]

  ggplot2::ggplot(plot_data, ggplot2::aes(x = Model, y = Term, size = -log10(P.value), color = -log10(P.value))) +
    ggplot2::geom_point() +
    ggplot2::scale_size(range = c(3, 10)) +
    ggplot2::scale_color_gradient(low = "blue", high = "red") +
    ggplot2::labs(title = title, x = "Stress Model", y = "GO Term", size = "-log10(p)", color = "-log10(p)") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 10),
                   legend.position = "right")
}
