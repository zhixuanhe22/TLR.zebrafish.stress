#' Create a bar plot comparing qPCR and RNA-seq results for TLR genes
#'
#' Generates a grouped bar chart that compares log2 fold changes measured
#' by qPCR with those from BGI RNA-seq for the 15 Toll-like receptor pathway
#' genes. Significant genes from the Mann-Whitney U test are marked with
#' asterisks.
#'
#' @param qpcr_data A data frame containing qPCR-derived statistics per gene,
#'   with columns Gene, log2FC, and Significant.
#' @param rnaseq_data A data frame containing RNA-seq-derived statistics per
#'   gene, with columns Gene and log2FC.
#' @param title Character. Title for the plot.
#'
#' @return A ggplot2 object (printed to the active graphics device).
#' @export
#' @importFrom ggplot2 ggplot aes geom_bar position_dodge labs theme_minimal theme element_text scale_fill_manual
#' @importFrom reshape2 melt
#'
#' @examples
#' \dontrun{
#' qpcr <- data.frame(Gene = c("TLR4","MyD88","FADD"),
#'                    log2FC = c(1.5, 2.1, 1.8),
#'                    Significant = c("Yes","Yes","Yes"))
#' rnaseq <- data.frame(Gene = c("TLR4","MyD88","FADD"),
#'                      log2FC = c(1.3, 1.9, 2.0))
#' plot_qpcr_vs_rnaseq(qpcr, rnaseq, "Comparison")
#' }
plot_qpcr_vs_rnaseq <- function(qpcr_data, rnaseq_data, title = "qPCR vs RNA-seq log2 Fold Change") {

  # Combine qPCR and RNA-seq datasets by gene
  plot_data <- merge(qpcr_data, rnaseq_data, by = "Gene", suffixes = c("_qPCR","_RNAseq"))

  # Reshape data into long format for ggplot
  plot_data_long <- reshape2::melt(plot_data, id.vars = c("Gene","Significant"),
                                   measure.vars = c("log2FC_qPCR","log2FC_RNAseq"),
                                   variable.name = "Platform", value.name = "log2FC")

  plot_data_long$Platform <- ifelse(plot_data_long$Platform == "log2FC_qPCR", "qPCR", "RNA-seq")
  plot_data_long$Significant <- ifelse(is.na(plot_data_long$Significant) | plot_data_long$Significant != "Yes (*)", "No", "Yes")

  # Build grouped bar plot
  p <- ggplot2::ggplot(plot_data_long, ggplot2::aes(x = Gene, y = log2FC, fill = Platform)) +
    ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge(width = 0.8), width = 0.7) +
    ggplot2::labs(title = title, x = "Gene", y = expression(log[2]~"Fold Change")) +
    ggplot2::scale_fill_manual(values = c("qPCR" = "#619CFF", "RNA-seq" = "#F8766D")) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   legend.position = "bottom")

  # Add asterisks for significantly changed genes
  sig_genes <- plot_data_long[plot_data_long$Significant == "Yes" & plot_data_long$Platform == "qPCR", ]
  if (nrow(sig_genes) > 0) {
    p <- p + ggplot2::geom_text(data = sig_genes,
                                ggplot2::aes(x = Gene, y = log2FC + 0.2, label = "*"),
                                size = 6, inherit.aes = FALSE)
  }

  print(p)
  return(invisible(p))
}
