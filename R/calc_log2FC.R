#' Calculate log2 Fold Change from qPCR Expression Data
#'
#' Computes treatment vs control group log2 fold change for specified
#' genes using raw Expression values from the bundled data file.
#' The Expression values are assumed to be normalized such that the
#' control group mean equals 1 (i.e., relative expression = 2^(-ddCt)
#' divided by mean control 2^(-ddCt)).
#'
#' @param genes Character vector of gene names. Must match column names
#'   in the expression matrix after correction (e.g., "NFkB" not "NF-kB").
#' @param treatment Character. Stress model to analyse:
#'   `"reserpine"` or `"chronic stress (CS)"`.
#'
#' @return A data.frame with columns `Gene` and `log2FC`.
#'   `log2FC` is rounded to 4 decimal places.
#'   If a gene has fewer than 2 valid observations in either group,
#'   its log2FC is returned as `NA`.
#'
#' @export
#' @importFrom readxl read_excel
#' @importFrom tidyr fill
#' @importFrom stats na.omit
#'
#' @examples
#' \dontrun{
#' # Compute log2FC for selected genes in the reserpine model
#' calc_log2FC(c("TLR4", "MyD88", "FADD"), "reserpine")
#'
#' # Compute log2FC for all 15 genes in the chronic stress model
#' gene_list <- c("TLR2","TLR3","TLR4","MyD88","IRAK1","IRAK4",
#'                "TRAF3","TRAF6","TAK1","NEMO","NFkB","TRIF",
#'                "TBK1","MAVS","FADD")
#' calc_log2FC(gene_list, "chronic stress (CS)")
#' }
calc_log2FC <- function(genes, treatment = c("reserpine", "chronic stress (CS)")) {
  treatment <- match.arg(treatment)

  # Step 1: Load raw expression data from the bundled Excel file
  data_path <- system.file("extdata", "expression_data_zhixuan_He.xlsx",
                           package = "TLR.zebrafish.stress")
  raw <- readxl::read_excel(data_path, sheet = 1)

  # Step 2: Pre-process data (fill merged cells, correct column names,
  #         remove rows that are almost entirely empty)
  raw <- tidyr::fill(raw, Treatment, .direction = "down")
  raw <- tidyr::fill(raw, Categories, .direction = "down")
  colnames(raw)[grep("NF-", colnames(raw))] <- "NFkB"
  raw <- raw[rowSums(is.na(raw[, 4:18])) < 15, ]

  # Step 3: Subset to the specified treatment model and split by group
  sub_data <- raw[raw[["Treatment"]] == treatment, ]
  ctrl  <- sub_data[sub_data[["Categories"]] == "control group", ]
  treat <- sub_data[sub_data[["Categories"]] == "experimental group", ]

  # Step 4: Compute log2FC per gene
  #         Because the control group is normalized to mean = 1,
  #         log2FC = log2(mean_treat / mean_ctrl) = log2(mean_treat)
  result_list <- lapply(genes, function(gene) {
    ctrl_vals  <- as.numeric(stats::na.omit(ctrl[[gene]]))
    treat_vals <- as.numeric(stats::na.omit(treat[[gene]]))

    # Require at least 2 observations in each group for a meaningful mean
    if (length(ctrl_vals) < 2 || length(treat_vals) < 2) {
      lfc <- NA_real_
    } else {
      lfc <- round(log2(mean(treat_vals)), 4)
    }
    data.frame(Gene = gene, log2FC = lfc, stringsAsFactors = FALSE)
  })

  do.call(rbind, result_list)
}
