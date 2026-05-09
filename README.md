# TLR.zebrafish.stress

An R package for analyzing the expression of 15 Toll-like receptor (TLR) 
signaling pathway genes in the zebrafish (*Danio rerio*) brain under two 
stress models: reserpine-induced monoamine depletion and chronic 
unpredictable stress (CS).

**Author:** Zhixuan He (2255189)  
**Supervisor:** Prof. Lee Wei Lim  
**Institution:** Department of Biosciences and Bioinformatics, XJTLU

## Installation

```r
if (!require("devtools")) install.packages("devtools")
devtools::install_github("ZhixuanHe22/TLR.zebrafish.stress", build_vignettes = TRUE)
```

## Functions

| Function | Description |
|----------|-------------|
| `run_mann_whitney()` | Mann-Whitney U test for comparing treatment vs control gene expression |
| `plot_qpcr_vs_rnaseq()` | Bar plot comparing qPCR and RNA-seq log2 fold changes |
| `plot_go_comparison()` | Bubble plot comparing GO enrichment between two stress models |
| `calc_log2FC()` | Compute log2 fold change from raw Expression data |

## File Structure

```
├── R/
│   ├── mann_whitney.R
│   ├── calc_log2FC.R
│   ├── compare_expression.R
│   └── plot_go_comparison.R
├── man/
├── inst/extdata/
│   ├── expression_data_zhixuan_He.xlsx
│   └── BGI/
├── vignettes/
├── DESCRIPTION
├── NAMESPACE
├── LICENSE
└── README.md
```

## Reproducing the Results

### 1. Mann-Whitney U Test

```r
library(TLR.zebrafish.stress)
library(readxl)
library(tidyr)

data_path <- system.file("extdata", "expression_data_zhixuan_He.xlsx", package = "TLR.zebrafish.stress")
raw <- read_excel(data_path, sheet = 1)

raw <- raw %>% fill(Treatment, .direction = "down") %>% fill(Categories, .direction = "down")
colnames(raw)[grep("NF-", colnames(raw))] <- "NFkB"
raw <- raw[rowSums(is.na(raw[, 4:18])) < 15, ]
raw <- raw[raw$Treatment %in% c("reserpine", "chronic stress (CS)"), ]

gene_list <- c("TLR2","TLR3","TLR4","MyD88","IRAK1","IRAK4",
               "TRAF3","TRAF6","TAK1","NEMO","NFkB","TRIF",
               "TBK1","MAVS","FADD")

# Reserpine model
res_data <- raw[raw$Treatment == "reserpine", ]
run_mann_whitney(res_data, gene_list, "Categories", "Reserpine Model")

# Chronic stress model
cs_data <- raw[raw$Treatment == "chronic stress (CS)", ]
run_mann_whitney(cs_data, gene_list, "Categories", "Chronic Stress Model")
```

### 2. GO Enrichment Comparison Plot

```r
res_go <- system.file("extdata/BGI/Reserpine_GO/GO_enrichment/Control-vs-Reserpine/Control-vs-Reserpine_gene_go_P_rich_term.xls", package = "TLR.zebrafish.stress")
cs_go <- system.file("extdata/BGI/CS_GO/GO_enrichment/Control-vs-Stress/Control-vs-Stress_gene_go_P_rich_term.xls", package = "TLR.zebrafish.stress")

plot_go_comparison(reserpine_file = res_go, cs_file = cs_go, title = "GO BP: Reserpine vs Chronic Stress")
```

### 3. qPCR vs RNA-seq Comparison Plot (Reserpine Model)

```r
# Compute log2FC dynamically from raw Expression data
qpcr_res <- calc_log2FC(gene_list, treatment = "reserpine")

# Merge significance from Mann-Whitney test
mw_res <- run_mann_whitney(res_data, gene_list, "Categories", "Reserpine Model")
qpcr_res <- merge(qpcr_res, mw_res[, c("Gene", "Significant")], by = "Gene", all.x = TRUE)
qpcr_res$Significant[is.na(qpcr_res$Significant)] <- "No"

bgi_path <- system.file("extdata/BGI/Reserpine_DEG/Diff_exp/gene_diff.xls", package = "TLR.zebrafish.stress")
bgi_raw <- read.delim(bgi_path, header = TRUE, sep = "\t")

name_map <- c("TLR2"="tlr2","TLR3"="tlr3","TLR4"="tlr4",
              "MyD88"="myd88","IRAK1"="irak1","IRAK4"="irak4",
              "TRAF3"="traf3","TRAF6"="traf6","TAK1"="tak1",
              "NEMO"="nemo","NFkB"="nfkb1","TRIF"="ticam1",
              "TBK1"="tbk1","MAVS"="mavs","FADD"="fadd")

rnaseq_res <- data.frame(Gene = gene_list, log2FC = NA)
for (i in seq_along(gene_list)) {
  bgi_name <- name_map[gene_list[i]]
  row <- bgi_raw[bgi_raw$gene_symbol == bgi_name, ]
  if (nrow(row) > 0) rnaseq_res$log2FC[i] <- row$diffexp_log2fc_Control.vs.Reserpine[1]
}


plot_qpcr_vs_rnaseq(qpcr_res, rnaseq_res, title = "qPCR vs RNA-seq: Reserpine Model")
```

### 4. qPCR vs RNA-seq Comparison Plot (Chronic Stress Model)

```r
# Compute log2FC dynamically from raw Expression data
qpcr_cs <- calc_log2FC(gene_list, treatment = "chronic stress (CS)")

# Merge significance from Mann-Whitney test
cs_mw <- run_mann_whitney(cs_data, gene_list, "Categories", "Chronic Stress Model")
qpcr_cs <- merge(qpcr_cs, cs_mw[, c("Gene", "Significant")], by = "Gene", all.x = TRUE)
qpcr_cs$Significant[is.na(qpcr_cs$Significant)] <- "No"

name_map <- c("TLR2"="tlr2","TLR3"="tlr3","TLR4"="tlr4",
              "MyD88"="myd88","IRAK1"="irak1","IRAK4"="irak4",
              "TRAF3"="traf3","TRAF6"="traf6","TAK1"="tak1",
              "NEMO"="nemo","NFkB"="nfkb1","TRIF"="ticam1",
              "TBK1"="tbk1","MAVS"="mavs","FADD"="fadd")
              
bgi_cs_path <- system.file("extdata/BGI/CS_DEG/Diff_exp/gene_diff.xls", package = "TLR.zebrafish.stress")
bgi_cs <- read.delim(bgi_cs_path, header = TRUE, sep = "\t")

rnaseq_cs <- data.frame(Gene = gene_list, log2FC = NA)
for (i in seq_along(gene_list)) {
  bgi_name <- name_map[gene_list[i]]
  match_row <- bgi_cs[bgi_cs$gene_symbol == bgi_name, ]
  if (nrow(match_row) > 0) {
    lfc_col <- grep("log2fc", colnames(bgi_cs), ignore.case = TRUE, value = TRUE)
    if (length(lfc_col) == 0) lfc_col <- colnames(bgi_cs)[4]
    rnaseq_cs$log2FC[i] <- as.numeric(match_row[[lfc_col[1]]])[1]
  }
}


plot_qpcr_vs_rnaseq(qpcr_cs, rnaseq_cs, title = "qPCR vs RNA-seq: Chronic Stress Model")
```

For the full tutorial, see `browseVignettes("TLR.zebrafish.stress")`.

## Key Findings
- **Reserpine Model:** *TLR4* and *MyD88* significantly upregulated; *FADD* significantly downregulated (p < 0.05)
- **Chronic Stress Model:** *FADD* significantly upregulated (p < 0.05)
- These results reveal stress-specific regulation of the TLR signaling pathway in the zebrafish brain.

## AI Use Declaration
In the development of this package, AI tools were used for debugging R code,
drafting initial versions of function documentation, and explaining error
messages. All AI-generated content was reviewed and modified by the author
to ensure accuracy and relevance to the project.

## Contact
Zhixuan.He22@student.xjtlu.edu.cn

## Session Information

R version 4.5.3 (2026-03-11 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 11 x64 (build 26200)

Matrix products: default
  LAPACK version 3.12.1

locale:
[1] LC_COLLATE=Chinese (Simplified)_China.utf8 
[2] LC_CTYPE=Chinese (Simplified)_China.utf8   
[3] LC_MONETARY=Chinese (Simplified)_China.utf8
[4] LC_NUMERIC=C                               
[5] LC_TIME=Chinese (Simplified)_China.utf8    

time zone: Asia/Hong_Kong
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils    
[5] datasets  methods   base     

other attached packages:
[1] TLR.zebrafish.stress_1.0.0
[2] tidyr_1.3.2               
[3] readxl_1.4.5              
[4] fs_2.1.0                  

loaded via a namespace (and not attached):
 [1] generics_0.1.4     xml2_1.5.2        
 [3] stringi_1.8.7      digest_0.6.39     
 [5] magrittr_2.0.4     evaluate_1.0.5    
 [7] grid_4.5.3         RColorBrewer_1.1-3
 [9] pkgload_1.5.2      fastmap_1.2.0     
[11] cellranger_1.1.0   rprojroot_2.1.1   
[13] plyr_1.8.9         processx_3.8.6    
[15] pkgbuild_1.4.8     sessioninfo_1.2.3 
[17] brio_1.1.5         ps_1.9.1          
[19] purrr_1.2.1        scales_1.4.0      
[21] cli_3.6.6          rlang_1.2.0       
[23] pak_0.9.5          ellipsis_0.3.3    
[25] remotes_2.5.0      withr_3.0.2       
[27] cachem_1.1.0       yaml_2.3.12       
[29] devtools_2.5.2     otel_0.2.0        
[31] tools_4.5.3        reshape2_1.4.5    
[33] memoise_2.0.1      dplyr_1.2.1       
[35] ggplot2_4.0.3      credentials_2.0.3 
[37] vctrs_0.7.2        R6_2.6.1          
[39] lifecycle_1.0.5    stringr_1.6.0     
[41] usethis_3.2.1      callr_3.7.6       
[43] pkgconfig_2.0.3    desc_1.4.3        
[45] pillar_1.11.1      gtable_0.3.6      
[47] glue_1.8.0         Rcpp_1.1.1        
[49] gert_2.3.1         xfun_0.57         
[51] tibble_3.3.1       tidyselect_1.2.1  
[53] rstudioapi_0.18.0  knitr_1.51        
[55] farver_2.1.2       htmltools_0.5.9   
[57] rmarkdown_2.31     labeling_0.4.3    
[59] testthat_3.3.2     compiler_4.5.3    
[61] roxygen2_8.0.0     S7_0.2.1          
[63] askpass_1.2.1      openssl_2.3.5   
