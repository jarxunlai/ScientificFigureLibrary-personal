# =============================================================================
# 同义替换率Ks密度直方图
# organized：线性脚本 + 中文分节；路径相对本条目目录
# =============================================================================
# =============================================================================


script_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg))) else normalizePath(getwd())
  }
)
root <- if (basename(script_dir) == "code") dirname(script_dir) else script_dir
# 数据路径：file.path(root, "data", ...)
out_dir <- file.path(root, "output", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# 载入R包：
library(ggplot2)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(stringr)
  library(forcats)
})


######## 构造模拟数据 ---------------
data <- data.frame(Cmo_Cse = rnorm(100, mean = 0.15, sd = 0.03),
                   Cmotri_Cmotri = rnorm(100, mean = 0.1, sd = 0.04),
                   Cse_Cse = rnorm(100, mean = 0.2, sd = 0.2),
                   Cmo_Cmo = rnorm(100, mean = 0.15, sd = 0.2),
                   Cna_Cna = rnorm(100, mean = 0.2, sd = 0.3),
                   Han_Han = rnorm(100, mean = 0.6, sd = 0.1),
                   Ccar_Ccar = rnorm(100, mean = 1, sd = 0.3),
                   Cmo_Ccar = rnorm(100, mean = 0.8, sd = 0.1))
head(data)
#     Cmo_Cse Cmotri_Cmotri     Cse_Cse   Cmo_Cmo    Cna_Cna   Han_Han Ccar_Ccar
# 1 0.1609222    0.05261231  0.16892730 0.1284194 -0.1312250 0.6101697 0.6331818
# 2 0.1404127    0.11737022  0.40848584 0.1803633 -0.2079286 0.6551174 0.9060431
# 3 0.1687347    0.09937649 -0.10444116 0.4073475  0.1882052 0.5268000 1.0351207
# 4 0.1203003    0.10224016  0.13193332 0.2518172 -0.1600887 0.8196790 1.8970021
# 5 0.1731269    0.04150252  0.08036157 0.1634842  0.4217774 0.6291792 1.4272142
# 6 0.1757105    0.10123956  0.18602537 0.6402353 -0.1095695 0.6001141 0.5667011

colors <- c("#d0b9d5", "#f4a7aa", "#fca1c2", "#fcd46e",
            "#c8e27f", "#fed2a1", "#94cf92", "#60b4e2")

######## 整体一张图 -- 普通 -----------------
p <- ggplot(data)+
  geom_histogram(aes(x = Cmo_Cse, y=after_stat(density)),
                 color = "white", fill = colors[1], bins = 60,
                 alpha = 0.9)+
  geom_density(aes(x = Cmo_Cse),
               color = colors[1], fill = NA)+
  ylab("Density")+
  xlab("Ks")+
  xlim(0,2)+
  scale_y_continuous(breaks = seq(0, 15, 2.5))+
  theme_classic()+
  theme(text = element_text(face = "bold"))

for (i in 2:ncol(data)) {
  tmp_data <- data
  colnames(tmp_data)[i] <- "X"
  if (i == 8) {
    p <- p +
      geom_histogram(data = tmp_data, aes(x = X, y=after_stat(density)),
                     color = "white", fill = colors[i], bins = 60,
                     alpha = 0.8)+
      geom_density(data = tmp_data, aes(x = X),
                   color = colors[i], fill = NA)

  } else {
    p <- p +
      geom_density(data = tmp_data, aes(x = X),
                   color = NA, fill = colors[i],
                   alpha = 0.8)
  }
}

p

ragg::agg_png(file.path(out_dir, "ks_density.png"), width = 6, height = 4, units = "in", res = 300)
print(last_plot())
dev.off()
file.copy(file.path(out_dir, "ks_density.png"), file.path(root, "preview.png"), overwrite = TRUE)
ggsave(file.path(out_dir, "ks_density.pdf"), last_plot(), width = 6, height = 4)








