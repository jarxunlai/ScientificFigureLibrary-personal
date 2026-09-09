# =============================================================================
# 脂质类别Log2FC气泡图
# organized：线性脚本 + 中文分节；路径相对本条目目录
# =============================================================================
# 来源：KS科研分享「生信绘图」高分SCI图表复现 053分组气泡图
# 本地复现：drafts/ks-shengxin-huitu-repro/053-grouped-bubble
# Pixi：项目根 default 环境；library(tidyverse) 已拆成 ggplot2/dplyr/tidyr 等。
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

########### 数据构建 --------------

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(stringr)
  library(forcats)
})


data <- as.data.frame(matrix(NA, nrow = 15, ncol = 13))

for (i in 1:13) {
  if (i < 10) {
    data[,i] <- c(runif(7, min = -2, max = -0.5),
                  runif(8, min = 0.5, max = 2.1))
  } else {
    data[,i] <- c(runif(7, min = -2, max = -0.5),
                  runif(8, min = 0.5, max = 4))
  }
}

colnames(data) <- c("TG", "Sphingolipid", "PUFA", "PS", "PHA",
                    "PG", "PE", "PC", "HBMP", "DAG", "BMP",
                    "Acyl-AA", "ACAR")

# 宽转长：
data_long <- pivot_longer(data,
                          cols = everything(),
                          names_to = "Class", values_to = "Log2FC")
data_long$Pvalue <- runif(nrow(data_long), 1e-08, 1e-02)

write.csv(data_long, file.path(root, "data", "data.csv"), row.names = FALSE)

############## 绘图 -----------------
library(ggplot2)

# 绘图：
ggplot(data_long)+
  # 散点：
  geom_point(aes(Log2FC, Class,
                 size = -log10(Pvalue), fill = Class),
             shape = 21, alpha = 0.6)+
  # 垂直虚线：
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey")+
  xlim(-2, 4) +
  # 主题调整：
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = c(0.98, 0.98),
        legend.justification = c(1, 1),
        text = element_text(face = "bold"))+
  guides(fill = "none")

ragg::agg_png(file.path(out_dir, "lipid_bubble.png"), width = 7, height = 3.5, units = "in", res = 300)
print(last_plot())
dev.off()
file.copy(file.path(out_dir, "lipid_bubble.png"), file.path(root, "preview.png"), overwrite = TRUE)
ggsave(file.path(out_dir, "lipid_bubble.pdf"), last_plot(), width = 7, height = 3.5)







