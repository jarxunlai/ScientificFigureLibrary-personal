# =============================================================================
# 堆积与百分比柱状图
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
out_dir <- file.path(root, "output", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(stringr)
  library(forcats)
})

# 构造模拟数据：
group1 <- sample(10:50, 40, replace = T)
group2 <- sample(5:30, 40, replace = T)

for (i in 3:15) {
  group <- sample(5:15, 40, replace = T)
  assign(paste0("group", i), group)
}

# 合并：
group_vars <- factor(ls(pattern = "group[0-9][0-9]*"),
                     levels = paste0("group", 1:15))
data <- as.data.frame(Reduce(cbind, lapply(as.character(sort(group_vars)), get)))
colnames(data) <- as.character(sort(group_vars))
data$x <- paste0("BRH148", 1:40)

# 宽转长:
# library(tidyverse) 已在文件头拆包

data_long <- data %>%
  pivot_longer(cols = -x,
               names_to = "group",
               values_to = "value") %>%
  group_by(x) %>%
  mutate(sum = sum(value)) %>%
  arrange(desc(sum))
data_long$x <- factor(data_long$x, levels = unique(data_long$x))
data_long$group <- factor(data_long$group, levels = paste0("group", 15:1))


# 画图：
library(RColorBrewer)

fills <- colorRampPalette(brewer.pal(8, "Set3"))(15)
p1 <- ggplot(data_long)+
  geom_bar(aes(x, value, fill = group),
           stat = "identity", position = "stack")+
  scale_fill_manual(values = rev(fills),
                    name = "",
                    guide = guide_legend(ncol = 5,
                                         keywidth = 0.5,
                                         keyheight = 0.5,
                                         reverse = T))+
  xlab("")+
  ylab("PAR(Peak Area Ratio)")+
  theme_classic()+
  theme(legend.position = c(0.6, 0.95),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p1
ggsave(file.path(out_dir, "plot1.pdf"), height = 4, width = 8)

p2 <- ggplot(data_long)+
  geom_bar(aes(x, value, fill = group),
           stat = "identity", position = "fill")+
  scale_fill_manual(values = rev(fills))+
  xlab("Subject ID")+
  ylab("% of total PAR per Subject")+
  theme_classic()+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90,vjust = 0.5))
p2
ggsave(file.path(out_dir, "plot2.pdf"), height = 4, width = 8)

library(cowplot)
plot_grid(p1,p2, ncol = 1)

ggsave(file.path(out_dir, "plot.pdf"), height = 8, width = 8)
