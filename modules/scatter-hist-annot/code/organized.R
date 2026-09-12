# 固定模拟示例的随机种子，便于重复生成预览。
set.seed(20260911)
# =============================================================================
# 分组散点直方图注释
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

# 创建示例数据(具体到你们自己的数据就是两种组学差异分析的FC值)
library(MASS)

covariance <- matrix(c(1,-0.2,-0.2,1), nrow=2, byrow=TRUE)
data <- mvrnorm(n=1000, mu=c(0.5, 3), covariance)
data <- as.data.frame(data)
colnames(data) <- c("Genome size", "GC content")

# 查看数据情况：
plot(data$`Genome size`, data$`GC content`)

# 添加分组信息：
data$group <- "Host unknown"
data$group[order(data$`Genome size`)[c(sample(10:30, 10),
                                       sample(500:1000, 300))]] <- "Bacteroidetes"
data$group[order(data$`Genome size`)[c(sample(30:60, 20),
                                       sample(200:800, 300))]] <- "Firmicutes"
data$group[order(data$`Genome size`)[c(sample(40:80, 20),
                                       sample(100:200, 20))]] <- "Actinobacteria"
data$group[order(data$`Genome size`)[sample(1:1000, 10)]] <- "Proteobacteria"

head(data)
# Genome size GC content          group
# 1   0.9136138   3.351591     Firmicutes
# 2   0.9569283   4.865890     Firmicutes
# 3   1.3357191   3.311626  Bacteroidetes
# 4  -1.3002559   4.492414 Actinobacteria
# 5  -0.2020433   2.203683   Host unknown
# 6  -0.5626978   2.839346   Host unknown

# 绘图：
# 方法一：使用psych包绘制
library(psych)
with(data, scatter.hist(`Genome size`, `GC content`))

data_tmp <- data

colnames(data_tmp) <- c("Genome_size", "GC_content", "group")
scatter.hist(Genome_size ~ GC_content, data = data_tmp)

scatter.hist(Genome_size ~ GC_content, #增加分组变量
             data=data_tmp,
             xlab="Genome size", #行坐标名
             ylab="GC content", #纵坐标名
             density = F,
             ab=F, #增加拟合直线
             correl=F,  #删除右上角的相关系数
             smooth=F,  #删除平滑曲线
             grid=F,   #删除网格线
             ellipse=F   #删除椭圆
)

# 方法二：使用ggplot2绘制：
library(ggplot2)

p1 <- ggplot(data)+
  geom_point(aes(`Genome size`, `GC content`, color = group))+
  scale_color_manual(values = c("Bacteroidetes" = "#ffbe5d",
                                "Firmicutes" = "#82a0c3",
                                "Actinobacteria" = "#f28c88",
                                "Proteobacteria" = "#77bda1",
                                "Host unknown" = "#d6d6d6"))+
  theme_classic()+
  theme(legend.position = c(0.99,0.99),
        legend.justification = c(1,1))+
  guides(colour = guide_legend(""))


# 上方直方图：
p2 <- ggplot(data)+
  geom_histogram(aes(`Genome size`), 
                 binwidth = 0.2, 
                 fill = "#9db1c1",
                 color = "#ffffff",
                 size = 1)+
  theme_classic()+
  xlab("")+
  ylab("# of vOTUs")


# 右侧直方图：
p3 <- ggplot(data)+
  geom_histogram(aes(`GC content`), 
                 binwidth = 0.2, 
                 fill = "#9db1c1",
                 color = "#ffffff",
                 size = 1)+
  theme_classic()+
  xlab("")+
  ylab("# of vOTUs")+
  coord_flip()


# 拼图：
library(cowplot)

# 共用坐标范围和面板列/行，使边缘分布与散点坐标对应。
library(patchwork)
xr <- range(data$`Genome size`)
yr <- range(data$`GC content`)
p1 <- p1 + coord_cartesian(xlim = xr, ylim = yr)
p2 <- p2 + coord_cartesian(xlim = xr)
p3 <- p3 + coord_flip(xlim = yr)
p <- wrap_plots(p2, plot_spacer(), p1, p3, ncol = 2, widths = c(4, 1), heights = c(1, 4))
ggsave(file.path(out_dir, "plot.pdf"), plot = p, height = 5, width = 8)

# 显式生成发布预览，不依赖外部图片转换。
ggplot2::ggsave(file.path(root, "preview.png"), plot = p, device = ragg::agg_png, width = 8, height = 5, units = "in", dpi = 300, bg = "white")
