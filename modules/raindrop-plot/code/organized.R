# =============================================================================
# 雨滴图
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
set.seed(20260912)
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

library(ggplot2)

# library(tidyverse) 已在文件头拆包

# 模拟数据构建：两列：降雨量+州名
rain_sum_state <- data.frame(State = paste0("State", 1:30),
                             meanRainfall = runif(30, 1, 10)*10^5)

# Calculate the quantiles of the mean rainfall
rain_quantiles <- quantile(x = rain_sum_state$meanRainfall,
                           # probabilities of 0.1, 0.3, 0.5, 0.7, 0.9
                           probs = seq(0.1, 0.9, by = 0.2))

# Linetypes to randomly select from
linetypes <- c("dashed", "dotted", "dotdash", "longdash", "twodash",
               "1F", "4C88C488", "12345678")


# Create a lookup table to add in State Abbreviations including DC
state_lookup <- data.frame(State = paste0("State", 1:30),
                           State.abbr = paste0("S", 1:30))


plotting_rain <- rain_sum_state |>
  mutate(image_size = case_when(meanRainfall <= rain_quantiles[1] ~ 0.005,
                                meanRainfall <= rain_quantiles[2] ~ 0.008,
                                meanRainfall <= rain_quantiles[3] ~ 0.011,
                                meanRainfall <= rain_quantiles[4] ~ 0.014,
                                meanRainfall <= rain_quantiles[5] ~ 0.017,
                                meanRainfall > rain_quantiles[5] ~ 0.02,
                                TRUE ~ NA_real_),
         # randomly sampling linetypes
         linetype = sample(linetypes, n(), replace = TRUE)) |>
  # Left-joining state abbreviations
  left_join(state_lookup, by = "State")










# 用数学轮廓直接绘制矢量雨滴，无外部图片或网络依赖。
plotting_rain$x <- match(plotting_rain$State, sort(plotting_rain$State))
theta <- seq(0,2*pi,length.out=80)
drop_vertices <- data.frame()
for(i in seq_len(nrow(plotting_rain))){
 drop_vertices <- rbind(drop_vertices,data.frame(
 x=plotting_rain$x[i]+sin(theta)*(1-cos(theta))*plotting_rain$image_size[i]*12,
 y=plotting_rain$meanRainfall[i]-cos(theta)*plotting_rain$image_size[i]*2000000, id=i))
}

# Create the baseplot
rain_plot <- ggplot(data = plotting_rain,
                    aes(x = x, y = meanRainfall)) +
  # Line graph that shows mean rainfall from top
  geom_linerange(aes(ymin = 0, ymax = meanRainfall, linetype = linetype),
                 color = "#7FCDEE") +
  # Scatter plot with raindrop images as points
  geom_polygon(data=drop_vertices, aes(x=x,y=y,group=id), inherit.aes=FALSE, fill="#38B5DF", colour=NA) +

  # Adding text to indicate which state is which. Added a bit of a buffer so that
  # labels will appear equidistant from images, even though images are all sized
  # differently.
  geom_text(aes(label = State.abbr, y = meanRainfall + (image_size*2000000)),
            size = 4,
            color = "#38B5DF",
            vjust = 2) +
  # Reversing y-scale so that the rain "falls" down from top
  scale_y_reverse() +
  ylab("Mean Annual Intercepted Rainfall (cubic meters)") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.y = element_text(size = 10, color = "#333333"),
        axis.title.y = element_text(size = 16, color = "#333333", margin = margin(r = 20)))

ggsave(file.path(out_dir, "ggimage.png"), plot = rain_plot, width = 15, height = 6, dpi = 300, bg = "white")

ggplot2::ggsave(file.path(root,"preview.png"), plot=rain_plot, width=15, height=6, dpi=200, bg="white")
