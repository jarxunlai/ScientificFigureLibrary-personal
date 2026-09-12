# 系统发育树与关联注释
# 路径由脚本位置确定；仅使用包内示例，不安装依赖。
# 已执行并生成预览；保持作者数据处理、函数、因子顺序、配色及布局，不自动安装依赖。
# 敏感/未测定/缺失以及WT/未列入突变必须结合原始数据字典再审查。

# 一、依赖和输出位置 -----------------------------------------------------------
library(ggtree)
library(ggtreeExtra)
library(ggplot2)
library(ggnewscale)
library(reshape2)
library(dplyr)
library(tidytree)
library(ggstar)
library(TDbook)
library(ragg)

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
out_dir <- file.path(root, "output", "figures-public")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


out_dir <- file.path(root, "output", "figures-public")
png_path <- file.path(out_dir, "cauris-resistance.png")
pdf_path <- file.path(out_dir, "cauris-resistance.pdf")
if (any(file.exists(c(png_path, pdf_path)))) {
  stop("输出已存在；停止以免覆盖。请先审核现有产物并明确选择新的输出名。")
}

# 二、读取作者包内的树和表，不重新建树 -----------------------------------------
# load tr and dat from the TDbook package 
dat <- df_Candidaauris_data
tr <- tree_Candidaauris

countries <- c("Canada", "United States",
               "Colombia", "Panama",
               "Venezuela", "France",
               "Germany", "Spain",
               "UK", "India",
               "Israel", "Pakistan",
               "Saudi Arabia", "United Arab Emirates",
               "Kenya", "South Africa",
               "Japan", "South Korea",
               "Australia")
# For the tip points
dat1 <- dat %>% select(c("ID", "COUNTRY", "COUNTRY__colour"))
dat1$COUNTRY <- factor(dat1$COUNTRY, levels=countries)
COUNTRYcolors <- dat1[match(countries,dat$COUNTRY),"COUNTRY__colour"]

# 三、药敏转长表；原文Not_合并规则保留，含义仍待数据字典核对 ----------------------
# For the heatmap layer
dat2 <- dat %>% select(c("ID", "FCZ", "AMB", "MCF"))
dat2 <- melt(dat2,id="ID", variable.name="Antifungal", value.name="type")
dat2$type <- paste(dat2$Antifungal, dat2$type)
dat2$type[grepl("Not_", dat2$type)] = "Susceptible"
dat2$Antifungal <- factor(dat2$Antifungal, levels=c("FCZ", "AMB", "MCF"))
dat2$type <- factor(dat2$type,
                    levels=c("FCZ Resistant",
                            "AMB Resistant",
                            "MCF Resistant",
                            "Susceptible"))

# 四、靶标变异转长表；WT不画点，未画点不能自行解释为已证实WT ----------------------
# For the points layer
dat3 <- dat %>% select(c("ID", "ERG11", "FKS1")) %>%
        melt(id="ID", variable.name="point", value.name="mutation")
dat3$mutation <- paste(dat3$point, dat3$mutation)
dat3$mutation[grepl("WT", dat3$mutation)] <- NA
dat3$mutation <- factor(dat3$mutation, 
                        levels=c("ERG11 Y132F", "ERG11 K143R",
                                 "ERG11 F126L", "FKS1 S639Y/P/F"))

# 五、作者的支系分组函数原样保留 -----------------------------------------------
# For the clade group
dat4 <- dat %>% select(c("ID", "CLADE"))
dat4 <- aggregate(.~CLADE, dat4, FUN=paste, collapse=",")
clades <- lapply(dat4$ID, function(x){unlist(strsplit(x,split=","))})
names(clades) <- dat4$CLADE

tr <- groupOTU(tr, clades, "Clade")
Clade <- NULL
p <- ggtree(tr=tr, layout="fan", open.angle=15, size=0.2, aes(colour=Clade)) +
     scale_colour_manual(
         values=c("black","#69B920","#9C2E88","#F74B00","#60C3DB"),
         labels=c("","I", "II", "III", "IV"),
         guide=guide_legend(keywidth=0.5,
                            keyheight=0.5,
                            order=1,
                            override.aes=list(linetype=c("0"=NA,
                                                         "Clade1"=1,
                                                         "Clade2"=1,
                                                         "Clade3"=1,
                                                         "Clade4"=1
                                                        )
                                             )
                           )
     ) + 
     new_scale_colour()

# 六、国家信息在叶端文字上；点alpha=0，不可见 -----------------------------------
p1 <- p %<+% dat1 +
     geom_tippoint(aes(colour=COUNTRY),
                   alpha=0) +
     geom_tiplab(aes(colour=COUNTRY),
                   align=TRUE,
                   linetype=3,
                   size=1,
                   linesize=0.2,
                   show.legend=FALSE
                   ) +
     scale_colour_manual(
         name="Country labels",
         values=COUNTRYcolors,
         guide=guide_legend(keywidth=0.5,
                            keyheight=0.5,
                            order=2,
                            override.aes=list(size=2,alpha=1))
     )

# 七、FCZ/AMB/MCF药敏热图 ------------------------------------------------------
p2 <- p1 +
      geom_fruit(
          data=dat2,
          geom=geom_tile,
          mapping=aes(x=Antifungal, y=ID, fill=type),
          width=0.1,
          color="white",
          pwidth=0.1,
          offset=0.15
      ) +
      scale_fill_manual(
           name="Antifungal susceptibility",
           values=c("#595959", "#B30000", "#020099", "#E6E6E6"),
           na.translate=FALSE,
           guide=guide_legend(keywidth=0.5,
                              keyheight=0.5,
                              order=3
                             )
      ) +
      new_scale_fill()

# 八、ERG11/FKS1靶标突变轨道 --------------------------------------------------
p3 <- p2 +
      geom_fruit(
          data=dat3,
          geom=geom_star,
          mapping=aes(x=mutation, y=ID, fill=mutation, starshape=point),
          size=1,
          starstroke=0,
          pwidth=0.1,
          inherit.aes = FALSE,
          grid.params=list(
                          linetype=3,
                          size=0.2
                      )

      ) +
      scale_fill_manual(
          name="Point mutations",
          values=c("#329901", "#0600FF", "#FF0100", "#9900CC"),
          guide=guide_legend(keywidth=0.5, keyheight=0.5, order=4,
                             override.aes=list(
                                    starshape=c("ERG11 Y132F"=15,
                                                "ERG11 K143R"=15,
                                                "ERG11 F126L"=15,
                                                "FKS1 S639Y/P/F"=1),
                                    size=2)
                            ),
          na.translate=FALSE,
      ) +
      scale_starshape_manual(
          values=c(15, 1),
          guide="none"
      ) +
      theme(
          legend.background=element_rect(fill=NA),
          legend.title=element_text(size=7), 
          legend.text=element_text(size=5.5),
          legend.spacing.y = unit(0.02, "cm")
      )
p3

# 九、显式导出；不替换原图、不自动设canonical preview ---------------------------
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
ggsave(png_path, plot=p3, width=12, height=12, units="in", dpi=300,
       bg="white", device=ragg::agg_png)
ggsave(pdf_path, plot=p3, width=12, height=12, units="in", bg="white",
       device=grDevices::cairo_pdf)
stopifnot(all(file.info(c(png_path, pdf_path))$size > 0))
packages <- c("ggtree", "ggtreeExtra", "ggplot2", "ggnewscale", "reshape2",
              "dplyr", "tidytree", "ggstar", "TDbook", "ragg")
# 原配色及小字号尚未通过可读性/色觉检查；出图成功也不等于能发布或科学结论成立。

stopifnot(file.copy(png_path,file.path(root,"preview.png"),overwrite=TRUE))
