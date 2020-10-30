#!/usr/bin/env Rscript

library(tidyverse)
library(reshape2)
library(dplyr)
library(grid)
library(gridExtra)
library(stringr)

results <- read.csv('results.csv', stringsAsFactors = FALSE)
rownames(results) <- results[, 1]
results <- results[, -1]
results$Secretion.high <- factor(results$Secretion.high)
results$Cytokines <- factor(results$Cytokines)

results$Genes <- floor(results$Genes / 10) * 10
results$Tissue.width <- floor(sqrt(results$Cells / 3))
results$Cells <- floor(results$Cells / 10) * 10

# The code below adds a suffix 'cells' to the cells column rows
#cells <- paste(results$Cells, 'cells')
#cells_levels <- paste(unique(sort(results$Cells)), 'cells')
#results$Cells <- factor(cells, levels = cells_levels)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  cat("Missing param to plot.")
  q()
}

results %>% ggplot(aes_string(x=args[1], y="Accuracy", group=1)) +
  #stat_summary(geom="ribbon", fun.data=mean_cl_normal,fun.args=list(conf.int=0.95), fill="lightblue", alpha=0.5) +
  stat_summary(geom="line", fun=mean, linetype="dashed") +
  stat_summary(geom="point", fun=mean, color="red") +
  geom_smooth(formula=y ~ x, method=stats::loess) +
  theme_bw() +
  #facet_wrap(~Window.size) +
  theme(text = element_text(size = 26), axis.title.y = element_text(angle = 0))
ggsave('plot.png', width = 20, height = 15)
