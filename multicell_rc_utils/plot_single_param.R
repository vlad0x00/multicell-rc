#!/usr/bin/env Rscript

library(tidyverse)
library(reshape2)
library(dplyr)
library(grid)
library(gridExtra)

results <- read.csv('results.csv', stringsAsFactors = FALSE)
rownames(results) <- results[, 1]
results <- results[, -1]
results$Secretion.high <- factor(results$Secretion.high)
results$Cytokines <- factor(results$Cytokines)

results$Cells <- floor(results$Cells / 10)
#cells <- paste(results$Cells, 'cells')
#cells_levels <- paste(unique(sort(results$Cells)), 'cells')
#results$Cells <- factor(cells, levels = cells_levels)

results %>% ggplot(aes(x=Cells, y=Accuracy, group=1)) +
  stat_summary(geom="ribbon", fun.data=mean_cl_normal,fun.args=list(conf.int=0.95), fill="lightblue", alpha=0.5) +
  stat_summary(geom="line", fun=mean, linetype="dashed") +
  stat_summary(geom="point", fun=mean, color="red") +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 0))
ggsave('plot.png', width = 20, height = 15)
