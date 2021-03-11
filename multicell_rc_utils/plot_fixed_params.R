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
results$ESMs <- factor(results$ESMs)

cells <- paste(results$Cells, 'cells')
cells_levels <- paste(unique(sort(results$Cells)), 'cells')
results$Cells <- factor(cells, levels = cells_levels)
window_levels <- nlevels(factor(results$Window.size))

results %>% ggplot(aes(x=Accuracy)) +
  geom_histogram(binwidth=0.025) +
  theme_bw() +
  theme(text = element_text(size = 26), axis.title.y = element_text(angle = 0)) +
  facet_wrap(~ Function + Window.size, ncol=window_levels)
ggsave('plot.png', width = 20, height = 15)
