#!/usr/bin/env Rscript

library(tidyverse)
library(reshape2)
library(dplyr)
library(grid)
library(gridExtra)
library(stringr)

results_2layers <- read.csv('results_2layers.csv')
rownames(results_2layers) <- results_2layers[, 1]
results_2layers <- results_2layers[, -1]

results_3layers <- read.csv('results_3layers.csv')
rownames(results_3layers) <- results_3layers[, 1]
results_3layers <- results_3layers[, -1]

results_4layers <- read.csv('results_4layers.csv')
rownames(results_4layers) <- results_4layers[, 1]
results_4layers <- results_4layers[, -1]

results_5layers <- read.csv('results_5layers.csv')
rownames(results_5layers) <- results_5layers[, 1]
results_5layers <- results_5layers[, -1]

results <- rbind(results_2layers, results_3layers, results_4layers, results_5layers)

results <- select(results, c(Accuracy, Cell.types, Cytokines, Tissue.depth))
results <- aggregate(results, by = list(results$Cell.types, results$Cytokines, results$Tissue.depth), FUN = mean)

results %>% ggplot(aes(x=Cell.types, y=Cytokines, fill=Accuracy)) +
  geom_tile() +
  theme_bw() +
  facet_wrap(~Tissue.depth) +
  theme(text = element_text(size = 26), axis.title.y = element_text(angle = 0))
ggsave('plot.png', width = 20, height = 15)
