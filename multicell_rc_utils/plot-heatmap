#!/usr/bin/env Rscript

library(tidyverse)
library(reshape2)
library(dplyr)
library(grid)
library(gridExtra)
library(stringr)
library(viridis)
library(viridisLite)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
 cat("Missing params to plot.")
 q()
}

results <- read.csv('results.csv')
results <- results[, -1]

results <- select(results, c(Accuracy, args[1], args[2]))
results <- aggregate(results, by = results[args], FUN = mean)
results <- results[,-length(colnames(results))]
results <- results[,-length(colnames(results))]
results$Accuracy <- round(results$Accuracy, 2)

results %>% ggplot(aes_string(x=args[1], y=args[2], fill="Accuracy", label="Accuracy")) +
  scale_fill_viridis(option = "inferno") +
  geom_tile() +
  geom_text() +
  theme_bw() +
  theme(text = element_text(size = 26), axis.title.y = element_text(angle = 0))
ggsave('plot.png', width = 20, height = 15)
