#!/usr/bin/env Rscript

library(tidyverse)
library(reshape2)
library(dplyr)
library(grid)
library(gridExtra)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  cat("Missing param to plot.")
  q()
}

results <- read.csv('results.csv')

results$Tissue.depth <- factor(results$Tissue.depth)

ggplot(results, aes_string(x=args[1], y="Accuracy")) +
  stat_summary(geom="line", fun=mean, linetype="dashed") +
  stat_summary(geom="point", fun=mean, color="red") +
  geom_smooth(formula=y ~ x, method=stats::loess) +
  theme_bw() +
  theme(text = element_text(size = 26), axis.title.y = element_text(angle = 0))
ggsave('plot.png', width = 20, height = 15)
