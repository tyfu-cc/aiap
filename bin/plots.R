#!/usr/bin/env Rscript

library(ggplot2)
library(cowplot)
library(dplyr)


args <- commandArgs(trailingOnly = T)
step <- args[1]
infile1 <- args[2]
infile2 <- args[3]
refpath <- "/atac_seq/ref/mm10_encode_pe"


# Theme to apply to each plot
custom_theme <- theme(
  plot.title = element_text(size = 14, family = "Tahoma", face = "bold", hjust = 0.5),
  text = element_text(size = 12, family = "Tahoma"),
  axis.title = element_text(face = "bold"),
  axis.text.x = element_text(size = 10, face = "bold"),
  axis.text.y = element_text(size = 10, face = "bold")
)


# ---- 2.2 ----
plot2_2 <- function(result) {
  chr <- read.table(result)
  # Move "chrM" to the last row
  if ("chrM" %in% chr$V1) {
    r <- which(chr$V1 == "chrM")
    chr <- rbind(chr[chr$V1 != "chrM", ], chr[chr$V1 == "chrM", ])
    rownames(chr) <- NULL
  }
  
  chr$V1 <- factor(nrow(chr):1, labels = rev(as.character(chr$V1)))
  chr$V2 <- chr$V2 / sum(chr$V2)
  chr$V3 <- chr$V3 / sum(chr$V3)
  colnames(chr) <- c("chromosome", "Total mapped reads", "Non-redundant uniquely mapped reads")
  chr <- data.frame(chr[1], stack(chr[2:3]))
  
  p <- 
    ggplot(chr, aes(x = ind, y = values, fill = chromosome)) +
    ggtitle("Stacked barplotm of reads percentage in each chromosome") +
    geom_bar(stat = "identity", color = "black") +
    scale_y_continuous(name = "Percentage of reads") +
    scale_x_discrete(name = "") +
    guides(fill = guide_legend(nrow = 1, byrow = T, reverse = T))+
    theme_bw() + theme_classic() + coord_flip() +
    custom_theme + 
    theme(legend.text = element_text(size = 10, face = "bold"),
          legend.title = element_text(size = 10, face = "bold"),
          legend.position = "bottom")
  
  png("plot2_2_reads_distri_in_chrom.png", height = 640, width = 6600, res = 300)
  print(p)
  dev.off()
}


# ---- 3.1 Insertion distribution ----
plot3_1 <- function(result) {
  df <- read.table(result, col.names = c("length", "freq"))
  df$freq <- df$freq / sum(df$freq)
  
  p <- 
    ggplot(df, aes(x = length, y = ..scaled.., weight = freq)) + 
    geom_density(linewidth = 1, adjust = 0.2) + 
    ggtitle("Density plot of insertion size distribution (Adjust=0.2)") + 
    theme_bw() + theme_classic() + 
    scale_y_continuous(name = "Density", limits = c(0, 1)) +
    scale_x_continuous(name = "Insertion size") + 
    custom_theme
  
  png("plot3_1_insertion_size.png", height = 1800, width = 2400, res = 300)
  print(p)
  dev.off()
}


# ---- 3.3 Peak length distribution ----
plot3_3 <- function(result) {
  df <- read.table(result, col.names = c("length", "freq"))

  df <- rbind(
    df[df$length < 1500, ], 
    c(1500, sum(df[df$length >= 1500, 2]))
  )
  df$freq <- df$freq / sum(df$freq)
  
  p <- 
    ggplot(df, aes(x = length, y = ..scaled.., weight = freq)) + 
    geom_density(linewidth = 1, adjust = 0.2) +
    ggtitle("Density plot of peaks length distribution (Adjust=0.2)") +
    theme_bw() + theme_classic() + expand_limits(x = 0, y = 0) +
    scale_y_continuous(name = "Density", limits = c(0, 1)) +
    scale_x_continuous(name = "Length of peaks", breaks = seq(0, max(df$length), by = 150)) +
    custom_theme
  
  png("plot3_3_peak_length.png", height = 2000, width = 3000, res = 300)
  print(p)
  dev.off()
}


# ---- 4.2 ----
plot4_2 <- function(promoter_result, subsample_result) {
  sample_ratio <- read.table(
    promoter_result, 
    header = T, sep = "\t"
  )$enrichment_ratio
  ref_ratios <- read.table(
    paste(refpath, "merged_coding_promoter_peak_enrichment.txt", sep = "/"), 
    header = T, sep = "\t"
  )$enrichment_ratio
  promoter_ratios <- data.frame(
    enrichment_ratio = c(ref_ratios, sample_ratio),
    class = c(rep("ENCODE PE", length(ref_ratios)), "Sample"),
    type = "Enrichment ratio in coding promoter regions"
  )
  
  sample_ratio <- read.table(
    subsample_result, 
    header = T, sep = "\t"
  )$enrichment
  ref_ratios <- read.table(
    paste(refpath, "merged_sub10M_enrichment.txt", sep = "/"),
    header = T, sep = "\t"
  )$sub10M_peak_enrichment
  subsample_ratios <- data.frame(
    enrichment_ratio = c(ref_ratios, sample_ratio),
    class = c(rep("ENCODE PE", length(ref_ratios)), "Sample"),
    type = "Subsampled 10M enrichment ratio"
  )
  
  df <- rbind(promoter_ratios, subsample_ratios)
  
  p <- 
    ggplot(df, aes(x = class, y = enrichment_ratio, fill = class)) +
    stat_boxplot(geom = "errorbar", linewidth = 1, width = 0.3, aes(colour = class)) +
    geom_boxplot(outlier.shape = NA, width = 0.3, lwd = 1, fatten = 1, aes(colour = class)) + 
    scale_x_discrete(name = "") + 
    ggtitle("Boxplot of peaks enrichment ratio") + 
    scale_fill_manual(values = c(`ENCODE PE` = "grey", Sample = "red")) + 
    scale_colour_manual(values=c(`ENCODE PE` = "black", Sample = "red")) + 
    facet_grid(. ~ type) + theme(legend.position = "none") + 
    custom_theme +
    theme(strip.text.x = element_text(size = 12, face = "bold"),
          axis.text.x = element_text(size = 10, face = "bold"),
          axis.text.y = element_text(size = 10, face= "bold"))
  
  # If the enrichment ratio of the sample is less than 130
  if (df[df$class == "Sample" & df$type == "Subsampled 10M enrichment ratio", 1] < 130) {
    p <- plot + 
      scale_y_continuous(name = "Enrichment ratio", limits = c(0, 100))
  } else {
    p <- plot + 
      scale_y_continuous(name = "Enrichment ratio")
  }
  
  png("plot4_2_peaks_enrichment_ratio.png", height = 1800, width = 2800, res = 300)
  print(p)
  dev.off()
}


# ---- main ----
if (step == "s3_1") {
  plot3_1(infile1)
} else if (step == "s3_3") {
  plot3_3(infile1)
} else if (step == "s4_2") {
  plot4_2(infile1, infile2)
}
