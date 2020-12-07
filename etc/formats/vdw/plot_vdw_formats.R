#!/usr/bin/env Rscript
library(ggplot2)
library(tidyr)
library(RColorBrewer)
# library(ggpmisc)

# Read data
data = read.csv("read_vdw_formats.csv", row.names=1)

stats = as.data.frame(t(data[6:7,]))
stats$formats = rownames(stats)

p <- ggplot(stats, aes(x=formats, y=AVERAGE*1000, fill=formats)) + 
    geom_bar(stat="identity", color="black", alpha=0.75, size=0.75) +
    geom_errorbar(aes(ymin=(AVERAGE*1000)-(STD*1000), ymax=(AVERAGE*1000)+(STD*1000)), width=0.2) +
    geom_text(aes(y=(AVERAGE*1000)+(STD*1000)+0.5, label=AVERAGE*1000), size=6) +
    scale_fill_brewer(palette="Dark2") +
    theme_classic() +
    labs(x=NULL, y="Time(ms)") +
    scale_y_continuous(breaks = c(seq(0, 15, by=1)), limits=c(0,15.1), expand=c(0,0)) +
    theme(legend.position="none",
        axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 18, hjust = 1),
        legend.text = element_text(size = 12)) +
    theme(axis.title.y = element_text(angle = 90, hjust = 0.5, size = 20),
        axis.title.x = element_text(size = 20),
        legend.title = element_text(size = 12))

png('read_vdw_formats.png', width = 4400, height = 2475, res = 300)
p
dev.off()
