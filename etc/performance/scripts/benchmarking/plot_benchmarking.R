#!/usr/bin/env Rscript
library(ggplot2)
library(tidyr)
library(RColorBrewer)

# args[1]: Working directory
# args[2]: Number of software under benchmarking

# Get command line arguments
args = commandArgs(trailingOnly=TRUE)

# Set working directory
wd = args[1]
setwd(wd)

# Create directory
dir.create('plots/')

# Read raw time data
raw_time = read.csv("time.csv", header=TRUE, row.names=1)

###
### [==> Processing time data
###

# Get number of replicates per software
replicates = ( ncol(raw_time) - 2 ) / as.numeric(args[2])

# Get number of threads
nthreads = colnames(raw_time)[3:(3+replicates-1)]
for (i in 1:replicates) {
    nthreads[i] = strsplit(nthreads[i], "_")[[1]][2]
}
nthreads = as.factor(as.numeric(nthreads))
nthreads = rep(rep(nthreads, each = 1000), as.numeric(args[2]))

# Get software
software = colnames(raw_time)[3:length(colnames(raw_time))]
for (i in 1:length(software)) {
    software[i] = strsplit(software[i], "_")[[1]][1]
}
software = unique(software)
software = as.factor(rep(software, each = 1000 * length(unique(nthreads))))

# Get average time
average_time = gather(raw_time[,3:dim(raw_time)[2]])$value

# Create time data
time = data.frame(nthreads, average_time, software)

###
### [==> Plotting: plots/boxplot_time_vs_nthreads.svg)
###

# Boxplot grouped by software: Time x Number of threads
p <- ggplot(time, aes(x=nthreads, y=average_time, fill=software)) +
        geom_boxplot() +
        theme_bw() +
        labs(x="Number of OpenMP threads", y="Time (s)", fill="Software") +    
        scale_y_continuous(breaks = c(seq(0, round(range(time$average_time)[2] / 10) * 10, by = 10)), expand = c(0,5), limits = c(0, round(range(time$average_time)[2] / 10) * 10)) +
        theme(axis.text.x = element_text(size = 14), 
            axis.text.y = element_text(size = 18, hjust = 1),
            legend.text = element_text(size = 12)) +
        theme(axis.title.y = element_text(angle = 90, 
                                        hjust = 0.5, 
                                        size = 20),
            axis.title.x = element_text(size = 20),
            legend.title = element_text(size = 12))

svg('plots/boxplot_time_vs_nthreads.svg')
p
dev.off()

png('plots/boxplot_time_vs_nthreads.png', width = 3200, height = 2100, res = 300)
p
dev.off()
