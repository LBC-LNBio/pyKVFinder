#!/usr/bin/env Rscript
library(ggplot2)
library(tidyr)
library(ggpmisc)
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

dir.create('plots/comparison/')
dir.create('plots/comparison/time')

# Get number of replicates per software
replicates = ( ncol(raw_time) - 2 ) / as.numeric(args[2])

# Get number of atoms
natoms = rep(raw_time$natoms, ncol(raw_time) - 2)

# Get number of threads
nthreads = colnames(raw_time)[3:(3+replicates-1)]
for (i in 1:replicates) {
    nthreads[i] = strsplit(nthreads[i], "_")[[1]][2]
}
nthreads = as.numeric(nthreads)
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
time = data.frame(natoms, nthreads, average_time, software)

###
### [==> Processing average and std of time data
###

# Get number of replicates per software
replicates = ( ncol(raw_time) - 2 ) / as.numeric(args[2])

# Get number of threads
nthreads = colnames(raw_time)[3:(3+replicates-1)]
for (i in 1:replicates) {
    nthreads[i] = strsplit(nthreads[i], "_")[[1]][2]
}
nthreads = as.numeric(nthreads)
nthreads = rep(nthreads, as.numeric(args[2]))

# Get software
software = colnames(raw_time)[3:length(colnames(raw_time))]
for (i in 1:length(software)) {
    software[i] = strsplit(software[i], "_")[[1]][1]
}
software = unique(software)
software = as.factor(rep(software, each = replicates))

# Get average, sd and median time
t = colMeans(raw_time[,3:dim(raw_time)[2]])
sd = apply(raw_time[,3:dim(raw_time)[2]], 2, sd)
median = apply(raw_time[,3:dim(raw_time)[2]], 2, median)

# Create time data
time2 = data.frame(nthreads, software, t, sd, median)

###
### [==> Plotting: plots/comparison/time/boxplot_time_vs_nthreads.svg)
###   Description: Boxplot with time for different softwares
###

pd=position_dodge(.75)

# Boxplot grouped by software: Time x Number of threads
p <- ggplot(time, aes(x=as.factor(nthreads), y=average_time, fill=software)) +
        geom_boxplot() +
        geom_point(data=time2, aes(x=as.factor(nthreads), y=t, group=software), shape=21, position=pd, fill='white', size=1) +
        
        theme_bw() +
        labs(x="Number of OpenMP threads", y="Time (s)", fill="Software") +    
        scale_y_continuous(breaks = c(seq(0, ceiling(range(time$average_time)[2] / 5) * 5, by = 2.5)), expand = c(0,0), limits = c(0, ceiling(range(time$average_time)[2] / 5) * 5 + 1)) +
        theme(axis.text.x = element_text(size = 20), 
            axis.text.y = element_text(size = 20, hjust = 1),
            legend.text = element_text(size = 20)) +
        theme(axis.title.y = element_text(angle = 90, 
                                        hjust = 0.5, 
                                        size = 20),
            axis.title.x = element_text(size = 20),
            legend.title = element_text(size = 20))

svg('plots/comparison/time/boxplot_time_vs_nthreads.svg', width = 16, height = 9)
p
dev.off()

png('plots/comparison/time/boxplot_time_vs_nthreads.png', width = 4400, height = 2475, res = 300)
p
dev.off()

###
### [==> Plotting: plots/comparison/time/lineplot_time_vs_atoms.svg)
###   Description: Lineplots with time against number of atoms for different number of threads
###

p <- ggplot(time, aes(x=natoms, y=average_time, group=software)) +

    geom_jitter(aes(group=as.factor(nthreads), color=as.factor(nthreads)), alpha=0.4, size=0.5) +

    geom_smooth(aes(group=as.factor(nthreads), color=as.factor(nthreads)), method='glm', se=FALSE, fullrange=TRUE) +

    geom_smooth(color='black', method='glm', se=FALSE, fullrange=TRUE) +
    
    facet_wrap(~software, ncol=as.numeric(args[2])) +

    theme_bw() +
    labs(x = "Number of atoms", y = "Time (s)", color = "OpenMP threads") +
    scale_y_continuous(breaks = c(seq(0, ceiling(range(time$average_time)[2] / 1) * 1, by = 2.5)), expand = c(0,0), limits = c(0, ceiling(range(time$average_time)[2] / 1) * 1)) +
    scale_x_continuous(breaks = c(seq(0, ceiling(range(time$natoms)[2] / 100) * 100, by = 1000)), expand = c(0,0), limits = c(0, ceiling(range(time$natoms)[2] / 100) * 100 + 100)) +
    scale_color_manual(values = c(brewer.pal(n=replicates, name = "Paired"), "black"), limits = c(unique(as.character(time$nthreads)), "Average")) +

    theme(axis.text.x = element_text(size = 16, angle=45, hjust = 1), 
        axis.text.y = element_text(size = 20, hjust = 1),
        legend.text = element_text(size = 20),
        strip.text = element_text(size=16)) +
    theme(axis.title.y = element_text(angle = 90, hjust = 0.5, size = 20),
        axis.title.x = element_text(size = 20),
        legend.title = element_text(size = 20))

svg('plots/comparison/time/lineplot_time_vs_atoms.svg', width = 16, height = 9)
p
dev.off()

png('plots/comparison/time/lineplot_time_vs_atoms.png', width = 4400, height = 2475, res = 300)
p
dev.off()

###
### [==> Plotting: plots/comparison/time/scatter_time_vs_nthreads.svg)
###   Description: Scatter plot with time against number of atoms for different number of threads
###

mid <- mean(range(time$natoms))

p <- ggplot(time, aes(x=as.factor(nthreads), y=average_time)) +
    
    geom_jitter(aes(group=as.factor(nthreads), fill=natoms), size=2, shape=21) +
    
    facet_wrap(~software, ncol=as.numeric(args[2])) +
    
    scale_fill_gradient2(midpoint=mid, low='black', mid='azure', high='red') +

    theme_bw() +

    labs(x = "OpenMP Threads", y = "Time (s)", fill = "Number of atoms") +
    scale_y_continuous(breaks = c(seq(0, ceiling(range(time$average_time)[2] / 1) * 1, by = 5.0)), expand = c(0,0), limits = c(0, ceiling(range(time$average_time)[2] / 1) * 1)) +
    scale_x_discrete(limits = c(unique(as.character(time$nthreads)))) +
    theme(axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20, hjust = 1),
        legend.text = element_text(size = 20),
        strip.text = element_text(size=20)) +
    theme(axis.title.y = element_text(angle = 90, hjust = 0.5, size = 20),
        axis.title.x = element_text(size = 20),
        legend.title = element_text(size = 20))

svg('plots/comparison/time/scatter_time_vs_nthreads.svg', width = 16, height = 9)
p
dev.off()

png('plots/comparison/time/scatter_time_vs_nthreads.png', width = 4400, height = 2475, res = 300)
p
dev.off()

###
### [==> Plotting: plots/comparison/time/average_time_vs_nthreads.svg)
###   Description: Lineplot plot with average time against number of threads
###

p <- ggplot(time2, aes(x=as.numeric(nthreads), y=t, color=software)) +
    geom_line(aes(group=software)) +
    geom_point(shape=21, fill='white') +
    theme_bw() + 
    labs(x = "OpenMP Threads", y = "Time (s)") +
    scale_y_continuous(breaks = c(seq(0, 2.5, by = .5)), expand = c(0,0), limits = c(0, 2.6)) + 
    scale_x_continuous(breaks = c(seq(1, max(unique(nthreads)), by=1)), limits=c(0.5, max(unique(nthreads)) + .5), expand=c(0,0)) +
    theme(axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20, hjust = 1),
        legend.text = element_text(size = 20),
        strip.text = element_text(size=20)) +
    theme(axis.title.y = element_text(angle = 90, hjust = 0.5, size = 20),
        axis.title.x = element_text(size = 20),
        legend.title = element_text(size = 20))

svg('plots/comparison/time/average_time_vs_nthreads.svg', width = 16, height = 9)
p
dev.off()

png('plots/comparison/time/average_time_vs_nthreads.png', width = 4400, height = 2475, res = 300)
p
dev.off()


###
### [==> Plotting: plots/comparison/time/median_time_vs_nthreads.svg)
###   Description: Lineplot plot with median time against number of threads
###

p <- ggplot(time2, aes(x=as.numeric(nthreads), y=median, color=software)) +
    geom_line(aes(group=software)) +
    geom_point(shape=21, fill='white') +
    theme_bw() + 
    labs(x = "OpenMP Threads", y = "Time (s)") +
    scale_y_continuous(breaks = c(seq(0, 2.5, by = .5)), expand = c(0,0), limits = c(0, 2.6)) + 
    scale_x_continuous(breaks = c(seq(1, max(unique(nthreads)), by=1)), limits=c(0.5, max(unique(nthreads)) + .5), expand=c(0,0)) +
    theme(axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20, hjust = 1),
        legend.text = element_text(size = 20),
        strip.text = element_text(size=16)) +
    theme(axis.title.y = element_text(angle = 90, hjust = 0.5, size = 20),
        axis.title.x = element_text(size = 20),
        legend.title = element_text(size = 20))

svg('plots/comparison/time/median_time_vs_nthreads.svg', width = 16, height = 9)
p
dev.off()

png('plots/comparison/time/median_time_vs_nthreads.png', width = 4400, height = 2475, res = 300)
p
dev.off()

###
### [==> Processing speedup data
###

dir.create('plots/comparison/speedup')

# Get number of replicates per software
replicates = ( ncol(raw_time) - 2 ) / as.numeric(args[2])

# Get number of atoms
natoms = rep(raw_time$natoms, replicates)

# Get number of threads
nthreads = colnames(raw_time)[3:(replicates + 2)]
for (i in 1:(replicates)) {
    nthreads[i] = strsplit(nthreads[i], "_")[[1]][2]
}
nthreads = as.numeric(nthreads)
nthreads = rep(rep(nthreads, each = 1000))

# Get speedup 
speedup = raw_time[3:(2 * replicates + 2)]
for (i in 1:replicates) {
    speedup[i+replicates] = speedup[i]/speedup[i+replicates]
}
speedup[1:replicates] = NULL
for (i in 1:replicates) {
    colnames(speedup)[i] = strsplit(colnames(speedup)[i], "_")[[1]][2]
}
speedup = gather(speedup)$value

# Create time data
speedup = data.frame(natoms, nthreads, speedup)

###
### [==> Plotting: plots/comparison/speedup/boxplot_speedup_vs_nthreads.svg)
###   Description: Boxplot with speedup against number of threads
###

# Boxplot grouped by software: Time x Number of threads
p <- ggplot(speedup, aes(x=as.factor(nthreads), y=speedup, fill=as.factor(nthreads))) +
        geom_boxplot() +
        theme_bw() +
        labs(x="Number of OpenMP threads", y="Speedup", fill=NULL) +    
        scale_y_continuous(breaks = c(seq(0, 8, by = 0.5)), expand = c(0,0), limits = c(0, 8.1)) +
        theme(legend.position="none",
            axis.text.x = element_text(size = 20), 
            axis.text.y = element_text(size = 20, hjust = 1),
            legend.text = element_text(size = 20)) +
        theme(axis.title.y = element_text(angle = 90, 
                                        hjust = 0.5, 
                                        size = 20),
            axis.title.x = element_text(size = 20),
            legend.title = element_text(size = 20))

svg('plots/comparison/speedup/boxplot_speedup_vs_nthreads.svg', width = 16, height = 9)
p
dev.off()

png('plots/comparison/speedup/boxplot_speedup_vs_nthreads.png', width = 4400, height = 2475, res = 300)
p
dev.off()

###
### [==> Plotting: plots/comparison/speedup/lineplot_speedup_vs_atoms.svg)
###   Description: Lineplots with speedup against number of atoms for different number of threads
###

p <- ggplot(speedup, aes(x=natoms, y=speedup)) +

    geom_jitter(aes(group=as.factor(nthreads), color=as.factor(nthreads)), alpha=0.4, size=0.5) +

    geom_smooth(aes(group=as.factor(nthreads), color=as.factor(nthreads)), method='lm', se=FALSE, fullrange=TRUE) +

    geom_smooth(color='black', method='lm', se=FALSE, fullrange=TRUE) +

    stat_poly_eq(formula = y ~ x, aes(color = as.factor(nthreads), label = paste(..eq.label.., ..rr.label.., sep = "~~~")), parse = TRUE, size=5, vstep=0.05) +

    theme_bw() +
    labs(x = "Number of atoms", y = "Speedup", color = "OpenMP threads") +
    scale_y_continuous(breaks = c(seq(0, ceiling(range(speedup$speedup)[2] / 1) * 2, by = 0.5)), expand = c(0,0), limits = c(0, ceiling(range(speedup$speedup)[2] / 1) * 1.25)) +
    scale_x_continuous(breaks = c(seq(0, ceiling(range(speedup$natoms)[2] / 100) * 100, by = 1000)), expand = c(0,0), limits = c(0, ceiling(range(speedup$natoms)[2] / 100) * 100 + 100)) +
    scale_color_manual(values = c(brewer.pal(n=replicates, name = "Paired"), "black"), limits = c(unique(as.character(speedup$nthreads)), "Average")) +

    theme(axis.text.x = element_text(size = 20, angle=45, hjust = 1), 
        axis.text.y = element_text(size = 20, hjust = 1),
        legend.text = element_text(size = 20),
        strip.text = element_text(size=20)) +
    theme(axis.title.y = element_text(angle = 90, hjust = 0.5, size = 20),
        axis.title.x = element_text(size = 20),
        legend.title = element_text(size = 20))

svg('plots/comparison/speedup/lineplot_speedup_vs_atoms.svg', width = 16, height = 9)
p
dev.off()

png('plots/comparison/speedup/lineplot_speedup_vs_atoms.png', width = 4400, height = 2475, res = 300)
p
dev.off()

###
### [==> Processing efficiency data
###

dir.create('plots/comparison/efficiency')

# Get number of replicates per software
replicates = ( ncol(raw_time) - 2 ) / as.numeric(args[2])

# Get number of atoms
natoms = rep(raw_time$natoms, replicates)

# Get number of threads
nthreads = colnames(raw_time)[3:(replicates + 2)]
for (i in 1:(replicates)) {
    nthreads[i] = strsplit(nthreads[i], "_")[[1]][2]
}
nthreads = as.numeric(nthreads)
nthreads = rep(rep(nthreads, each = 1000))

# Get efficiency 
efficiency = raw_time[3:(2 * replicates + 2)]
for (i in 1:replicates) {
    efficiency[i+replicates] = efficiency[i]/efficiency[i+replicates]
}
efficiency[1:replicates] = NULL
for (i in 1:replicates) {
    colnames(efficiency)[i] = strsplit(colnames(efficiency)[i], "_")[[1]][2]
}
efficiency = gather(efficiency)$value

# Create time data
efficiency = data.frame(natoms, nthreads, efficiency)
efficiency$efficiency = efficiency$efficiency/efficiency$nthreads

###
### [==> Plotting: plots/comparison/efficiency/boxplot_time_vs_nthreads.svg)
###   Description: Boxplot with efficiency against number of threads
###

# Boxplot grouped by software: Time x Number of threads
p <- ggplot(efficiency, aes(x=as.factor(nthreads), y=efficiency*100, fill=as.factor(nthreads))) +
        geom_boxplot() +
        theme_bw() +
        labs(x="Number of OpenMP threads", y="Efficiency (%)", fill=NULL) +    
        scale_y_continuous(breaks = c(seq(0, 200, by = 10)), expand = c(0,0), limits = c(0, 200)) +
        theme(legend.position="none",
            axis.text.x = element_text(size = 20), 
            axis.text.y = element_text(size = 20, hjust = 1),
            legend.text = element_text(size = 20)) +
        theme(axis.title.y = element_text(angle = 90, 
                                        hjust = 0.5, 
                                        size = 20),
            axis.title.x = element_text(size = 20),
            legend.title = element_text(size = 20))

svg('plots/comparison/efficiency/boxplot_efficiency_vs_nthreads.svg', width = 16, height = 9)
p
dev.off()

png('plots/comparison/efficiency/boxplot_efficiency_vs_nthreads.png', width = 4400, height = 2475, res = 300)
p
dev.off()

###
### [==> Plotting: plots/comparison/efficiency/lineplot_efficiency_vs_atoms.svg)
###   Description: Lineplots with efficiency against number of atoms for different number of threads
###

p <- ggplot(efficiency, aes(x=natoms, y=efficiency*100)) +

    geom_jitter(aes(group=as.factor(nthreads), color=as.factor(nthreads)), alpha=0.1, size=0.5) +

    geom_smooth(aes(group=as.factor(nthreads), color=as.factor(nthreads)), method='lm', se=FALSE, fullrange=TRUE) +

    geom_smooth(color='black', method='lm', se=FALSE, fullrange=TRUE) +

    stat_poly_eq(formula = y ~ x, aes(color = as.factor(nthreads), label = paste(..eq.label.., ..rr.label.., sep = "~~~")), parse = TRUE, size=5, vstep=0.05) +

    theme_bw() +
    labs(x = "Number of atoms", y = "Efficiency (%)", color = "OpenMP threads") +
    scale_y_continuous(breaks = c(seq(0, 200, by=10)), expand = c(0,0), limits = c(0, 200)) +
    scale_x_continuous(breaks = c(seq(0, ceiling(range(efficiency$natoms)[2] / 100) * 100, by = 1000)), expand = c(0,0), limits = c(0, ceiling(range(efficiency$natoms)[2] / 100) * 100 + 100)) +
    scale_color_manual(values = c(brewer.pal(n=replicates, name = "Paired"), "black"), limits = c(unique(as.character(efficiency$nthreads)), "Average")) +

    theme(axis.text.x = element_text(size = 20, angle=45, hjust = 1), 
        axis.text.y = element_text(size = 20, hjust = 1),
        legend.text = element_text(size = 20),
        strip.text = element_text(size=20)) +
    theme(axis.title.y = element_text(angle = 90, hjust = 0.5, size = 20),
        axis.title.x = element_text(size = 20),
        legend.title = element_text(size = 20))

svg('plots/comparison/efficiency/lineplot_efficiency_vs_atoms.svg', width = 16, height = 9)
p
dev.off()

png('plots/comparison/efficiency/lineplot_efficiency_vs_atoms.png', width = 4400, height = 2475, res = 300)
p
dev.off()

###
### [==> pyKVFinder analysis: time, speedup, efficiency
###

dir.create('plots/pyKVFinder')
dir.create('plots/pyKVFinder/time')

pyKVFinder = time[time$software == 'pyKVFinder',]
pyKVFinder$speedup = rep(pyKVFinder$average_time[pyKVFinder$nthreads == 1], replicates) / pyKVFinder$average_time
pyKVFinder$efficiency = pyKVFinder$speedup / pyKVFinder$nthreads

###
### [==> Plotting: plots/pyKVFinder/time/lineplot_time_vs_atoms.svg)
###   Description: Lineplots with pyKVFinder time against number of atoms for different number of threads
###

p <- ggplot(pyKVFinder, aes(x=natoms, y=average_time)) +

    geom_jitter(aes(group=as.factor(nthreads), color=as.factor(nthreads)), alpha=0.4, size=0.5) +

    geom_smooth(aes(group=as.factor(nthreads), color=as.factor(nthreads)), method='glm', se=FALSE, fullrange=TRUE) +

    geom_smooth(color='black', method='glm', se=FALSE, fullrange=TRUE) +

    stat_poly_eq(formula = y ~ x, aes(color = as.factor(nthreads), label = paste(..eq.label.., ..rr.label.., sep = "~~~")), parse = TRUE, size=6, vstep=0.05) +
    
    theme_bw() +
    labs(x = "Number of atoms", y = "Time (s)", color = "OpenMP threads") +
    scale_y_continuous(breaks = c(seq(0, ceiling(range(pyKVFinder$average_time)[2] / 1) * 1, by = 1.0)), expand = c(0,0), limits = c(0, ceiling(range(pyKVFinder$average_time)[2] / 1) * 1.0 + 0.1)) +
    scale_x_continuous(breaks = c(seq(0, ceiling(range(pyKVFinder$natoms)[2] / 100) * 100, by = 1000)), expand = c(0,0), limits = c(0, ceiling(range(pyKVFinder$natoms)[2] / 100) * 100 + 100)) +
    scale_color_manual(values = c(brewer.pal(n=replicates, name = "Paired"), "black"), limits = c(unique(as.character(pyKVFinder$nthreads)), "Average")) +

    theme(axis.text.x = element_text(size = 20, angle=45, hjust = 1), 
        axis.text.y = element_text(size = 20, hjust = 1),
        legend.text = element_text(size = 20),
        strip.text = element_text(size=20)) +
    theme(axis.title.y = element_text(angle = 90, hjust = 0.5, size = 20),
        axis.title.x = element_text(size = 20),
        legend.title = element_text(size = 20))

svg('plots/pyKVFinder/time/lineplot_time_vs_atoms.svg', width = 16, height = 9)
p
dev.off()

png('plots/pyKVFinder/time/lineplot_time_vs_atoms.png', width = 4400, height = 2475, res = 300)
p
dev.off()

###
### [==> Plotting: plots/pyKVFinder/speedup/boxplot_speedup_vs_nthreads.svg)
###   Description: Boxplot with speedup against number of threads
###

dir.create('plots/pyKVFinder/speedup')

# Boxplot grouped by software: Time x Number of threads
p <- ggplot(pyKVFinder, aes(x=as.factor(nthreads), y=speedup, fill=as.factor(nthreads))) +
        geom_boxplot() +
        theme_bw() +
        labs(x="Number of OpenMP threads", y="Speedup", fill=NULL) +    
        scale_y_continuous(breaks = c(seq(0, ceiling(range(pyKVFinder$speedup)[2]), by = 0.5)), expand = c(0,0), limits = c(0, ceiling(range(pyKVFinder$speedup)[2]) + 0.1)) +
        theme(legend.position="none",
            axis.text.x = element_text(size = 20), 
            axis.text.y = element_text(size = 20, hjust = 1),
            legend.text = element_text(size = 20)) +
        theme(axis.title.y = element_text(angle = 90, 
                                        hjust = 0.5, 
                                        size = 20),
            axis.title.x = element_text(size = 20),
            legend.title = element_text(size = 20))

svg('plots/pyKVFinder/speedup/boxplot_speedup_vs_nthreads.svg', width = 16, height = 9)
p
dev.off()

png('plots/pyKVFinder/speedup/boxplot_speedup_vs_nthreads.png', width = 4400, height = 2475, res = 300)
p
dev.off()

###
### [==> Plotting: plots/pyKVFinder/speedup/lineplot_speedup_vs_nthreads.svg)
###   Description: Lineplot with speedup against number of threads
###

p <- ggplot(pyKVFinder, aes(x = nthreads, y = speedup)) +
    stat_summary(fun = mean, geom = "line", color = 4, linetype = 1, size = 1) +
    stat_summary(fun = mean, fun.min = function(x) mean(x) - sd(x), fun.max = function(x) mean(x) + sd(x), geom = "errorbar", width=0.2) +
    stat_summary(fun = mean, geom = "point", fill = "white", shape = 21, size = 2) +
    theme_bw() + 
    labs(x="OpenMP threads", y="Speedup", fill=NULL) +    
    scale_y_continuous(breaks = c(seq(0, ceiling(range(pyKVFinder$speedup)[2]), by = 0.5)), expand = c(0,0), limits = c(0, ceiling(range(pyKVFinder$speedup)[2]) + 0.1)) +
    scale_x_continuous(breaks = c(seq(1, max(unique(nthreads)), by=1)), limits=c(0.8, max(unique(nthreads)) + .5), expand=c(0,0)) +
    theme(legend.position="none",
        axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20, hjust = 1),
        legend.text = element_text(size = 20)) +
    theme(axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        legend.title = element_text(size = 20))

svg('plots/pyKVFinder/speedup/lineplot_speedup_vs_nthreads.svg', width = 16, height = 9)
p
dev.off()

png('plots/pyKVFinder/speedup/lineplot_speedup_vs_nthreads.png', width = 4400, height = 2475, res = 300)
p
dev.off()

###
### [==> Plotting: plots/pyKVFinder/speedup/boxplot_efficiency_vs_nthreads.svg)
###   Description: Boxplot with efficiency against number of threads
###

dir.create('plots/pyKVFinder/efficiency')

# Boxplot grouped by software: Time x Number of threads
p <- ggplot(pyKVFinder, aes(x=as.factor(nthreads), y=efficiency*100, fill=as.factor(nthreads))) +
        geom_boxplot() +
        theme_bw() +
        labs(x="Number of OpenMP threads", y="Efficiency (%)", fill=NULL) +    
        scale_y_continuous(breaks = c(seq(0, ceiling(range(pyKVFinder$efficiency*100)[2]), by = 10)), expand = c(0,0), limits = c(0, ceiling(range(pyKVFinder$efficiency*100)[2]) + 1)) +
        theme(legend.position="none",
            axis.text.x = element_text(size = 20), 
            axis.text.y = element_text(size = 20, hjust = 1),
            legend.text = element_text(size = 20)) +
        theme(axis.title.y = element_text(angle = 90, 
                                        hjust = 0.5, 
                                        size = 20),
            axis.title.x = element_text(size = 20),
            legend.title = element_text(size = 20))

svg('plots/pyKVFinder/efficiency/boxplot_efficiency_vs_nthreads.svg', width = 16, height = 9)
p
dev.off()

png('plots/pyKVFinder/efficiency/boxplot_efficiency_vs_nthreads.png', width = 4400, height = 2475, res = 300)
p
dev.off()

###
### [==> Plotting: plots/pyKVFinder/efficiency/lineplot_efficiency_vs_nthreads.svg)
###   Description: Lineplot with speedup against number of threads
###

p <- ggplot(pyKVFinder, aes(x = nthreads, y = efficiency*100)) +
    stat_summary(fun = mean, geom = "line", color = 4, linetype = 1, size = 1) +
    stat_summary(fun = mean, fun.min = function(x) mean(x) - sd(x), fun.max = function(x) mean(x) + sd(x), geom = "errorbar", width=0.2) +
    stat_summary(fun = mean, geom = "point", fill = "white", shape = 21, size = 2) +
    theme_bw() + 
    labs(x="OpenMP threads", y="Efficiency (%)", fill=NULL) +    
    scale_y_continuous(breaks = c(seq(0, ceiling(range(pyKVFinder$efficiency*100)[2]), by = 10)), expand = c(0,0), limits = c(0, ceiling(range(pyKVFinder$efficiency*100)[2]) + 1)) +
    scale_x_continuous(breaks = c(seq(1, max(unique(nthreads)), by=1)), limits=c(0.8, max(unique(nthreads)) + .5), expand=c(0,0)) +
    theme(legend.position="none",
        axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20, hjust = 1),
        legend.text = element_text(size = 20)) +
    theme(axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        legend.title = element_text(size = 20))

svg('plots/pyKVFinder/efficiency/lineplot_efficiency_vs_nthreads.svg', width = 16, height = 9)
p
dev.off()

png('plots/pyKVFinder/efficiency/lineplot_efficiency_vs_nthreads.png', width = 4400, height = 2475, res = 300)
p
dev.off()

###
### [==> pyKVFinder --depth analysis: time, speedup, efficiency
###

dir.create('plots/pyKVFinder-depth')
dir.create('plots/pyKVFinder-depth/time')

pyKVFinder.depth = time[time$software == 'pyKVFinder.depth',]
pyKVFinder.depth$speedup = pyKVFinder$average_time / pyKVFinder.depth$average_time
pyKVFinder.depth$efficiency = pyKVFinder.depth$speedup / pyKVFinder.depth$nthreads

###
### [==> Plotting: plots/pyKVFinder-depth/time/lineplot_time_vs_atoms.svg)
###   Description: Lineplots with pyKVFinder-depth time against number of atoms for different number of threads
###

p <- ggplot(pyKVFinder.depth, aes(x=natoms, y=average_time)) +

    geom_jitter(aes(group=as.factor(nthreads), color=as.factor(nthreads)), alpha=0.4, size=0.5) +

    geom_smooth(aes(group=as.factor(nthreads), color=as.factor(nthreads)), method='glm', se=FALSE, fullrange=TRUE) +

    geom_smooth(color='black', method='glm', se=FALSE, fullrange=TRUE) +

    stat_poly_eq(formula = y ~ x, aes(color = as.factor(nthreads), label = paste(..eq.label.., ..rr.label.., sep = "~~~")), parse = TRUE, size=5) +
    
    theme_bw() +
    labs(x = "Number of atoms", y = "Time (s)", color = "OpenMP threads") +
    scale_y_continuous(breaks = c(seq(0, ceiling(range(pyKVFinder.depth$average_time)[2] / 1) * 1, by = 1)), expand = c(0,0), limits = c(0, ceiling(range(pyKVFinder.depth$average_time)[2] / 1) * 1)) +
    scale_x_continuous(breaks = c(seq(0, ceiling(range(pyKVFinder.depth$natoms)[2] / 100) * 100, by = 1000)), expand = c(0,0), limits = c(0, ceiling(range(pyKVFinder.depth$natoms)[2] / 100) * 100 + 100)) +
    scale_color_manual(values = c(brewer.pal(n=replicates, name = "Paired"), "black"), limits = c(unique(as.character(pyKVFinder.depth$nthreads)), "Average")) +

    theme(axis.text.x = element_text(size = 20, angle=45, hjust = 1), 
        axis.text.y = element_text(size = 20, hjust = 1),
        legend.text = element_text(size = 20),
        strip.text = element_text(size=20)) +
    theme(axis.title.y = element_text(angle = 90, hjust = 0.5, size = 20),
        axis.title.x = element_text(size = 20),
        legend.title = element_text(size = 20))

svg('plots/pyKVFinder-depth/time/lineplot_time_vs_atoms.svg', width = 16, height = 9)
p
dev.off()

png('plots/pyKVFinder-depth/time/lineplot_time_vs_atoms.png', width = 4400, height = 2475, res = 300)
p
dev.off()

###
### [==> Plotting: plots/pyKVFinder-depth/speedup/boxplot_speedup_vs_nthreads.svg)
###   Description: Boxplot with speedup against number of threads
###

dir.create('plots/pyKVFinder-depth/speedup')

# Boxplot grouped by software: Time x Number of threads
p <- ggplot(pyKVFinder.depth, aes(x=as.factor(nthreads), y=speedup, fill=as.factor(nthreads))) +
        geom_boxplot() +
        theme_bw() +
        labs(x="Number of OpenMP threads", y="Speedup", fill=NULL) +    
        scale_y_continuous(breaks = c(seq(0, ceiling(range(pyKVFinder.depth$speedup)[2]), by = 0.5)), expand = c(0,0), limits = c(0, ceiling(range(pyKVFinder.depth$speedup)[2]) + 0.1)) +
        theme(legend.position="none",
            axis.text.x = element_text(size = 20), 
            axis.text.y = element_text(size = 20, hjust = 1),
            legend.text = element_text(size = 20)) +
        theme(axis.title.y = element_text(angle = 90, 
                                        hjust = 0.5, 
                                        size = 20),
            axis.title.x = element_text(size = 20),
            legend.title = element_text(size = 20))

svg('plots/pyKVFinder-depth/speedup/boxplot_speedup_vs_nthreads.svg', width = 16, height = 9)
p
dev.off()

png('plots/pyKVFinder-depth/speedup/boxplot_speedup_vs_nthreads.png', width = 4400, height = 2475, res = 300)
p
dev.off()

###
### [==> Plotting: plots/pyKVFinder-depth/speedup/lineplot_speedup_vs_nthreads.svg)
###   Description: Lineplot with speedup against number of threads
###

p <- ggplot(pyKVFinder.depth, aes(x = nthreads, y = speedup)) +
    stat_summary(fun = mean, geom = "line", color = 4, linetype = 1, size = 1) +
    stat_summary(fun = mean, fun.min = function(x) mean(x) - sd(x), fun.max = function(x) mean(x) + sd(x), geom = "errorbar", width=0.2) +
    stat_summary(fun = mean, geom = "point", fill = "white", shape = 21, size = 2) +
    theme_bw() + 
    labs(x="OpenMP threads", y="Speedup", fill=NULL) +    
    scale_y_continuous(breaks = c(seq(0, ceiling(range(pyKVFinder.depth$speedup)[2] * 1.5), by = 0.1)), expand = c(0,0), limits = c(0, ceiling(range(pyKVFinder.depth$speedup)[2]) + 0.1)) +
    scale_x_continuous(breaks = c(seq(1, max(unique(nthreads)), by=1)), limits=c(0.8, max(unique(nthreads)) + .5), expand=c(0,0)) +
    theme(legend.position="none",
        axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20, hjust = 1),
        legend.text = element_text(size = 20)) +
    theme(axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        legend.title = element_text(size = 20))

svg('plots/pyKVFinder-depth/speedup/lineplot_speedup_vs_nthreads.svg', width = 16, height = 9)
p
dev.off()

png('plots/pyKVFinder-depth/speedup/lineplot_speedup_vs_nthreads.png', width = 4400, height = 2475, res = 300)
p
dev.off()

###
### [==> Plotting: plots/pyKVFinder-depth/speedup/boxplot_efficiency_vs_nthreads.svg)
###   Description: Boxplot with efficiency against number of threads
###

dir.create('plots/pyKVFinder-depth/efficiency')

# Boxplot grouped by software: Time x Number of threads
p <- ggplot(pyKVFinder.depth, aes(x=as.factor(nthreads), y=efficiency*100, fill=as.factor(nthreads))) +
        geom_boxplot() +
        theme_bw() +
        labs(x="Number of OpenMP threads", y="Efficiency (%)", fill=NULL) +    
        scale_y_continuous(breaks = c(seq(0, ceiling(range(pyKVFinder.depth$efficiency*100)[2]), by = 10)), expand = c(0,0), limits = c(0, ceiling(range(pyKVFinder.depth$efficiency*100)[2]) + 1)) +
        theme(legend.position="none",
            axis.text.x = element_text(size = 20), 
            axis.text.y = element_text(size = 20, hjust = 1),
            legend.text = element_text(size = 20)) +
        theme(axis.title.y = element_text(angle = 90, 
                                        hjust = 0.5, 
                                        size = 20),
            axis.title.x = element_text(size = 20),
            legend.title = element_text(size = 20))

svg('plots/pyKVFinder-depth/efficiency/boxplot_efficiency_vs_nthreads.svg', width = 16, height = 9)
p
dev.off()

png('plots/pyKVFinder-depth/efficiency/boxplot_efficiency_vs_nthreads.png', width = 4400, height = 2475, res = 300)
p
dev.off()

###
### [==> Plotting: plots/pyKVFinder-depth/efficiency/lineplot_efficiency_vs_nthreads.svg)
###   Description: Lineplot with speedup against number of threads
###

p <- ggplot(pyKVFinder.depth, aes(x = nthreads, y = efficiency*100)) +
    stat_summary(fun = mean, geom = "line", color = 4, linetype = 1, size = 1) +
    stat_summary(fun = mean, fun.min = function(x) mean(x) - sd(x), fun.max = function(x) mean(x) + sd(x), geom = "errorbar", width=0.2) +
    stat_summary(fun = mean, geom = "point", fill = "white", shape = 21, size = 2) +
    theme_bw() + 
    labs(x="OpenMP threads", y="Efficiency (%)", fill=NULL) +    
    scale_y_continuous(breaks = c(seq(0, ceiling(range(pyKVFinder.depth$efficiency*100)[2]), by = 10)), expand = c(0,0), limits = c(0, ceiling(range(pyKVFinder.depth$efficiency*100)[2]) + 1)) +
    scale_x_continuous(breaks = c(seq(1, max(unique(nthreads)), by=1)), limits=c(0.8, max(unique(nthreads)) + .5), expand=c(0,0)) +
    theme(legend.position="none",
        axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20, hjust = 1),
        legend.text = element_text(size = 20)) +
    theme(axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        legend.title = element_text(size = 20))

svg('plots/pyKVFinder-depth/efficiency/lineplot_efficiency_vs_nthreads.svg', width = 16, height = 9)
p
dev.off()

png('plots/pyKVFinder-depth/efficiency/lineplot_efficiency_vs_nthreads.png', width = 4400, height = 2475, res = 300)
p
dev.off()

###
### [==> pyKVFinder (hydropathy) analysis: time, speedup, efficiency
###

dir.create('plots/pyKVFinder-hydropathy')
dir.create('plots/pyKVFinder-hydropathy/time')

pyKVFinder.hydropathy = time[time$software == 'pyKVFinder.hydropathy',]
pyKVFinder.hydropathy$speedup = pyKVFinder$average_time / pyKVFinder.hydropathy$average_time
pyKVFinder.hydropathy$efficiency = pyKVFinder.hydropathy$speedup / pyKVFinder.hydropathy$nthreads

###
### [==> Plotting: plots/pyKVFinder-hydropathy/time/lineplot_time_vs_atoms.svg)
###   Description: Lineplots with pyKVFinder-hydropathy time against number of atoms for different number of threads
###

p <- ggplot(pyKVFinder.hydropathy, aes(x=natoms, y=average_time)) +

    geom_jitter(aes(group=as.factor(nthreads), color=as.factor(nthreads)), alpha=0.4, size=0.5) +

    geom_smooth(aes(group=as.factor(nthreads), color=as.factor(nthreads)), method='glm', se=FALSE, fullrange=TRUE) +

    geom_smooth(color='black', method='glm', se=FALSE, fullrange=TRUE) +

    stat_poly_eq(formula = y ~ x, aes(color = as.factor(nthreads), label = paste(..eq.label.., ..rr.label.., sep = "~~~")), parse = TRUE, size=5) +
    
    theme_bw() +
    labs(x = "Number of atoms", y = "Time (s)", color = "OpenMP threads") +
    scale_y_continuous(breaks = c(seq(0, ceiling(range(pyKVFinder.hydropathy$average_time)[2] / 1) * 1, by = 1)), expand = c(0,0), limits = c(0, ceiling(range(pyKVFinder.hydropathy$average_time)[2] / 1) * 1)) +
    scale_x_continuous(breaks = c(seq(0, ceiling(range(pyKVFinder.hydropathy$natoms)[2] / 100) * 100, by = 1000)), expand = c(0,0), limits = c(0, ceiling(range(pyKVFinder.hydropathy$natoms)[2] / 100) * 100 + 100)) +
    scale_color_manual(values = c(brewer.pal(n=replicates, name = "Paired"), "black"), limits = c(unique(as.character(pyKVFinder.hydropathy$nthreads)), "Average")) +

    theme(axis.text.x = element_text(size = 20, angle=45, hjust = 1), 
        axis.text.y = element_text(size = 20, hjust = 1),
        legend.text = element_text(size = 20),
        strip.text = element_text(size=20)) +
    theme(axis.title.y = element_text(angle = 90, hjust = 0.5, size = 20),
        axis.title.x = element_text(size = 20),
        legend.title = element_text(size = 20))

svg('plots/pyKVFinder-hydropathy/time/lineplot_time_vs_atoms.svg', width = 16, height = 9)
p
dev.off()

png('plots/pyKVFinder-hydropathy/time/lineplot_time_vs_atoms.png', width = 4400, height = 2475, res = 300)
p
dev.off()

###
### [==> Plotting: plots/pyKVFinder-hydropathy/speedup/boxplot_speedup_vs_nthreads.svg)
###   Description: Boxplot with speedup against number of threads
###

dir.create('plots/pyKVFinder-hydropathy/speedup')

# Boxplot grouped by software: Time x Number of threads
p <- ggplot(pyKVFinder.hydropathy, aes(x=as.factor(nthreads), y=speedup, fill=as.factor(nthreads))) +
        geom_boxplot() +
        theme_bw() +
        labs(x="Number of OpenMP threads", y="Speedup", fill=NULL) +    
        # scale_y_continuous(breaks = c(seq(0, ceiling(range(pyKVFinder.hydropathy$speedup)[2]), by = 0.5)), expand = c(0,0), limits = c(0, ceiling(range(pyKVFinder.hydropathy$speedup)[2]) + 0.1)) +
        theme(legend.position="none",
            axis.text.x = element_text(size = 20), 
            axis.text.y = element_text(size = 20, hjust = 1),
            legend.text = element_text(size = 20)) +
        theme(axis.title.y = element_text(angle = 90, 
                                        hjust = 0.5, 
                                        size = 20),
            axis.title.x = element_text(size = 20),
            legend.title = element_text(size = 20))

svg('plots/pyKVFinder-hydropathy/speedup/boxplot_speedup_vs_nthreads.svg', width = 16, height = 9)
p
dev.off()

png('plots/pyKVFinder-hydropathy/speedup/boxplot_speedup_vs_nthreads.png', width = 4400, height = 2475, res = 300)
p
dev.off()

###
### [==> Plotting: plots/pyKVFinder-hydropathy/speedup/lineplot_speedup_vs_nthreads.svg)
###   Description: Lineplot with speedup against number of threads
###

p <- ggplot(pyKVFinder.hydropathy, aes(x = nthreads, y = speedup)) +
    stat_summary(fun = mean, geom = "line", color = 4, linetype = 1, size = 1) +
    stat_summary(fun = mean, fun.min = function(x) mean(x) - sd(x), fun.max = function(x) mean(x) + sd(x), geom = "errorbar", width=0.2) +
    stat_summary(fun = mean, geom = "point", fill = "white", shape = 21, size = 2) +
    theme_bw() + 
    labs(x="OpenMP threads", y="Speedup", fill=NULL) +    
    scale_y_continuous(breaks = c(seq(0, ceiling(range(pyKVFinder.hydropathy$speedup)[2]), by = 0.1)), expand = c(0,0), limits = c(0, ceiling(range(pyKVFinder.hydropathy$speedup)[2]) + 0.1)) +
    scale_x_continuous(breaks = c(seq(1, max(unique(nthreads)), by=1)), limits=c(0.8, max(unique(nthreads)) + .5), expand=c(0,0)) +
    theme(legend.position="none",
        axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20, hjust = 1),
        legend.text = element_text(size = 20)) +
    theme(axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        legend.title = element_text(size = 20))

svg('plots/pyKVFinder-hydropathy/speedup/lineplot_speedup_vs_nthreads.svg', width = 16, height = 9)
p
dev.off()

png('plots/pyKVFinder-hydropathy/speedup/lineplot_speedup_vs_nthreads.png', width = 4400, height = 2475, res = 300)
p
dev.off()

###
### [==> Plotting: plots/pyKVFinder-hydropathy/speedup/boxplot_efficiency_vs_nthreads.svg)
###   Description: Boxplot with efficiency against number of threads
###

dir.create('plots/pyKVFinder-hydropathy/efficiency')

# Boxplot grouped by software: Time x Number of threads
p <- ggplot(pyKVFinder.hydropathy, aes(x=as.factor(nthreads), y=efficiency*100, fill=as.factor(nthreads))) +
        geom_boxplot() +
        theme_bw() +
        labs(x="Number of OpenMP threads", y="Efficiency (%)", fill=NULL) +    
        scale_y_continuous(breaks = c(seq(0, ceiling(range(pyKVFinder.hydropathy$efficiency*100)[2]), by = 10)), expand = c(0,0), limits = c(0, ceiling(range(pyKVFinder.hydropathy$efficiency*100)[2]) + 1)) +
        theme(legend.position="none",
            axis.text.x = element_text(size = 20), 
            axis.text.y = element_text(size = 20, hjust = 1),
            legend.text = element_text(size = 20)) +
        theme(axis.title.y = element_text(angle = 90, 
                                        hjust = 0.5, 
                                        size = 20),
            axis.title.x = element_text(size = 20),
            legend.title = element_text(size = 20))

svg('plots/pyKVFinder-hydropathy/efficiency/boxplot_efficiency_vs_nthreads.svg', width = 16, height = 9)
p
dev.off()

png('plots/pyKVFinder-hydropathy/efficiency/boxplot_efficiency_vs_nthreads.png', width = 4400, height = 2475, res = 300)
p
dev.off()

###
### [==> Plotting: plots/pyKVFinder-hydropathy/efficiency/lineplot_efficiency_vs_nthreads.svg)
###   Description: Lineplot with speedup against number of threads
###

p <- ggplot(pyKVFinder.hydropathy, aes(x = nthreads, y = efficiency*100)) +
    stat_summary(fun = mean, geom = "line", color = 4, linetype = 1, size = 1) +
    stat_summary(fun = mean, fun.min = function(x) mean(x) - sd(x), fun.max = function(x) mean(x) + sd(x), geom = "errorbar", width=0.2) +
    stat_summary(fun = mean, geom = "point", fill = "white", shape = 21, size = 2) +
    theme_bw() + 
    labs(x="OpenMP threads", y="Efficiency (%)", fill=NULL) +    
    scale_y_continuous(breaks = c(seq(0, ceiling(range(pyKVFinder.hydropathy$efficiency*100)[2]), by = 10)), expand = c(0,0), limits = c(0, ceiling(range(pyKVFinder.hydropathy$efficiency*100)[2]) + 1)) +
    scale_x_continuous(breaks = c(seq(1, max(unique(nthreads)), by=1)), limits=c(0.8, max(unique(nthreads)) + .5), expand=c(0,0)) +
    theme(legend.position="none",
        axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20, hjust = 1),
        legend.text = element_text(size = 20)) +
    theme(axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        legend.title = element_text(size = 20))

svg('plots/pyKVFinder-hydropathy/efficiency/lineplot_efficiency_vs_nthreads.svg', width = 16, height = 9)
p
dev.off()

png('plots/pyKVFinder-hydropathy/efficiency/lineplot_efficiency_vs_nthreads.png', width = 4400, height = 2475, res = 300)
p
dev.off()
