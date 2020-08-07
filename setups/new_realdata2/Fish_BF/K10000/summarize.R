library(ggplot2)
library(Rmisc) 
library(Sleuth2)
library(ggpubr)


PGSBF1 <- data.frame(read.csv("results1/results.csv"))

PGSBF2 <- data.frame(read.csv("results2/results.csv"))

PGSBF3 <- data.frame(read.csv("results3/results.csv"))

PGSBF4 <- data.frame(read.csv("results4/results.csv"))

PGSBF5 <- data.frame(read.csv("results5/results.csv"))

PGSBF6 <- data.frame(read.csv("results6/results.csv"))

PGSBF7 <- data.frame(read.csv("results7/results.csv"))

PGSBF8 <- data.frame(read.csv("results8/results.csv"))

PGSBF9 <- data.frame(read.csv("results9/results.csv"))

PGSBF10 <- data.frame(read.csv("results10/results.csv"))

PGSBF11 <- data.frame(read.csv("results11/results.csv"))

PGSBF12 <- data.frame(read.csv("results12/results.csv"))


PGSBF <- rbind(PGSBF1, PGSBF2, PGSBF3, PGSBF4, PGSBF5, PGSBF6, PGSBF7, PGSBF8, PGSBF9, PGSBF10, PGSBF11, PGSBF12)

index <- which(PGSBF$Metric == 'ConsensusLogLL')
CLL <- PGSBF$Value[index]
mean(CLL)
sd(CLL)
