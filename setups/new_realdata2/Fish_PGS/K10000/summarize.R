library(ggplot2)
library(Rmisc) 
library(Sleuth2)
library(ggpubr)


PGS1 <- data.frame(read.csv("results1/results.csv"))

PGS2 <- data.frame(read.csv("results2/results.csv"))

PGS3 <- data.frame(read.csv("results3/results.csv"))

#PGS4 <- data.frame(read.csv("results4/results.csv"))

PGS5 <- data.frame(read.csv("results5/results.csv"))

PGS6 <- data.frame(read.csv("results6/results.csv"))

PGS7 <- data.frame(read.csv("results7/results.csv"))

PGS8 <- data.frame(read.csv("results8/results.csv"))

#PGS9 <- data.frame(read.csv("results9/results.csv"))

PGS10 <- data.frame(read.csv("results10/results.csv"))

#PGS11 <- data.frame(read.csv("results11/results.csv"))

PGS12 <- data.frame(read.csv("results12/results.csv"))


#PGS <- rbind(PGS1, PGS2, PGS3, PGS4, PGS5, PGS6, PGS7, PGS8, PGS9, PGS10, PGS11, PGS12)
PGS <- rbind(PGS1, PGS2, PGS3, PGS5, PGS6, PGS7, PGS8, PGS10, PGS12)

index <- which(PGS$Metric == 'ConsensusLogLL')
CLL <- PGS$Value[index]
mean(CLL)
sd(CLL)
