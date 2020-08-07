library(ggplot2)
library(Rmisc) 
library(Sleuth2)
library(ggpubr)


IPGS1 <- data.frame(read.csv("results1/results.csv"))

IPGS2 <- data.frame(read.csv("results2/results.csv"))

IPGS3 <- data.frame(read.csv("results3/results.csv"))

IPGS4 <- data.frame(read.csv("results4/results.csv"))

IPGS5 <- data.frame(read.csv("results5/results.csv"))

IPGS6 <- data.frame(read.csv("results6/results.csv"))

IPGS7 <- data.frame(read.csv("results7/results.csv"))

IPGS8 <- data.frame(read.csv("results8/results.csv"))

IPGS9 <- data.frame(read.csv("results9/results.csv"))

IPGS10 <- data.frame(read.csv("results10/results.csv"))

IPGS11 <- data.frame(read.csv("results11/results.csv"))

IPGS12 <- data.frame(read.csv("results12/results.csv"))


IPGS <- rbind(IPGS1, IPGS2, IPGS3, IPGS4, IPGS5, IPGS6, IPGS7, IPGS8, IPGS9, IPGS10, IPGS11, IPGS12)

index <- which(IPGS$Metric == 'ConsensusLogLL')
CLL <- IPGS$Value[index]
mean(CLL)
sd(CLL)
