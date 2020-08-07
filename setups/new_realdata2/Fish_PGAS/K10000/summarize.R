library(ggplot2)
library(Rmisc) 
library(Sleuth2)
library(ggpubr)


PGAS1 <- data.frame(read.csv("results1/results.csv"))

PGAS2 <- data.frame(read.csv("results2/results.csv"))

PGAS3 <- data.frame(read.csv("results3/results.csv"))

PGAS4 <- data.frame(read.csv("results4/results.csv"))

PGAS5 <- data.frame(read.csv("results5/results.csv"))

PGAS6 <- data.frame(read.csv("results6/results.csv"))

PGAS7 <- data.frame(read.csv("results7/results.csv"))

PGAS8 <- data.frame(read.csv("results8/results.csv"))

PGAS9 <- data.frame(read.csv("results9/results.csv"))

PGAS10 <- data.frame(read.csv("results10/results.csv"))

PGAS11 <- data.frame(read.csv("results11/results.csv"))

PGAS12 <- data.frame(read.csv("results12/results.csv"))


PGAS <- rbind(PGAS1, PGAS2, PGAS3, PGAS4, PGAS5, PGAS6, PGAS7, PGAS8, PGAS9, PGAS10, PGAS11, PGAS12)

index <- which(PGAS$Metric == 'ConsensusLogLL')
CLL <- PGAS$Value[index]
mean(CLL)
sd(CLL)



