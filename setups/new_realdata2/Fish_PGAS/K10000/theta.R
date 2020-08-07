library(ggplot2)
library(Rmisc) 
library(Sleuth2)
library(ggpubr)


PGAS1 <- mean(data.frame(read.csv("results1/PGAS4K2PBF.csv"))$trans2tranv[-(1:2500)])

PGAS2 <- mean(data.frame(read.csv("results2/PGAS4K2PBF.csv"))$trans2tranv[-(1:2500)])

PGAS3 <- mean(data.frame(read.csv("results3/PGAS4K2PBF.csv"))$trans2tranv[-(1:2500)])

PGAS4 <- mean(data.frame(read.csv("results4/PGAS4K2PBF.csv"))$trans2tranv[-(1:2500)])

PGAS5 <- mean(data.frame(read.csv("results5/PGAS4K2PBF.csv"))$trans2tranv[-(1:2500)])

PGAS6 <- mean(data.frame(read.csv("results6/PGAS4K2PBF.csv"))$trans2tranv[-(1:2500)])

PGAS7 <- mean(data.frame(read.csv("results7/PGAS4K2PBF.csv"))$trans2tranv[-(1:2500)])

PGAS8 <- mean(data.frame(read.csv("results8/PGAS4K2PBF.csv"))$trans2tranv[-(1:2500)])

PGAS9 <- mean(data.frame(read.csv("results9/PGAS4K2PBF.csv"))$trans2tranv[-(1:2500)])

PGAS10 <- mean(data.frame(read.csv("results10/PGAS4K2PBF.csv"))$trans2tranv[-(1:2500)])

PGAS11 <- mean(data.frame(read.csv("results11/PGAS4K2PBF.csv"))$trans2tranv[-(1:2500)])

PGAS12 <- mean(data.frame(read.csv("results12/PGAS4K2PBF.csv"))$trans2tranv[-(1:2500)])


PGAS <- rbind(PGAS1, PGAS2, PGAS3, PGAS4, PGAS5, PGAS6, PGAS7, PGAS8, PGAS9, PGAS10, PGAS11, PGAS12)


mean(PGAS)
sd(PGAS)



