library(ggplot2)
library(Rmisc) 
library(Sleuth2)
library(ggpubr)


core1 <- data.frame(read.csv("1/results/results.csv"))
core1$cores = 1

core2 <- data.frame(read.csv("2/results/results.csv"))
core2$cores = 2

core3 <- data.frame(read.csv("3/results/results.csv"))
core3$cores = 3

core4 <- data.frame(read.csv("4/results/results.csv"))
core4$cores = 4

core5 <- data.frame(read.csv("5/results/results.csv"))
core5$cores = 5

core6 <- data.frame(read.csv("6/results/results.csv"))
core6$cores = 6


cores = rbind(core1, core2, core3, core4, core5, core6)
cores$core = as.factor(cores$core)

coressummary <- summarySE(cores, measurevar="Time", groupvars=c("core"))


gname = c("cores.eps",sep="")  
postscript(gname,width=6,height=3,horizontal = FALSE, onefile = FALSE, paper = "special")
par(mfrow=c(1,1),oma=c(0.2,1.5,0.2,1.5),mar=c(3,2,0.2,2),cex.axis=1,las=1,mgp=c(1,0.5,0),adj=0.5)

ggplot(coressummary, aes(x=core, y=Time)) + 
  geom_errorbar(aes(ymin=Time-ci, ymax=Time+ci), width=.1) +
  geom_line() +
  geom_point()+ theme_bw() + ylab("Time (milliseconds)")+ xlab("Number of Threads")
dev.off()

