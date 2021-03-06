library(ggplot2)
library(Rmisc) 
library(Sleuth2)
library(ggpubr)


CSMC1 <- data.frame(read.csv("CSMC/1/results.csv"))
CSMC1$gen.nTaxa <- read.csv("CSMC/1/options.map", sep = '')[57,2]
CSMC2 <- data.frame(read.csv("CSMC/2/results.csv"))
CSMC2$gen.nTaxa <- read.csv("CSMC/2/options.map", sep = '')[57,2]
CSMC3 <- data.frame(read.csv("CSMC/3/results.csv"))
CSMC3$gen.nTaxa <- read.csv("CSMC/3/options.map", sep = '')[57,2]
CSMC4 <- data.frame(read.csv("CSMC/4/results.csv"))
CSMC4$gen.nTaxa <- read.csv("CSMC/4/options.map", sep = '')[57,2]

CSMC5 <- data.frame(read.csv("CSMC/5/results.csv"))
CSMC5$gen.nTaxa <- read.csv("CSMC/5/options.map", sep = '')[57,2]
CSMC6 <- data.frame(read.csv("CSMC/6/results.csv"))
CSMC6$gen.nTaxa <- read.csv("CSMC/6/options.map", sep = '')[57,2]
CSMC7 <- data.frame(read.csv("CSMC/7/results.csv"))
CSMC7$gen.nTaxa <- read.csv("CSMC/7/options.map", sep = '')[57,2]
CSMC8 <- data.frame(read.csv("CSMC/8/results.csv"))
CSMC8$gen.nTaxa <- read.csv("CSMC/8/options.map", sep = '')[57,2]
CSMC9 <- data.frame(read.csv("CSMC/9/results.csv"))
CSMC9$gen.nTaxa <- read.csv("CSMC/9/options.map", sep = '')[57,2]
CSMC10 <- data.frame(read.csv("CSMC/10/results.csv"))
CSMC10$gen.nTaxa <- read.csv("CSMC/10/options.map", sep = '')[57,2]
CSMC11 <- data.frame(read.csv("CSMC/11/results.csv"))
CSMC11$gen.nTaxa <- read.csv("CSMC/11/options.map", sep = '')[57,2]
CSMC12 <- data.frame(read.csv("CSMC/12/results.csv"))
CSMC12$gen.nTaxa <- read.csv("CSMC/12/options.map", sep = '')[57,2]

CSMC <- rbind(CSMC1, CSMC2, CSMC3, CSMC4, CSMC5, CSMC6, CSMC7, CSMC8, CSMC9, CSMC10, CSMC11, CSMC12)

CSMC$IterScale[which(CSMC$IterScale == 1.5)] <- 1
CSMC$IterScale[which(CSMC$IterScale == 4.5)] <- 3
CSMC$IterScale[which(CSMC$IterScale == 15)] <- 10
CSMC$IterScale[which(CSMC$IterScale == 45)] <- 30


CSMCBF1 <- data.frame(read.csv("CSMCBF/results1/results.csv"))
CSMCBF1$gen.nTaxa <- read.csv("CSMCBF/results1/options.map", sep = '')[57,2]
CSMCBF2 <- data.frame(read.csv("CSMCBF/results2/results.csv"))
CSMCBF2$gen.nTaxa <- read.csv("CSMCBF/results2/options.map", sep = '')[57,2]
CSMCBF3 <- data.frame(read.csv("CSMCBF/results3/results.csv"))
CSMCBF3$gen.nTaxa <- read.csv("CSMCBF/results3/options.map", sep = '')[57,2]
CSMCBF4 <- data.frame(read.csv("CSMCBF/results4/results.csv"))
CSMCBF4$gen.nTaxa <- read.csv("CSMCBF/results4/options.map", sep = '')[57,2]

CSMCBF5 <- data.frame(read.csv("CSMCBF/results5/results.csv"))
CSMCBF5$gen.nTaxa <- read.csv("CSMCBF/results5/options.map", sep = '')[57,2]
CSMCBF6 <- data.frame(read.csv("CSMCBF/results6/results.csv"))
CSMCBF6$gen.nTaxa <- read.csv("CSMCBF/results6/options.map", sep = '')[57,2]
CSMCBF7 <- data.frame(read.csv("CSMCBF/results7/results.csv"))
CSMCBF7$gen.nTaxa <- read.csv("CSMCBF/results7/options.map", sep = '')[57,2]
CSMCBF8 <- data.frame(read.csv("CSMCBF/results8/results.csv"))
CSMCBF8$gen.nTaxa <- read.csv("CSMCBF/results8/options.map", sep = '')[57,2]
CSMCBF9 <- data.frame(read.csv("CSMCBF/results9/results.csv"))
CSMCBF9$gen.nTaxa <- read.csv("CSMCBF/results9/options.map", sep = '')[57,2]
CSMCBF10 <- data.frame(read.csv("CSMCBF/results10/results.csv"))
CSMCBF10$gen.nTaxa <- read.csv("CSMCBF/results10/options.map", sep = '')[57,2]
CSMCBF11 <- data.frame(read.csv("CSMCBF/results11/results.csv"))
CSMCBF11$gen.nTaxa <- read.csv("CSMCBF/results11/options.map", sep = '')[57,2]
CSMCBF12 <- data.frame(read.csv("CSMCBF/results12/results.csv"))
CSMCBF12$gen.nTaxa <- read.csv("CSMCBF/results12/options.map", sep = '')[57,2]

CSMCBF <- rbind(CSMCBF1, CSMCBF2, CSMCBF3, CSMCBF4, CSMCBF5, CSMCBF6, CSMCBF7, CSMCBF8, CSMCBF9, CSMCBF10, CSMCBF11, CSMCBF12)



CSMCRM <- CSMCBF
CSMC$method <- 'CSMC'
CSMCRM$method <- 'CSMC-RDouP'

results <- rbind(CSMC, CSMCRM)
results_RF <- results[which(results$Metric == 'RobinsonFouldsMetric'),]

unique(results_RF$gen.nTaxa)

results_RF10 <- results_RF[which(results_RF$gen.nTaxa == 10),]
results_RF40 <- results_RF[which(results_RF$gen.nTaxa == 40),]
results_RF100 <- results_RF[which(results_RF$gen.nTaxa == 100),]

results_RF_CSMC <- results_RF[which(results_RF$method == 'CSMC'),]
results_RF_CSMCRM <- results_RF[which(results_RF$method == 'CSMC-RDouP'),]

results_RF10_CSMC <- results_RF_CSMC[which(results_RF_CSMC$gen.nTaxa == 10),]
results_RF40_CSMC <- results_RF_CSMC[which(results_RF_CSMC$gen.nTaxa == 40),]
results_RF100_CSMC <- results_RF_CSMC[which(results_RF_CSMC$gen.nTaxa == 100),]

results_RF10_CSMCRM <- results_RF_CSMCRM[which(results_RF_CSMCRM$gen.nTaxa == 10),]
results_RF40_CSMCRM <- results_RF_CSMCRM[which(results_RF_CSMCRM$gen.nTaxa == 40),]
results_RF100_CSMCRM <- results_RF_CSMCRM[which(results_RF_CSMCRM$gen.nTaxa == 100),]



RFsummary10_CSMCRM <- summarySE(results_RF10_CSMCRM, measurevar="Value", groupvars=c("IterScale","method"))
RFsummary10_CSMCRM$IterScale = as.factor(RFsummary10_CSMCRM$IterScale)
RFsummary10_CSMC <- summarySE(results_RF10_CSMC, measurevar="Value", groupvars=c("IterScale","method"))
RFsummary10_CSMC$IterScale = as.factor(RFsummary10_CSMC$IterScale)
RFsummary10 <- rbind(RFsummary10_CSMCRM, RFsummary10_CSMC)
p2 <- ggplot(RFsummary10, aes(x = IterScale, y=Value, group = method, colour=method))



RFsummary40_CSMCRM <- summarySE(results_RF40_CSMCRM, measurevar="Value", groupvars=c("IterScale","method"))
RFsummary40_CSMCRM$IterScale = as.factor(RFsummary40_CSMCRM$IterScale)
RFsummary40_CSMC <- summarySE(results_RF40_CSMC, measurevar="Value", groupvars=c("IterScale","method"))
RFsummary40_CSMC$IterScale = as.factor(RFsummary40_CSMC$IterScale)
RFsummary40 <- rbind(RFsummary40_CSMCRM, RFsummary40_CSMC)
p3 <- ggplot(RFsummary40, aes(x = IterScale, y=Value, group = method, colour=method))


RFsummary100_CSMCRM <- summarySE(results_RF100_CSMCRM, measurevar="Value", groupvars=c("IterScale","method"))
RFsummary100_CSMCRM$IterScale = as.factor(RFsummary100_CSMCRM$IterScale)
RFsummary100_CSMC <- summarySE(results_RF100_CSMC, measurevar="Value", groupvars=c("IterScale","method"))
RFsummary100_CSMC$IterScale = as.factor(RFsummary100_CSMC$IterScale)
RFsummary100 <- rbind(RFsummary100_CSMCRM, RFsummary100_CSMC)
p4 <- ggplot(RFsummary100, aes(x = IterScale, y=Value, group = method, colour=method))


pd <- position_dodge(.2)
gname = c("CSMCRF.eps",sep="")  
postscript(gname,width=8,height=2.5,horizontal = FALSE, onefile = FALSE, paper = "special")
par(mfrow=c(1,1),oma=c(0.2,1.5,0.2,1.5),mar=c(3,2,0.2,2),cex.axis=1,las=1,mgp=c(1,0.5,0),adj=0.5)
ggarrange(  p2 + geom_line(position=pd) +
            geom_point(position=pd) +  geom_errorbar(aes(ymin=Value-ci, ymax=Value+ci), width=.3, position=pd) +
            rremove("ylab") + xlab('#taxa = 10')+ theme_bw()+ylab(NULL)  +
            scale_color_manual(values=c('#999999','#000000')),
          p3 + geom_line(position=pd) +
            geom_point(position=pd) +  geom_errorbar(aes(ymin=Value-ci, ymax=Value+ci), width=.3, position=pd) +
            rremove("ylab") + xlab('#taxa = 40')+ theme_bw()+ylab(NULL)  +
            scale_color_manual(values=c('#999999','#000000')),
          p4 + geom_line(position=pd) +
            geom_point(position=pd) +  geom_errorbar(aes(ymin=Value-ci, ymax=Value+ci), width=.3, position=pd) +
            rremove("ylab") + xlab('#taxa = 100')+ theme_bw()+ylab(NULL)  +
            scale_color_manual(values=c('#999999','#000000')),
          ncol = 3, nrow = 1, common.legend = TRUE)

dev.off()





