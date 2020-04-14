library(ggplot2)
library(Rmisc) 
library(Sleuth2)
library(ggpubr)


CSMC <- data.frame(read.csv("CSMC.csv"))
CSMC$IterScale[which(CSMC$IterScale == 1.5)] <- 1
CSMC$IterScale[which(CSMC$IterScale == 4.5)] <- 3
CSMC$IterScale[which(CSMC$IterScale == 15)] <- 10
CSMC$IterScale[which(CSMC$IterScale == 45)] <- 30
CSMCRM <- data.frame(read.csv("CSMCRM.csv"))
CSMC$method <- 'CSMC'
CSMCRM$method <- 'CSMC-RDouP'

results <- rbind(CSMC, CSMCRM)
results_RF <- results[which(results$Metric == 'RobinsonFouldsMetric'),]

unique(results_RF$gen.nTaxa)
results_RF5 <- results_RF[which(results_RF$gen.nTaxa == 5),]
results_RF10 <- results_RF[which(results_RF$gen.nTaxa == 10),]
results_RF20 <- results_RF[which(results_RF$gen.nTaxa == 20),]
results_RF40 <- results_RF[which(results_RF$gen.nTaxa == 40),]

results_RF_CSMC <- results_RF[which(results_RF$method == 'CSMC'),]
results_RF_CSMCRM <- results_RF[which(results_RF$method == 'CSMC-RDouP'),]

results_RF5_CSMC <- results_RF_CSMC[which(results_RF_CSMC$gen.nTaxa == 5),]
results_RF10_CSMC <- results_RF_CSMC[which(results_RF_CSMC$gen.nTaxa == 10),]
results_RF20_CSMC <- results_RF_CSMC[which(results_RF_CSMC$gen.nTaxa == 20),]
results_RF40_CSMC <- results_RF_CSMC[which(results_RF_CSMC$gen.nTaxa == 40),]

results_RF5_CSMCRM <- results_RF_CSMCRM[which(results_RF_CSMCRM$gen.nTaxa == 5),]
results_RF10_CSMCRM <- results_RF_CSMCRM[which(results_RF_CSMCRM$gen.nTaxa == 10),]
results_RF20_CSMCRM <- results_RF_CSMCRM[which(results_RF_CSMCRM$gen.nTaxa == 20),]
results_RF40_CSMCRM <- results_RF_CSMCRM[which(results_RF_CSMCRM$gen.nTaxa == 40),]

RFsummary5_CSMCRM <- summarySE(results_RF5_CSMCRM, measurevar="Value", groupvars=c("IterScale","method"))
RFsummary5_CSMCRM$IterScale = as.factor(RFsummary5_CSMCRM$IterScale)
RFsummary5_CSMC <- summarySE(results_RF5_CSMC, measurevar="Value", groupvars=c("IterScale","method"))
RFsummary5_CSMC$IterScale = as.factor(RFsummary5_CSMC$IterScale)
RFsummary5 <- rbind(RFsummary5_CSMCRM, RFsummary5_CSMC)
p1 <- ggplot(RFsummary5, aes(x = IterScale, y=Value, group = method, colour=method))


RFsummary10_CSMCRM <- summarySE(results_RF10_CSMCRM, measurevar="Value", groupvars=c("IterScale","method"))
RFsummary10_CSMCRM$IterScale = as.factor(RFsummary10_CSMCRM$IterScale)
RFsummary10_CSMC <- summarySE(results_RF10_CSMC, measurevar="Value", groupvars=c("IterScale","method"))
RFsummary10_CSMC$IterScale = as.factor(RFsummary10_CSMC$IterScale)
RFsummary10 <- rbind(RFsummary10_CSMCRM, RFsummary10_CSMC)
p2 <- ggplot(RFsummary10, aes(x = IterScale, y=Value, group = method, colour=method))

RFsummary20_CSMCRM <- summarySE(results_RF20_CSMCRM, measurevar="Value", groupvars=c("IterScale","method"))
RFsummary20_CSMCRM$IterScale = as.factor(RFsummary20_CSMCRM$IterScale)
RFsummary20_CSMC <- summarySE(results_RF20_CSMC, measurevar="Value", groupvars=c("IterScale","method"))
RFsummary20_CSMC$IterScale = as.factor(RFsummary20_CSMC$IterScale)
RFsummary20 <- rbind(RFsummary20_CSMCRM, RFsummary20_CSMC)
p3 <- ggplot(RFsummary20, aes(x = IterScale, y=Value, group = method, colour=method))

RFsummary40_CSMCRM <- summarySE(results_RF40_CSMCRM, measurevar="Value", groupvars=c("IterScale","method"))
RFsummary40_CSMCRM$IterScale = as.factor(RFsummary40_CSMCRM$IterScale)
RFsummary40_CSMC <- summarySE(results_RF40_CSMC, measurevar="Value", groupvars=c("IterScale","method"))
RFsummary40_CSMC$IterScale = as.factor(RFsummary40_CSMC$IterScale)
RFsummary40 <- rbind(RFsummary40_CSMCRM, RFsummary40_CSMC)
p4 <- ggplot(RFsummary40, aes(x = IterScale, y=Value, group = method, colour=method))


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
            rremove("ylab") + xlab('#taxa = 20')+ theme_bw()+ylab(NULL)  +
            scale_color_manual(values=c('#999999','#000000')),
          p4 + geom_line(position=pd) +
            geom_point(position=pd) +  geom_errorbar(aes(ymin=Value-ci, ymax=Value+ci), width=.3, position=pd) +
            rremove("ylab") + xlab('#taxa = 40')+ theme_bw()+ylab(NULL)  +
            scale_color_manual(values=c('#999999','#000000')),
          ncol = 3, nrow = 1, common.legend = TRUE)

dev.off()





