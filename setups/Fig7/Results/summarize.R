library(ggplot2)
library(Rmisc) 
library(Sleuth2)
library(ggpubr)


PGS1 <- data.frame(read.csv("results1/results.csv"))
PGS2 <- data.frame(read.csv("results2/results.csv"))
PGS3 <- data.frame(read.csv("results3/results.csv"))
#PGS4 <- data.frame(read.csv("results4/results.csv"))

#results <- rbind(PGS1, PGS2, PGS3, PGS4)
results <- rbind(PGS1, PGS2, PGS3)
results_RF <- results[which(results$Metric == 'RobinsonFouldsMetric'),]
results_RF$K <-2000000/results_RF$IterScale

RFsummary_PGSRM <- summarySE(results_RF, measurevar="Value", groupvars=c("K"))
RFsummary_PGSRM$K = as.factor(RFsummary_PGSRM$K)

results_PM <- results[which(results$Metric == 'PartitionMetric'),]
results_PM$K <-results_PM$IterScale

PMsummary_PGSRM <- summarySE(results_PM, measurevar="Value", groupvars=c("K"))
PMsummary_PGSRM$K = as.factor(PMsummary_PGSRM$K)

results_ConsensusLL <- results[which(results$Metric == 'ConsensusLogLL'),]
results_ConsensusLL$K <-2000000/results_ConsensusLL$IterScale

CLLsummary_PGSRM <- summarySE(results_ConsensusLL, measurevar="Value", groupvars=c("K"))
CLLsummary_PGSRM$K = as.factor(CLLsummary_PGSRM$K)

p3 <- ggplot(RFsummary_PGSRM, aes(x = K, y=Value, group = 1))
#p2 <- ggplot(PMsummary_PGSRM, aes(x = K, y=Value, group = 1))
p1 <- ggplot(CLLsummary_PGSRM, aes(x = K, y=Value, group = 1))


pd <- position_dodge(.2)
gname = c("PGSsquare.eps",sep="")  
postscript(gname,width=5,height=1.5,horizontal = FALSE, onefile = FALSE, paper = "special")
par(mfrow=c(1,1),oma=c(0.2,1.5,0.2,1.5),mar=c(3,2,0.2,2),cex.axis=1,las=1,mgp=c(1,0.5,0),adj=0.5)
ggarrange(p1 + geom_line(position=pd) +
            geom_point(position=pd) +  geom_errorbar(aes(ymin=Value-ci, ymax=Value+ci), width=.3) +
            rremove("ylab") + xlab('K') + theme_bw()+ ylab('ConsensusLL')  +
            scale_color_manual(values=c('#999999','#000000')),
          p3 + geom_line(position=pd) +
            geom_point(position=pd) +  geom_errorbar(aes(ymin=Value-ci, ymax=Value+ci), width=.3) +
            rremove("ylab") + xlab('K')+ theme_bw()+ylab('RFmetric')  +
            scale_color_manual(values=c('#999999','#000000')),
          
          ncol = 2, nrow = 1)

dev.off()




