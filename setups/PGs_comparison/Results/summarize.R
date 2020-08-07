library(ggplot2)
library(Rmisc) 
library(Sleuth2)
library(ggpubr)


PGS200_1 <- data.frame(read.csv("PGs/K200/results1.csv"))
PGS200_2 <- data.frame(read.csv("PGs/K200/results2.csv"))
PGS200_3 <- data.frame(read.csv("PGs/K200/results3.csv"))
PGS200_4 <- data.frame(read.csv("PGs/K200/results4.csv"))
PGS200_5 <- data.frame(read.csv("PGs/K200/results5.csv"))
PGS200 <- rbind(PGS200_1, PGS200_2, PGS200_3, PGS200_4, PGS200_5)
PGS200$K = 200

PGS500_1 <- data.frame(read.csv("PGs/K500/results1.csv"))
PGS500_2 <- data.frame(read.csv("PGs/K500/results2.csv"))
PGS500_3 <- data.frame(read.csv("PGs/K500/results3.csv"))
PGS500_4 <- data.frame(read.csv("PGs/K500/results4.csv"))
PGS500_5 <- data.frame(read.csv("PGs/K500/results5.csv"))
PGS500 <- rbind(PGS500_1, PGS500_2, PGS500_3, PGS500_4, PGS500_5)
PGS500$K = 500

PGS1000_1 <- data.frame(read.csv("PGs/K1000/results1.csv"))
PGS1000_2 <- data.frame(read.csv("PGs/K1000/results2.csv"))
PGS1000_3 <- data.frame(read.csv("PGs/K1000/results3.csv"))
PGS1000_4 <- data.frame(read.csv("PGs/K1000/results4.csv"))
PGS1000_5 <- data.frame(read.csv("PGs/K1000/results5.csv"))
PGS1000 <- rbind(PGS1000_1, PGS1000_2, PGS1000_3, PGS1000_4, PGS1000_5)
PGS1000$K = 1000

PGS <- rbind(PGS200, PGS500, PGS1000)
PGS$method <- 'PGs'


PGSRM200_1 <- data.frame(read.csv("PGsRM/K200/results1.csv"))
PGSRM200_2 <- data.frame(read.csv("PGsRM/K200/results2.csv"))
PGSRM200_3 <- data.frame(read.csv("PGsRM/K200/results3.csv"))
PGSRM200_4 <- data.frame(read.csv("PGsRM/K200/results4.csv"))
PGSRM200_5 <- data.frame(read.csv("PGsRM/K200/results5.csv"))
PGSRM200 <- rbind(PGSRM200_1, PGSRM200_2, PGSRM200_3, PGSRM200_4, PGSRM200_5)
PGSRM200$K = 200

PGSRM500_1 <- data.frame(read.csv("PGsRM/K500/results1.csv"))
PGSRM500_2 <- data.frame(read.csv("PGsRM/K500/results2.csv"))
PGSRM500_3 <- data.frame(read.csv("PGsRM/K500/results3.csv"))
PGSRM500_4 <- data.frame(read.csv("PGsRM/K500/results4.csv"))
PGSRM500_5 <- data.frame(read.csv("PGsRM/K500/results5.csv"))
PGSRM500 <- rbind(PGSRM500_1, PGSRM500_2, PGSRM500_3, PGSRM500_4, PGSRM500_5)
PGSRM500$K = 500

PGSRM1000_1 <- data.frame(read.csv("PGsRM/K1000/results1.csv"))
PGSRM1000_2 <- data.frame(read.csv("PGsRM/K1000/results2.csv"))
PGSRM1000_3 <- data.frame(read.csv("PGsRM/K1000/results3.csv"))
PGSRM1000_4 <- data.frame(read.csv("PGsRM/K1000/results4.csv"))
PGSRM1000_5 <- data.frame(read.csv("PGsRM/K1000/results5.csv"))
PGSRM1000 <- rbind(PGSRM1000_1, PGSRM1000_2, PGSRM1000_3, PGSRM1000_4, PGSRM1000_5)
PGSRM1000$K = 1000

PGSRM <- rbind(PGSRM200, PGSRM500, PGSRM1000)
PGSRM$method <- 'PGs-RDouP'
#PGS4 <- data.frame(read.csv("results4/results.csv"))

IPGS200_1 <- data.frame(read.csv("IPGs/K200/results1.csv"))
IPGS200_2 <- data.frame(read.csv("IPGs/K200/results2.csv"))
IPGS200_3 <- data.frame(read.csv("IPGs/K200/results3.csv"))
IPGS200_4 <- data.frame(read.csv("IPGs/K200/results4.csv"))
IPGS200_5 <- data.frame(read.csv("IPGs/K200/results5.csv"))
IPGS200 <- rbind(IPGS200_1, IPGS200_2, IPGS200_3, IPGS200_4, IPGS200_5)
IPGS200$K = 200

IPGS500_1 <- data.frame(read.csv("IPGs/K500/results1.csv"))
IPGS500_2 <- data.frame(read.csv("IPGs/K500/results2.csv"))
IPGS500_3 <- data.frame(read.csv("IPGs/K500/results3.csv"))
IPGS500_4 <- data.frame(read.csv("IPGs/K500/results4.csv"))
IPGS500_5 <- data.frame(read.csv("IPGs/K500/results5.csv"))
IPGS500 <- rbind(IPGS500_1, IPGS500_2, IPGS500_3, IPGS500_4, IPGS500_5)
IPGS500$K = 500

IPGS1000_1 <- data.frame(read.csv("IPGs/K1000/results1.csv"))
IPGS1000_2 <- data.frame(read.csv("IPGs/K1000/results2.csv"))
IPGS1000_3 <- data.frame(read.csv("IPGs/K1000/results3.csv"))
IPGS1000_4 <- data.frame(read.csv("IPGs/K1000/results4.csv"))
IPGS1000_5 <- data.frame(read.csv("IPGs/K1000/results5.csv"))
IPGS1000 <- rbind(IPGS1000_1, IPGS1000_2, IPGS1000_3, IPGS1000_4, IPGS1000_5)
IPGS1000$K = 1000

IPGS <- rbind(IPGS200, IPGS500, IPGS1000)
IPGS$method <- 'IPGs'

IPGSRM200_1 <- data.frame(read.csv("IPGsRM/K200/results1.csv"))
IPGSRM200_2 <- data.frame(read.csv("IPGsRM/K200/results2.csv"))
IPGSRM200_3 <- data.frame(read.csv("IPGsRM/K200/results3.csv"))
IPGSRM200_4 <- data.frame(read.csv("IPGsRM/K200/results4.csv"))
IPGSRM200_5 <- data.frame(read.csv("IPGsRM/K200/results5.csv"))
IPGSRM200 <- rbind(IPGSRM200_1, IPGSRM200_2, IPGSRM200_3, IPGSRM200_4, IPGSRM200_5)
IPGSRM200$K = 200

IPGSRM500_1 <- data.frame(read.csv("IPGsRM/K500/results1.csv"))
IPGSRM500_2 <- data.frame(read.csv("IPGsRM/K500/results2.csv"))
IPGSRM500_3 <- data.frame(read.csv("IPGsRM/K500/results3.csv"))
IPGSRM500_4 <- data.frame(read.csv("IPGsRM/K500/results4.csv"))
IPGSRM500_5 <- data.frame(read.csv("IPGsRM/K500/results5.csv"))
IPGSRM500 <- rbind(IPGSRM500_1, IPGSRM500_2, IPGSRM500_3, IPGSRM500_4, IPGSRM500_5)
IPGSRM500$K = 500

IPGSRM1000_1 <- data.frame(read.csv("IPGsRM/K1000/results1.csv"))
IPGSRM1000_2 <- data.frame(read.csv("IPGsRM/K1000/results2.csv"))
IPGSRM1000_3 <- data.frame(read.csv("IPGsRM/K1000/results3.csv"))
IPGSRM1000_4 <- data.frame(read.csv("IPGsRM/K1000/results4.csv"))
IPGSRM1000_5 <- data.frame(read.csv("IPGsRM/K1000/results5.csv"))
IPGSRM1000 <- rbind(IPGSRM1000_1, IPGSRM1000_2, IPGSRM1000_3, IPGSRM1000_4, IPGSRM1000_5)
IPGSRM1000$K = 1000

IPGSRM <- rbind(IPGSRM200, IPGSRM500, IPGSRM1000)
IPGSRM$method <- 'IPGs-RDouP'

PGAS200_1 <- data.frame(read.csv("PGAs/K200/results1.csv"))
PGAS200_2 <- data.frame(read.csv("PGAs/K200/results2.csv"))
PGAS200_3 <- data.frame(read.csv("PGAs/K200/results3.csv"))
PGAS200_4 <- data.frame(read.csv("PGAs/K200/results4.csv"))
PGAS200_5 <- data.frame(read.csv("PGAs/K200/results5.csv"))
PGAS200 <- rbind(PGAS200_1, PGAS200_2, PGAS200_3, PGAS200_4, PGAS200_5)
PGAS200$K = 200

PGAS500_1 <- data.frame(read.csv("PGAs/K500/results1.csv"))
PGAS500_2 <- data.frame(read.csv("PGAs/K500/results2.csv"))
PGAS500_3 <- data.frame(read.csv("PGAs/K500/results3.csv"))
PGAS500_4 <- data.frame(read.csv("PGAs/K500/results4.csv"))
PGAS500_5 <- data.frame(read.csv("PGAs/K500/results5.csv"))
PGAS500 <- rbind(PGAS500_1, PGAS500_2, PGAS500_3, PGAS500_4, PGAS500_5)
PGAS500$K = 500

PGAS1000_1 <- data.frame(read.csv("PGAs/K1000/results1.csv"))
PGAS1000_2 <- data.frame(read.csv("PGAs/K1000/results2.csv"))
PGAS1000_3 <- data.frame(read.csv("PGAs/K1000/results3.csv"))
PGAS1000_4 <- data.frame(read.csv("PGAs/K1000/results4.csv"))
PGAS1000_5 <- data.frame(read.csv("PGAs/K1000/results5.csv"))
PGAS1000 <- rbind(PGAS1000_1, PGAS1000_2, PGAS1000_3, PGAS1000_4, PGAS1000_5)
PGAS1000$K = 1000

PGAS <- rbind(PGAS200, PGAS500, PGAS1000)
PGAS$method <- 'PGAs'

#results <- rbind(PGS1, PGS2, PGS3, PGS4)
results <- rbind(PGS, PGSRM, IPGS, IPGSRM, PGAS)
results_RF <- results[which(results$Metric == 'RobinsonFouldsMetric'),]
#results_RF$K <-results_RF$IterScale

RFsummary <- summarySE(results_RF, measurevar="Value", groupvars=c("K", "method"))
RFsummary$K = as.factor(RFsummary$K)

results_ConsensusLL <- results[which(results$Metric == 'ConsensusLogLL'),]


CLLsummary <- summarySE(results_ConsensusLL, measurevar="Value", groupvars=c("K", "method"))
CLLsummary$K = as.factor(CLLsummary$K)

p2 <- ggplot(RFsummary, aes(x = K, y=Value, group = method, colour=method))
#p2 <- ggplot(PMsummary_PGSRM, aes(x = K, y=Value, group = 1))
p1 <- ggplot(CLLsummary, aes(x = K, y=Value, group = method, colour=method))



pd <- position_dodge(.2)
gname = c("PGS.eps",sep="")  
postscript(gname,width=5,height=2,horizontal = FALSE, onefile = FALSE, paper = "special")
par(mfrow=c(1,1),oma=c(0.2,1.5,0.2,1.5),mar=c(3,2,0.2,2),cex.axis=1,las=1,mgp=c(1,0.5,0),adj=0.5)
ggarrange(  p1 + geom_line(position=pd) +
              geom_point(position=pd) +  geom_errorbar(aes(ymin=Value-ci, ymax=Value+ci), width=.3, position=pd) +
              rremove("ylab") + xlab('K')+ theme_bw() + ylab('ConsensusLL')  +
              scale_color_manual(values=c('red','#000000', "#333BFF", "#E69F00", '#A4A4A4')),
            p2 + geom_line(position=pd) +
              geom_point(position=pd) +  geom_errorbar(aes(ymin=Value-ci, ymax=Value+ci), width=.3, position=pd) +
              rremove("ylab") + xlab('K')+ theme_bw()+ylab('RFmetric')  +
              scale_color_manual(values=c('red','#000000', "#333BFF", "#E69F00", '#A4A4A4')),
            ncol = 2, nrow = 1, common.legend = TRUE)

dev.off()


###only compare PGs-RM and IPGS-RM
results2 <- rbind(IPGSRM, PGAS, PGSRM)
results2_RF <- results2[which(results2$Metric == 'RobinsonFouldsMetric'),]
#results_RF$K <-results_RF$IterScale

RF2summary <- summarySE(results2_RF, measurevar="Value", groupvars=c("K", "method"))
RF2summary$K = as.factor(RF2summary$K)

results2_ConsensusLL <- results2[which(results2$Metric == 'ConsensusLogLL'),]


CLL2summary <- summarySE(results2_ConsensusLL, measurevar="Value", groupvars=c("K", "method"))
CLL2summary$K = as.factor(CLL2summary$K)

p4 <- ggplot(RF2summary, aes(x = K, y=Value, group = method, colour=method))
#p2 <- ggplot(PMsummary_PGSRM, aes(x = K, y=Value, group = 1))
p3 <- ggplot(CLL2summary, aes(x = K, y=Value, group = method, colour=method))


pd <- position_dodge(.2)
gname = c("PGSRM.eps",sep="")
postscript(gname,width=5,height=2,horizontal = FALSE, onefile = FALSE, paper = "special")
par(mfrow=c(1,1),oma=c(0.2,1.5,0.2,1.5),mar=c(3,2,0.2,2),cex.axis=1,las=1,mgp=c(1,0.5,0),adj=0.5)
ggarrange(p3 + geom_line(position=pd) +
            geom_point(position=pd) +  geom_errorbar(aes(ymin=Value-ci, ymax=Value+ci), width=.3) +
            rremove("ylab") + xlab('K') + theme_bw()+ ylab('ConsensusLL')  +
            scale_color_manual(values=c('#000000', "#333BFF", '#A4A4A4')),
          p4 + geom_line(position=pd) +
            geom_point(position=pd) +  geom_errorbar(aes(ymin=Value-ci, ymax=Value+ci), width=.3) +
            rremove("ylab") + xlab('K')+ theme_bw()+ylab('RFmetric')  +
            scale_color_manual(values=c('#000000', "#333BFF", '#A4A4A4')),
          ncol = 2, nrow = 1, common.legend = TRUE)

dev.off()


