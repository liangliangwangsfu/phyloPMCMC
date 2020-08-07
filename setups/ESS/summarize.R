library(ggplot2)
library(Rmisc) 
library(Sleuth2)
library(ggpubr)

CSMC <- data.frame(ESS = c(6.13849598163571,
                           2.243468525711404,
                           22.04765302478196,
                           37.876167828609375,
                           51.26777866766022,
                           36.573951511780216,
                           54.98937108620675,
                           10.508014355802366,
                           71.68418612560117), rank = 1:9)
CSMC_RDouP <- data.frame(ESS = c(57.75011835613956,
                                 1.0000429461439693,
                                 1.0000000046398205,
                                 1.0183892264064913,
                                 2.7662116017628424,
                                 7.490147311729196,
                                 4.110124691135874,
                                 1.506029738085357,
                                 19.604299046940877), rank = 1:9)
CSMC$method <- 'CSMC'
CSMC_RDouP$method <- 'CSMC-RDouP'

DESS <- rbind(CSMC, CSMC_RDouP)
p0 <- ggplot(DESS, aes(rank, ESS, group = method)) + geom_line(aes(linetype=method)) + theme_classic()




gname = c("ESS.eps",sep="")  
postscript(gname,width=4,height=3,horizontal = FALSE, onefile = FALSE, paper = "special")
par(mfrow=c(1,1),oma=c(0.2,1.5,0.2,1.5),mar=c(3,2,0.2,2),cex.axis=1,las=1,mgp=c(1,0.5,0),adj=0.5)

ggarrange(p0,
          ncol = 1, nrow = 1, common.legend = TRUE)

dev.off()






