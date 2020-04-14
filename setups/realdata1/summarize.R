library(ggplot2)
library(Rmisc) 
library(Sleuth2)
library(ggpubr)


PGS1000 <- data.frame(read.csv("PGs/prime_SMC1/results/results.csv"))

PGS2000 <- data.frame(read.csv("PGs/prime_SMC2/results/results.csv"))

PGS5000 <- data.frame(read.csv("PGs/prime_SMC3/results2/results.csv"))

for(i in 4:10){
  dir <- paste0("PGs/prime_SMC3/results", i, "/results.csv")
  file <- data.frame(read.csv(dir))
  PGS5000 <- rbind(PGS5000, file)
}

PGS10000 <- data.frame(read.csv("PGs/prime_SMC4/results1/results.csv"))

for(i in 2:10){
  dir <- paste0("PGs/prime_SMC4/results", i, "/results.csv")
  file <- data.frame(read.csv(dir))
  PGS10000 <- rbind(PGS10000, file)
}

PGS1000$K = 1000
PGS2000$K = 2000
PGS5000$K = 5000
PGS10000$K = 10000

results <- rbind(PGS1000, PGS2000, PGS5000, PGS10000)



results_ConsensusLL <- results[which(results$Metric == 'ConsensusLogLL'),]
CLLsummary_PGS <- summarySE(results_ConsensusLL, measurevar="Value", groupvars=c("K"))


##IPMCMC
IPGS1000 <- data.frame(read.csv("IPMCMC/prime1_IPGS/results1/results.csv"))
for(i in 2:10){
  dir <- paste0("IPMCMC/prime1_IPGS/results", i, "/results.csv")
  file <- data.frame(read.csv(dir))
  IPGS1000 <- rbind(IPGS1000, file)
}

IPGS2000 <- data.frame(read.csv("IPMCMC/prime2_IPGS/results1/results.csv"))
for(i in 2:10){
  dir <- paste0("IPMCMC/prime2_IPGS/results", i, "/results.csv")
  file <- data.frame(read.csv(dir))
  IPGS2000 <- rbind(IPGS2000, file)
}

IPGS5000 <- data.frame(read.csv("IPMCMC/prime4_IPGS/results1/results.csv"))
for(i in 2:10){
  dir <- paste0("IPMCMC/prime4_IPGS/results", i, "/results.csv")
  file <- data.frame(read.csv(dir))
  IPGS5000 <- rbind(IPGS5000, file)
}

IPGS10000 <- data.frame(read.csv("IPMCMC/prime3_IPGS/results1/results.csv"))
for(i in 2:10){
  dir <- paste0("IPMCMC/prime3_IPGS/results", i, "/results.csv")
  file <- data.frame(read.csv(dir))
  IPGS10000 <- rbind(IPGS10000, file)
}

IPGS1000$K = 1000
IPGS2000$K = 2000
IPGS5000$K = 5000
IPGS10000$K = 10000

results <- rbind(IPGS1000, IPGS2000, IPGS5000, IPGS10000)



results_ConsensusLL <- results[which(results$Metric == 'ConsensusLogLL'),]
CLLsummary_IPGS <- summarySE(results_ConsensusLL, measurevar="Value", groupvars=c("K"))

##theta ###
IPGS_theta <- data.frame(read.csv("IPMCMC/prime4_IPGS/results1/PGS.csv"))$trans2tranv
IPGS_theta_burn <- IPGS_theta[-(1:2000)]
mean_theta <- mean(IPGS_theta_burn)
for(i in 2:10){
  dir <- paste0("IPMCMC/prime4_IPGS/results", i, "/PGS.csv")
  file <- data.frame(read.csv(dir))$trans2tranv
  IPGS_theta_burn <- file[-(1:2000)]
  mean_theta <- c(mean_theta, mean(IPGS_theta_burn))
  
  #IPGS5000 <- rbind(IPGS5000, file)
}

library(ggplot2)
library(Rmisc) 
library(Sleuth2)
library(ggpubr)
library(stringr)

mean_p1 <- mean_theta

data_p11 <- data.frame(p = as.vector(IPGS_theta), Iteration = as.vector(1:length(IPGS_theta)) )

g1 <- ggplot(data_p11, aes(x=Iteration, y=p)) + geom_line() + theme_bw()+ xlab("Iteration") + ylab(expression(theta)) + scale_color_manual(values = c('#595959', 'blue'))+ geom_hline(yintercept=mean_p1[1], linetype="dashed", color = "yellow")

gname = c("traceplot.eps",sep="")  
postscript(gname,width=5,height=3,horizontal = FALSE, onefile = FALSE, paper = "special")
par(mfrow=c(1,1),oma=c(0.2,1.5,0.2,1.5),mar=c(3,2,0.2,2),cex.axis=1,las=1,mgp=c(1,0.5,0),adj=0.5)
ggarrange(g1,
          ncol = 1, nrow = 1, common.legend = TRUE)
dev.off()



