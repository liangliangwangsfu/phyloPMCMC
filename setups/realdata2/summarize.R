library(ggplot2)
library(Rmisc) 
library(Sleuth2)
library(ggpubr)

##theta ###
IPGS_theta <- data.frame(read.csv("Fish_PGAS/results1/PGAS4K2PBF.csv"))$trans2tranv
IPGS_theta_burn <- IPGS_theta[-(1:2000)]
mean_theta <- mean(IPGS_theta_burn)

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



