
library(treespace)
library(phytools)
# from SMC sampler (deterministic after the adaptive scheme)
tr1=read.tree(text="(((Callochromis_macrops:0.0216,(Chalinochromis_popeleni:0.0012,Xenotilapia_sima:0.0012)internal1:0.0204)internal2:0.02105,(Ophthalmotilapia_ventralis:0.0116,(Xenotilapia_flavipinnus:0.00165,Cardiopharynx_schoutedeni:0.00165)internal3:0.00995)internal4:0.03105)internal5:0.00885,(((Lepidiolamprologus_elongatus:0.00975,Lamprologus_callipterus:0.00975)internal9:0.00425,(Julidochromis_marlieri:0.0012,Neolamprologus_brichardi:0.0012)internal0:0.0128)internal8:0.04455,Telmatochromis_temporalis:0.05685)internal7:0.00535,Neolamprologus_tetracanthus:0.0515)internal6;")
# from MCMC running for a really long time. 
tr2=read.tree(text="(((Ophthalmotilapia_ventralis:0.0026,Neolamprologus_brichardi:0.0026)internal1:0.0538,((Callochromis_macrops:0.00715,Cardiopharynx_schoutedeni:0.00715)internal2:0.0137,(Xenotilapia_flavipinnus:0.0196,Xenotilapia_sima:0.0196)internal3:0.00125)internal4:0.03555)internal5:0.0051,(((Chalinochromis_popeleni:0.0366,Lamprologus_callipterus:0.0366)internal9:0.0279,Telmatochromis_temporalis:0.0645)internal8:0.01835,Neolamprologus_tetracanthus:0.07745)internal7:0.01625,(Lepidiolamprologus_elongatus:0.03895,Julidochromis_marlieri:0.03895)internal0:0.02255)internal6;")

tr1 <- reroot(tr1, node.number = 2)
tr2 <- reroot(tr2, node.number = 1)



tipDiff <- tipDiff(tr1,tr2)
#plotTreeDiff(tr1,tr2,tipDiff=tipDiff, baseCol="black", colourMethod="palette", 
#             edge.width=2, type="cladogram", cex=0.5, font=2)
postscript("consensusChildFish.eps",width=5,height=4,horizontal = FALSE, onefile = FALSE, paper = "special")
plotTreeDiff(tr1,tr2,tipDiff=tipDiff, baseCol="black", colourMethod="palette", 
                 edge.width=2,  cex=0.5, font=2)
dev.off()             
