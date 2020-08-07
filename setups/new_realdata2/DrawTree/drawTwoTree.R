
library(treespace)
library(phytools)
# PGS
tr1=read.tree(text="(((Callochromis_macrops:0.0321,(Xenotilapia_sima:0.0047,Cardiopharynx_schoutedeni:0.0047)internal0:0.0274)internal1:0.0085,((Chalinochromis_popeleni:0.00575,Telmatochromis_temporalis:0.00575)internal2:0.00855,(Lepidiolamprologus_elongatus:0.01085,Julidochromis_marlieri:0.01085)internal3:0.00345)internal4:0.0263)internal5:0.0559,(((Lamprologus_callipterus:0.00515,Neolamprologus_brichardi:0.00515)internal9:0.02845,Ophthalmotilapia_ventralis:0.0336)internal8:0.0092,Neolamprologus_tetracanthus:0.0428)internal7:0.0556,Xenotilapia_flavipinnus:0.0965)internal6;")
# PGAS
tr2=read.tree(text="(((Lepidiolamprologus_elongatus:0.04195,Julidochromis_marlieri:0.04195)internal4:0.0179,(Xenotilapia_flavipinnus:0.03165,(Xenotilapia_sima:0.01885,(Callochromis_macrops:0.01345,Cardiopharynx_schoutedeni:0.01345)internal1:0.0054)internal2:0.0128)internal3:0.0282)internal5:0.0021,(((Chalinochromis_popeleni:0.02215,Lamprologus_callipterus:0.02215)internal9:0.0416,Telmatochromis_temporalis:0.06375)internal8:0.0195,Neolamprologus_tetracanthus:0.07895)internal7:0.0162,(Ophthalmotilapia_ventralis:0.00575,Neolamprologus_brichardi:0.00575)internal0:0.0562)internal6;")

tr1 <- reroot(tr1, node.number = 2)
tr2 <- reroot(tr2, node.number = 1)



tipDiff <- tipDiff(tr1,tr2)
#plotTreeDiff(tr1,tr2,tipDiff=tipDiff, baseCol="black", colourMethod="palette", 
#             edge.width=2, type="cladogram", cex=0.5, font=2)
postscript("consensusChildFish.eps",width=5,height=4,horizontal = FALSE, onefile = FALSE, paper = "special")
plotTreeDiff(tr1,tr2,tipDiff=tipDiff, baseCol="black", colourMethod="palette", 
                 edge.width=2,  cex=0.5, font=2)
dev.off()             
