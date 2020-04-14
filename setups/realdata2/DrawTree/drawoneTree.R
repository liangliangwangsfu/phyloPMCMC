library(ape)
library(treespace)
library (rpart)
library(phytools)

text.string<-
  "((ce-macaque:0.1748,(gibbon:0.1213,(orang-utan:0.09615,(gorilla:0.056,(chimpanzee:0.04225,human:0.04225)internal1:0.01375)internal2:0.04015)internal3:0.02515)internal4:0.0535)internal5:0.0393,s-monkey:0.2139,(lemur:0.16355,tarsier:0.16355)internal0:0.09335)internal6;"
vert.tree<-read.tree(text=text.string)

#vert.tree$root.edge <- 0

postscript("primate.eps",width=5,height=4,horizontal = FALSE, onefile = FALSE, paper = "special")
plot(vert.tree,no.margin=TRUE,edge.width=2)
dev.off() 




plot(root(vert.tree, "tarsier"))
