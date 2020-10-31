# Code for GMYC analysis of phylogenetic tree ----------------------------------
install.packages("ape")
install.packages("paran")
install.packages("splits", repos = "http://R-Forge.R-project.org")
library(ape)
library(paran)
library(splits)

# Creating ultrametric tree from Bayesian tree constructed in MrBayes
Adineta1 <- read.tree("C:\\R\\Adineta20.tre") #the tree was converted to Newick format using TreeGraph4
Adineta2 <- drop.tip(Adineta1, "FLM") #removing outgroup
Adineta3 <- chronopl(Adineta2, lambda = 100)
Adineta4 <- multi2di(Adineta3) #resolving multichotomies in the ultrametric tree
plot(Adineta4)

# Optimizing COX1 clusters using GMYC, same method is used for ultrametric trees
# converted from Bayesian, or generated in BEAST and converted to Newick format
gmyc.ad.single <- gmyc(Adineta4, method = "single", interval = c(0, 10)) #optimizing single threshhold version of GMYC
summary(gmyc.ad.single)
plot(gmyc.ad.single)
gmyc.ad.multiple <- gmyc(Adineta4, method = "multiple", interval = c(0, 10)) #optimizing multiple threshhold version of GMYC
summary(gmyc.ad.multiple)
plot(gmyc.ad.multiple)
