## Package installations are included here for convenience

install.packages("MPsychoR")
install.packages("qgraph")
install.packages("smacof")
install.packages("wordcloud")
install.packages("psych")
install.packages("eigenmodel")
install.packages("networktools")

## Note: The following R code is identical to code found in the manuscript

library("MPsychoR")
data(Rogers)
dim(Rogers)

data(Rogers_Adolescent)
dim(Rogers_Adolescent)

colnames(Rogers) <- colnames(Rogers_Adolescent) <- 1:26

library("qgraph")
adult_zeroorder <- cor(Rogers)
qgraph(adult_zeroorder, layout="spring",
       groups = list(Depression = 1:16, "OCD" = 17:26), 
       color = c("lightblue", "lightsalmon"))

adult_zeroorder <- cor(Rogers)

library("smacof")
dissimilarity_adult <- sim2diss(adult_zeroorder)

adult_MDS <- mds(dissimilarity_adult)
head(round(adult_MDS$conf, 2)) # top of configuration matrix

adult_MDS_ordinal <- mds(dissimilarity_adult, type="ordinal")
plot(adult_MDS_ordinal, plot.type = "Shepard", main="Ordinal")
text(1.1,0.3, paste("Stress =", round(adult_MDS_ordinal$stress,2))) 

adult_MDS_ratio <- mds(dissimilarity_adult, type="ratio")
plot(adult_MDS_ratio, plot.type = "Shepard", main="Ratio")
text(1.1,0.3, paste("Stress =", round(adult_MDS_ratio$stress,2))) 

adult_MDS_interval <- mds(dissimilarity_adult, type="interval")
plot(adult_MDS_interval, plot.type = "Shepard", main="Interval")
text(1.1,0.3, paste("Stress =", round(adult_MDS_interval$stress,2))) 

adult_MDS_mspline <- mds(dissimilarity_adult, type="mspline")
plot(adult_MDS_mspline, plot.type = "Shepard", main="Spline")
text(1.1,0.3, paste("Stress =", round(adult_MDS_mspline$stress,2)))

adult_MDS_mspline$stress

qgraph(adult_zeroorder, layout=adult_MDS_mspline$conf, 
       groups = list(Depression = 1:16, "OCD" = 17:26), 
       color = c("lightblue", "lightsalmon"), vsize=4)
text(-1,-1, paste("Stress=", round(adult_MDS_mspline$stress,2)))

library("wordcloud")
qgraph(adult_zeroorder, layout=adult_MDS_mspline$conf, 
       groups = list(Depression = 1:16, "OCD" = 17:26), 
       color = c("lightblue", "lightsalmon"),
       vsize=0, rescale=FALSE, labels=FALSE)
points(adult_MDS_mspline$conf, pch=16)
textplot(adult_MDS_mspline$conf[,1]+.03,
         adult_MDS_mspline$conf[,2]+.03,
         colnames(adult_zeroorder),
         new=F)

adult_glasso <- EBICglasso(cor(Rogers), n=408)
qgraph(adult_glasso, layout=adult_MDS_mspline$conf, 
       groups = list(Depression = 1:16, "OCD" = 17:26), 
       color = c("lightblue", "lightsalmon"), vsize=4)
text(-1,-1, paste("Stress=", round(adult_MDS_mspline$stress,2))) 

adolescent_zeroorder <- cor(Rogers_Adolescent)
dissimilarity_adolescent <- sim2diss(adolescent_zeroorder)
adolescent_MDS <- mds(dissimilarity_adolescent, type="mspline")

fit_procrustes <- Procrustes(adult_MDS_mspline$conf, adolescent_MDS$conf)

adolescent_glasso <- EBICglasso(cor(Rogers_Adolescent), n=87, gamma=0)

qgraph(adult_glasso, layout=fit_procrustes$X, groups = list(Depression = 1:16, "OCD" = 17:26),
       color = c("lightblue", "lightsalmon"), title= "Adults, n=408", vsize=4)
text(-1,-1, paste("Stress=", round(adult_MDS_mspline$stress,2)))
qgraph(adolescent_glasso, layout=fit_procrustes$Yhat, 
       groups = list(Depression = 1:16, "OCD" = 17:26),
       color = c("lightblue", "lightsalmon"), title="Adolescents, n=87", vsize=4)
text(-1,-1, paste("Stress=", round(adolescent_MDS$stress,2)))

round(fit_procrustes$congcoef, 3)

library("psych")
PCA_adult <- principal(cor(Rogers), nfactors = 2)
qgraph(adult_glasso, layout=PCA_adult$loadings, groups = list(Depression = 1:16, "OCD" = 17:26), 
       color = c("lightblue", "lightsalmon"), title= "Adults, n=408", layoutOffset=c(.3,.1), vsize=4)

text(1.5,-.8, paste("% var=", round(sum(PCA_adult$values[1:2]/length(PCA_adult$values)),2)))
title(xlab="Component 1", ylab= "Component 2")

library("eigenmodel")

diag(adult_glasso) <- NA   ## the function needs NA diagonals
p <- 2               		## 2-dimensional solution
fitEM <- eigenmodel_mcmc(Y = adult_glasso, R = p, S = 1000, burn = 200, seed = 123)
EVD <- eigen(fitEM$ULU_postmean) 
evecs <- EVD$vec[, 1:p]      ## eigenvectors (coordinates)

qgraph(adult_glasso, layout=evecs, groups = list(Depression = 1:16, "OCD" = 17:26), 
       color = c("lightblue", "lightsalmon"), title= "Adults, n=408", vsize=4)
title(xlab="Dimension 1", ylab= "Dimension 2")

library("networktools")
adult_glasso <- EBICglasso(cor(Rogers), n=408)
adult_qgraph <- qgraph(adult_glasso)
MDSnet(adult_qgraph, MDSadj=cor(Rogers))
PCAnet(adult_qgraph, cormat = cor(Rogers))
EIGENnet(adult_qgraph)

