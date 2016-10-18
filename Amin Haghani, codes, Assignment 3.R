######################################
# Amin Haghani
# Assignment 3
# BISC 481
#############################################################
#(1) Application of an open-source and distributed revision control project: 
#(a) Create a public repository with a README at GitHub https://github.com. 
#(b) Write your name in the file of README.md. 
#(c) Add a collaborator (username: TsuPeiChiu) to the repository. You are required to push your report and R scripts to the repository. The example and file template are shown in https://github.com/TsuPeiChiu/BISC481. 5 pts.

#############################################################

#(2) High-throughput binding assays: Briefly describe 
#(a) the in vitro experiments SELEX-seq and PBM, and 
#(b) the in vivo experiment ChIP-seq. 
#(c) Compare and discuss the advantage and disadvantage of these methods. 10 pts.

#############################################################
#Question 3
#(3) Preparation of high-throughput in vitro data analysis: 
#(a) Download and install R (version >= 3.3.0) from https://www.r-project.org. 
#(b) Install Bioconductor on your R platform. The installation instruction can be found at https://www.bioconductor.org/install/. 
#(c) Install package DNAshapeR on your R platform. The installation instruction can be found at https://www.bioconductor.org/packages/devel/bioc/html/DNAshapeR.html 
#(d) Install the machine learning package caret on your R platform. The installation instruction can be found at https://github.com/topepo/caret 
#(e) Download the gcPBM in vitro experimental data of Mad, Max and Myc from https://github.com/TsuPeiChiu/BISC481/tree/master/gcPBM. 10 pts.

## Install packages
# Bioconductor
source("https://bioconductor.org/biocLite.R")
biocLite()
# DNAshapeR
biocLite("DNAshapeR")
# Caret
install.packages("caret")

## Initialization
library(DNAshapeR)
library(caret)

################################################
#Question 4
#(4) Build prediction models for in vitro data: 
#(a) Use the DNAshapeR package to generate a feature vector for “1-mer” sequence model and a feature vector for “1-mer+shape” model with respect to the datasets of Mad, Max and Myc. 
#(b) Use the caret package to build L2-regularized MLR models for “1-mer” and “1-mer+shape” features with 10-fold cross validation, and print out the average R 2 (coefficient of determination) for these two models with respect to the Mad, Max and Myc datasets. 20 pts.

## Predict DNA shapes
Mad.fa <- "gcPBM/Mad.txt.fa"
Max.fa <- "gcPBM/Max.txt.fa"
Myc.fa <- "gcPBM/Myc.txt.fa"
pred.Mad <- getShape(Mad.fa)
pred.Max <- getShape(Max.fa)
pred.Myc <- getShape(Myc.fa)

## Encode feature vectors
featuretype <- c("1-mer")
featuretype2 <- c("1-mer", "1-shape")
featureVector.Mad <- encodeSeqShape(Mad.fa, pred.Mad, featuretype)
featureVector.Max <- encodeSeqShape(Max.fa, pred.Max, featuretype)
featureVector.Myc <- encodeSeqShape(Myc.fa, pred.Myc, featuretype)

featureVector.Mad2 <- encodeSeqShape(Mad.fa, pred.Mad, featuretype2)
featureVector.Max2 <- encodeSeqShape(Max.fa, pred.Max, featuretype2)
featureVector.Myc2 <- encodeSeqShape(Myc.fa, pred.Myc, featuretype2)


## Build MLR model by using Caret
# Data preparation
fn_exp.mad <- "gcPBM/Mad.txt"
fn_exp.max <- "gcPBM/Max.txt"
fn_exp.myc <- "gcPBM/Myc.txt"
exp_data.mad <- read.table(fn_exp.mad)
exp_data.max <- read.table(fn_exp.max)
exp_data.myc <- read.table(fn_exp.myc)

df1 <- data.frame(affinity=exp_data.mad$V2, featureVector.Mad)
df2 <- data.frame(affinity=exp_data.mad$V2, featureVector.Mad2)

df3 <- data.frame(affinity=exp_data.max$V2, featureVector.Max)
df4 <- data.frame(affinity=exp_data.max$V2, featureVector.Max2)

df5 <- data.frame(affinity=exp_data.myc$V2, featureVector.Myc)
df6 <- data.frame(affinity=exp_data.myc$V2, featureVector.Myc2)

# Arguments setting for Caret
# "cv" mean cross validation
trainControl <- trainControl(method = "cv", number = 10, savePredictions = TRUE)

# Prediction without L2-regularized
model <- train (affinity~., data = df1, trControl=trainControl, 
                method = "lm", preProcess=NULL)
summary(model)

# Prediction with L2-regularized
model.mad.1 <- train(affinity~., data = df1, trControl=trainControl, 
                method = "glmnet", tuneGrid = data.frame(alpha = 0, lambda = c(2^c(-15:15))))
model.mad.2 <- train(affinity~., data = df2, trControl=trainControl, 
                     method = "glmnet", tuneGrid = data.frame(alpha = 0, lambda = c(2^c(-15:15))))
model.max.1 <- train(affinity~., data = df3, trControl=trainControl, 
                     method = "glmnet", tuneGrid = data.frame(alpha = 0, lambda = c(2^c(-15:15))))
model.max.2 <- train(affinity~., data = df4, trControl=trainControl, 
                     method = "glmnet", tuneGrid = data.frame(alpha = 0, lambda = c(2^c(-15:15))))
model.myc.1 <- train(affinity~., data = df5, trControl=trainControl, 
                     method = "glmnet", tuneGrid = data.frame(alpha = 0, lambda = c(2^c(-15:15))))
model.myc.2 <- train(affinity~., data = df6, trControl=trainControl, 
                     method = "glmnet", tuneGrid = data.frame(alpha = 0, lambda = c(2^c(-15:15))))


results <- data.frame(rownames = c("Mad seq", "Mad seq+shape", "Max seq", "Max seq+shape", "Myc seq", "Myc seq+shape")
                      , Rsquared = rbind(model.mad.1$results$Rsquared[1]
                      , model.mad.2$results$Rsquared[1]
                      , model.max.1$results$Rsquared[1]
                      , model.max.2$results$Rsquared[1]
                      , model.myc.1$results$Rsquared[1]
                      , model.myc.2$results$Rsquared[1])
                      , data = rbind("Mad", "Mad", "Max", "Max", "Myc", "Myc")
                      , model = rbind("seq", "seq+shape", "seq", "seq+shape", "seq", "seq+shape")
                      )
write.table(results, file = "results, in vitro.csv", sep = "\t", quote = FALSE)

########################################################
#Question 5
#(5) High-throughput in vitro data analysis: 
#(a) Draw a plot for the comparison of two different models (1mer vs. 1mer+shape) as shown in Figure 1B of Zhou et al. PNAS 2015. 
#(b) Briefly discuss what you have learned from the results. 15 pts.

## Theme
my.theme <- theme(
  plot.margin = unit(c(0.1, 0.5, 0.1, 0.1), "cm"),
  axis.text = element_text(colour="black", size=12),
  axis.title.x = element_text(colour="black", size=12),
  axis.title.y = element_text(colour="black", size=12),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  axis.line = element_line(colour = "black"),
  axis.text = element_text(colour ="black"),
  axis.ticks = element_line(colour = "black")
)

## Data preparation
Seq <- results$Rsquared[results$model=="seq"]
Seq_shape <- results$Rsquared[results$model=="seq+shape"]

## Ploting
ggplot() +
  geom_point(aes(x = Seq, y = Seq_shape), color = "red", size=1) +
  geom_abline(slope=1) + geom_vline(xintercept=0) + geom_hline(yintercept=0) +
  coord_fixed(ratio = 1, xlim = c(0,1), ylim = c(0,1)) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  my.theme  

##plotshape and heatshape of the results (not included in the question)
layout(matrix(1:4, nrow = 2))
plotShape(pred.Mad$MGW, main= "Mad Minor groove width")
plotShape(pred.Mad$HelT, main= "Mad Helix twist")
plotShape(pred.Mad$ProT, main="Mad Propellar twist")
plotShape(pred.Mad$Roll, main= "Mad Roll")

plotShape(pred.Max$MGW, main= "Max Minor groove width")
plotShape(pred.Max$HelT, main= "Max Helix twist")
plotShape(pred.Max$ProT, main="Max Propellar twist")
plotShape(pred.Max$Roll, main= "Max Roll")

plotShape(pred.Myc$MGW, main= "Myc Minor groove width")
plotShape(pred.Myc$HelT, main= "Myc Helix twist")
plotShape(pred.Myc$ProT, main="Myc Propellar twist")
plotShape(pred.Myc$Roll, main= "Myc Roll")

heatShape(pred.Mad$MGW, main= "Mad Minor groove width",4)
heatShape(pred.Mad$HelT, main= "Mad Helix twist",5)
heatShape(pred.Mad$ProT, main="Mad Propellar twist",4)
heatShape(pred.Mad$Roll, main= "Mad Roll",5)

heatShape(pred.Max$MGW, main= "Max Minor groove width",4)
heatShape(pred.Max$HelT, main= "Max Helix twist",5)
heatShape(pred.Max$ProT, main="Max Propellar twist",4)
heatShape(pred.Max$Roll, main= "Max Roll",5)

heatShape(pred.Myc$MGW, main= "Myc Minor groove width",4)
heatShape(pred.Myc$HelT, main= "Myc Helix twist",5)
heatShape(pred.Myc$ProT, main="Myc Propellar twist",4)
heatShape(pred.Myc$Roll, main= "Myc Roll",5)

####################################################
#Question 6 (Logistic regression)
#(6) Preparation of high-throughput in vivo data analysis: 
#(a) Download the ChIP-seq data (including “bound” and “non-bound” data) of CTCF transcription factor of Mus musculus from https://github.com/TsuPeiChiu/BISC481/tree/master/CTCF. 
#(b) Install the R packages mentioned in question (3). 5 pts.

## Install packages
# Bioconductor
source("https://bioconductor.org/biocLite.R")
biocLite()
# DNAshapeR
biocLite("DNAshapeR")
# Caret
install.packages("caret")
install.packages("e1071")
install.packages("ROCR")

## Initialization
library(DNAshapeR)
library(caret)
library(ROCR)

## Prepare data
# Predict DNA shapes
fn_fasta_bound <- "CTCF/bound_500.fa"
fn_fasta_unbound <- "CTCF/unbound_500.fa"

##############################################
#Question 7
#(7) High-throughput in vivo data analysis: 
#(a) Use plotShape() or heatShape() functions of DNAshapeR to generate ensemble plots for the DNA shape parameters of minor groove width (MGW), propeller twist (ProT), Roll, and helix twist (HelT) based on the sequences downloaded for question (6). 
#(b) Briefly discuss what you have learned from the results. 15 pts.

pred1 <- getShape(fn_fasta_bound)
pred2 <- getShape(fn_fasta_unbound)

layout(matrix(1:4, nrow = 2))
plotShape(pred1$MGW, main= "Bound Minor groove width")
plotShape(pred1$HelT, main= "Bound Helix twist")
plotShape(pred1$ProT, main="Bound Propellar twist")
plotShape(pred1$Roll, main= "Bound Roll")

heatShape(pred1$MGW, main= "Bound Minor groove width",500)
heatShape(pred1$HelT, main= "Bound Helix twist",499)
heatShape(pred1$ProT, main="Bound Propellar twist",500)
heatShape(pred1$Roll, main= "Bound Roll",499)

plotShape(pred2$MGW, main= "Unbound Minor groove width")
plotShape(pred2$HelT, main= "Unbound Helix twist")
plotShape(pred2$ProT, main="Unbound Propellar twist")
plotShape(pred2$Roll, main= "Unbound Roll")

heatShape(pred2$MGW, main= "Unbound Minor groove width",500)
heatShape(pred2$HelT, main= "Unbound Helix twist",499)
heatShape(pred2$ProT, main="Unbound Propellar twist",500)
heatShape(pred2$Roll, main= "Unbound Roll",499)

###############################################
#Question 8
#(8) Build prediction models for in vitro data: 
#(a) Build logistic regression models for “1-mer” and “1- mer+shape” features, draw a plot of the ROC curves, and calculate the AUC score for each curve. 
#(b) Briefly discuss what you have learned from the results. 20 pts.

library(Biostrings)
library(caret)
library(DNAshapeR)
library(ROCR)

## Generate data for the classifcation (assign Y to bound and N to non-bound)
# bound
workingPath <- "~/Documents/Rstudio/BISC481 R project/Assignment 3/CTCF/"
boundFasta <- readDNAStringSet(paste0(workingPath, "bound_30.fa"))
sequences <- paste(boundFasta)
boundTxt <- data.frame(seq=sequences, isBound="Y")

# non-bound
nonboundFasta <- readDNAStringSet(paste0(workingPath, "unbound_30.fa"))
sequences <- paste(nonboundFasta)
nonboundTxt <- data.frame(seq=sequences, isBound="N")

# merge two datasets
writeXStringSet( c(boundFasta, nonboundFasta), paste0(workingPath, "ctcf.fa"))
exp_data <- rbind(boundTxt, nonboundTxt)


## DNAshapeR prediction
pred <- getShape("ctcf.fa")


## Encode feature vectors
featuretype <- c("1-mer")
featuretype2 <- c("1-mer", "1-shape")
featureVector1 <- encodeSeqShape(paste0(workingPath, "ctcf.fa"), pred, featuretype)
featureVector2 <- encodeSeqShape(paste0(workingPath, "ctcf.fa"), pred, featuretype2)
df7 <- data.frame(isBound = exp_data$isBound, featureVector1)
df8 <- data.frame(isBound = exp_data$isBound, featureVector2)


## Logistic regression
# Set parameters for Caret
trainControl <- trainControl(method = "cv", number = 10, 
                             savePredictions = TRUE, classProbs = TRUE)
# Perform prediction
model1 <- train(isBound~ ., data = df7, trControl = trainControl,
               method = "glm", family = binomial, metric ="ROC")
model2 <- train(isBound~ ., data = df8, trControl = trainControl,
                method = "glm", family = binomial, metric ="ROC")
summary(model1)
summary(model2)

## Plot AUROC
prediction <- prediction( model1$pred$Y, model1$pred$obs )
performance <- performance( prediction, "tpr", "fpr" )
prediction2 <- prediction( model2$pred$Y, model2$pred$obs )
performance2 <- performance( prediction2, "tpr", "fpr" )
layout(cbind(c(1), c(2)))
plot(performance, main= "Sequence model")
text(x = 0.5, y = 0.5,labels = "auc=0.8406051" )
plot(performance2, main= "Sequence+Shape model")
text(x = 0.5, y = 0.5,labels = "auc=0.8408944" )
## Caluculate AUROC
# any Auc more than 0.5 shows a successful model
auc <- performance(prediction, "auc")
auc <- unlist(slot(auc, "y.values"))
auc2 <- performance(prediction2, "auc")
auc2 <- unlist(slot(auc2, "y.values"))
auc
auc2
auc_summary <- data.frame(Models= c("Sequence", "Sequence+Shape")
                          , auc=c(auc, auc2))







