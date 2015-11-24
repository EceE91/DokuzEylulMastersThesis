# Convert txt files to csv files
#tab.gbm.survival = read.delim("GBM-survival-class.txt")
#write.table(tab.gbm.survival, file="GBM-survival-class.csv",sep=",",col.names=TRUE,row.names=FALSE)

#tab = read.delim("GBM-RPPA-CommonPatients-Normalized-v2.txt")
#write.table(tab, file="GBM-RPPA-CommonPatients-Normalized-v2.csv",sep=",",col.names=TRUE,row.names=FALSE)


balanced.sample <- function(pData, n) {
   
    all.good = which(pData$z == "C1")
    all.poor = which(pData$z == "C0")
      
    sample.good = sample(all.good, round(length(all.good)*n), replace=FALSE)
    sample.poor = sample(all.poor, round(length(all.poor)*n), replace=FALSE)

    training.set = c(pData[sample.good, 1], pData[sample.poor, 1])  # return patient ids -> TCGA.06.2567
    return(training.set)
}

TCGA.start = Gbm.RPPA[5:39]

# read GBM RPPA
Gbm.RPPA = read.csv("GBM-RPPA-CommonPatients-Normalized-v2.csv",header=T)
Gbm.RPPA.patient = colnames(TCGA.start)

# get only TCGA columns
Gbm.RPPA.TCGA.DATA <- TCGA.start
# get entrezids
#Gbm.RPPA.TCGA.EntrezIds <- as.vector(Gbm.RPPA[3])

# read GBM Survival 
Gbm.Survival = read.csv("GBM-survival-class.csv",header=T)
# replace "-" with "." to match RPPA TCGA data (TCGA-12-1098 => TCGA.12.1098)
Gbm.Survival$feature = gsub("-",".",as.factor(Gbm.Survival$feature))
#rownames(Gbm.Survival) = Gbm.Survival$feature
#Gbm.Survival.TCGA = rownames(Gbm.Survival)
Gbm.Survival.TCGA= Gbm.Survival$feature


# find same TCGA's = patients
Gbm.Patients = intersect(c(Gbm.Survival.TCGA), c(Gbm.RPPA.patient))

# create data frame for new TCGA-SURVIVAL STATUS data
Gbm.Survival.TCGA.Data <- data.frame(x = numeric(), z = character(), stringsAsFactors = FALSE)

# add data to data frame
for(i in 1:length(Gbm.Patients))
{
  pID = Gbm.Patients[i]
  rowNumber = which(Gbm.Survival$feature == pID) # get rownumber of matching data
  Gbm.Survival.TCGA.Data[i,1] = pID
  Gbm.Survival.TCGA.Data[i,2] = paste("C", Gbm.Survival$OS_vital_status[rowNumber], sep="")
}



# change colnames
#colnames(Gbm.Survival.TCGA.Data) <- c("Feature", "Vital_Status_Class_Label")
# write data frame to csv file
#write.table(Gbm.Survival.TCGA.Data, "C:/Users/Ece/Desktop/GBM/GBM_TCGA_SURVIVAL.csv", sep="\t", row.names=FALSE, col.names = TRUE)

# install genefilter library to use rowttest() funtion
library(genefilter)
# install packages for cross validation
library(caret)
library(pROC)


set.seed(123)

trainPatients <- balanced.sample(Gbm.Survival.TCGA.Data, 0.7) 
testPatients <- setdiff(Gbm.Patients, trainPatients)

trainSet <- Gbm.RPPA.TCGA.DATA[ ,match (trainPatients, colnames(Gbm.RPPA.TCGA.DATA))]
testSet <- Gbm.RPPA.TCGA.DATA[ ,match (testPatients, colnames(Gbm.RPPA.TCGA.DATA))]


#create label for ttest
trainSet.label  <- factor(Gbm.Survival.TCGA.Data[match (colnames(trainSet), Gbm.Survival.TCGA.Data[,1]), 2])
testSet.label  <- factor(Gbm.Survival.TCGA.Data[match (colnames(testSet), Gbm.Survival.TCGA.Data[,1]), 2])

#apply ttest
trainSet.tt <- rowttests(as.matrix(trainSet), trainSet.label)
trainSet.pval <- cbind(trainSet, trainSet.tt$p.value)

trainSet.tt.statistic <- abs(trainSet.tt["statistic"])
trainSet.tt.dm <- abs(trainSet.tt["dm"])


install.packages("prob")
require("prob")
trainSet.tt.probability <- probspace(trainSet.tt.statistic)["probs"]

#testSet.tt <- rowttests(as.matrix(testSet), testSet.label)
#testSet.pval <- cbind(testSet, testSet.tt$p.value)



#CORRECT THE REST OF THE CODE


trainSet.filtered <- subset(trainSet.pval, trainSet.tt$p.value < 0.05)
testSet.filtered <- subset(testSet.pval, trainSet.tt$p.value < 0.05)

lastCol <- dim(trainSet.filtered)[2]
preTrain <- as.data.frame(trainSet.filtered [,-lastCol])
transposeTrain <- t(preTrain)
svmTrainData <- data.frame(transposeTrain, survivalClass = trainSet.label)
#names(svmTrainData) <- gsub("X", "", names(svmTrainData))

lastCol <- dim(testSet.filtered)[2]
preTest <- testSet.filtered [,-lastCol] 
transposeTest <- t(preTest)
svmTestData <- data.frame(transposeTest, survivalClass = testSet.label)
#names(svmTestData) <- gsub("X", "", names(svmTestData))


# "repeatedcv", is used to specify
# repeated K–fold cross–validation (and the argument repeats controls the number of repetitions).
# K is controlled by the number argument and defaults to 10.
# to change k use this: trainControl(number=100)
ctrl <- trainControl(method = "repeatedcv",repeats = 10,classProbs = TRUE,summaryFunction = twoClassSummary)
plsFit <- train(survivalClass ~ .,data = svmTrainData, method = "pls", tuneLength = 15,trControl = ctrl,metric = "ROC", preProc = c("center", "scale"))
plsClasses <- predict(plsFit, newdata = svmTestData)
str(plsClasses)
plsProbs <- predict(plsFit, newdata = svmTestData, type = "prob")
head(plsProbs)

confMatrix = confusionMatrix(data = plsClasses, svmTestData$survivalClass)
print(confMatrix$table)
print(confMatrix$overall)
print(confMatrix$byClass)

plot(plsFit) 

setwd("C:/Users/Ece/Desktop/GBM")
write.table(confMatrix$overall, file = "results.csv", sep = ",",col.names = NA)
write.table(confMatrix$table, file = "results.csv", sep = ",",col.names = NA,append=TRUE) 

