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

#trainSet.tt.statistic <- abs(trainSet.tt["statistic"])
#trainSet.tt.dm <- abs(trainSet.tt["dm"])
#install.packages("prob")
#require("prob")
#trainSet.tt.probability <- probspace(trainSet.tt.statistic)["probs"]

#testSet.tt <- rowttests(as.matrix(testSet), testSet.label)
#testSet.pval <- cbind(testSet, testSet.tt$p.value)


trainSet.tt.stat <- abs(trainSet.tt$statistic)

trainSet.tt.dm <- abs(trainSet.tt$dm)
prob.data.dm <- as.data.frame(cbind(Gbm.RPPA$EntrezID, trainSet.tt.stat))
colnames(prob.data) <- c("EntrezID","tt.stat")


prob.data <- as.data.frame(cbind(Gbm.RPPA$EntrezID, trainSet.tt.stat))
colnames(prob.data) <- c("EntrezID","tt.stat")

library(prob)
library(igraph)

gS <- read.graph("graph.txt", format="ncol",directed = FALSE)

length(V(gS))
# 10579
length(E(gS))
# 200091

dg = degree(gS)
degree <- as.data.frame(sort(dg, decreasing=TRUE))

vertex.name <- as.numeric(as.vector(V(gS)$name))
inp <- rep(0,length(V(gS)))

for(i in 1:dim(prob.data)[1])
{  
  a = which ((prob.data$EntrezID[i]) == vertex.name)
  if (length(a)>0){
    inp[a] <- prob.data[i,2]
  }
}

#create a probability vector for Personalized PageRank
prob.vector <- probspace(inp, inp)$probs


gS.rank <- page_rank(gS, vids = V(gS), directed=FALSE, damping=0.1)$vector      # returns numeric vector with rank scores
rank.top100 <- as.data.frame(head(sort(gS.rank, decreasing=TRUE), 100))

for(i in 1:dim(rank.top100)[1])
{	
  id = as.numeric(rownames(rank.top100)[i])
  a = which (id == as.numeric(rownames(degree)))
  write.table(cbind(id, degree[a,1], rank.top100[i,1]), file = "rank.top100.d-0.1.txt", append=TRUE, quote=FALSE, sep ="\t", row.names=FALSE,col.names=FALSE)
}


gS.person <- page_rank(gS, vids = V(gS), directed=FALSE, damping=0.1, personalized=prob.vector)$vector      
person.top100 <- as.data.frame(head(sort(gS.person, decreasing=TRUE), 100))

for(i in 1:dim(person.top100)[1])
{	
  id = as.numeric(rownames(person.top100)[i])
  a = which (id == as.numeric(rownames(degree)))
  write.table(cbind(id, degree[a,1], person.top100[i,1]), file = "person.top100.d-0.1.txt", append=TRUE, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
}

# only takr top-100 most significant proteins base on tt.statistic

sorted.rppa.100 <- head(prob.data[order(prob.data[,2], decreasing = TRUE),],100)

#------------------------------------------------------------------------------------------------
rank.data <- read.table("rank.top100.d-0.1.txt", header = FALSE, sep = "\t")
person.data <- read.table("person.top100.d-0.1.txt", header = FALSE, sep = "\t")

common.sig.rank <- intersect (rank.data[,1], prob.data$EntrezID)
length(common.sig.rank)
# 26

common.sig.person <- intersect (person.data[,1], prob.data$EntrezID)
length(common.sig.person)
# 100


common.sig.rank <- intersect (rank.data[,1], sorted.rppa.100$EntrezID)
length(common.sig.rank)
# 17

common.sig.person <- intersect (person.data[,1], sorted.rppa.100$EntrezID)
length(common.sig.person)
# 97
# 97 rppa sig. proteins are also listed in the top-ranked proteins of Personalized PageRank

#------------------------------------------------------------------------------------------------
rank.data <- read.table("rank.top100.d-0.5.txt", header = FALSE, sep = "\t")
person.data <- read.table("person.top100.d-0.5.txt", header = FALSE, sep = "\t")

common.sig.rank <- intersect (rank.data[,1], prob.data$EntrezID)
length(common.sig.rank)
# 30

common.sig.person <- intersect (person.data[,1], prob.data$EntrezID)
length(common.sig.person)
# 98


common.sig.rank <- intersect (rank.data[,1], sorted.rppa.100$EntrezID)
length(common.sig.rank)
# 18

common.sig.person <- intersect (person.data[,1], sorted.rppa.100$EntrezID)
length(common.sig.person)
# 93

#------------------------------------------------------------------------------------------------
rank.data <- read.table("rank.top100.d-0.8.txt", header = FALSE, sep = "\t")
person.data <- read.table("person.top100.d-0.8.txt", header = FALSE, sep = "\t")

common.sig.rank <- intersect (rank.data[,1], prob.data$EntrezID)
length(common.sig.rank)
# 32

common.sig.person <- intersect (person.data[,1], prob.data$EntrezID)
length(common.sig.person)
# 90


common.sig.rank <- intersect (rank.data[,1], sorted.rppa.100$EntrezID)
length(common.sig.rank)
# 21

common.sig.person <- intersect (person.data[,1], sorted.rppa.100$EntrezID)
length(common.sig.person)
# 81


## How many proteins are in common on the top-10 sig. proteins?

sorted.rppa.10 <- head(prob.data[order(prob.data[,2], decreasing = TRUE),],10)

common.sig.rank <- intersect (rank.data[1:10,1], sorted.rppa.10$EntrezID)
length(common.sig.rank)
# 0

common.sig.person <- intersect (person.data[1:10,1], sorted.rppa.10$EntrezID)
length(common.sig.person)
# 7

## How many proteins are in common on the top-20 sig. proteins?
sorted.rppa.20 <- head(prob.data[order(prob.data[,2], decreasing = TRUE),],20)

common.sig.rank <- intersect (rank.data[1:20,1], sorted.rppa.20$EntrezID)
length(common.sig.rank)
# 1

common.sig.person <- intersect (person.data[1:20,1], sorted.rppa.20$EntrezID)
length(common.sig.person)
# 15


## How many proteins are in common between the top-ranked-100 and all sig. proteins?

common.sig.rank <- intersect (rank.data[1:10,1], sorted.rppa$EntrezID)
length(common.sig.rank)
#2

common.sig.person <- intersect (person.data[1:10,1], sorted.rppa$EntrezID)
length(common.sig.person)
#9

##### For tt.dm
trainSet.tt.dm <- abs(trainSet.tt$dm)
prob.data.dm <- as.data.frame(cbind(Gbm.RPPA$EntrezID, trainSet.tt.dm))
colnames(prob.data.dm) <- c("EntrezID","tt.dm")

library(prob)
library(igraph)

setwd("C:/Users/Ece/Desktop")
gS <- read.graph("graph.txt", format="ncol",directed = FALSE)

length(V(gS))
# 10579
length(E(gS))
# 200091

dg = degree(gS)
degree <- as.data.frame(sort(dg, decreasing=TRUE))

vertex.name <- as.numeric(as.vector(V(gS)$name))
#zero or dm value
inp.dm <- rep(0,length(V(gS)))

for(i in 1:dim(prob.data.dm)[1])
{  
  match = which ((prob.data.dm$EntrezID[i]) == vertex.name)
  if (length(match)>0){
    inp.dm[match] <- prob.data.dm[i,2]
  }
}

#create a probability vector for Personalized PageRank
prob.vector.dm <- probspace(inp.dm, inp.dm)$probs


gS.rank <- page_rank(gS, vids = V(gS), directed=FALSE, damping=0.1)$vector      # returns numeric vector with rank scores
rank.top100 <- as.data.frame(head(sort(gS.rank, decreasing=TRUE), 100))

for(i in 1:dim(rank.top100)[1])
{  
  id = as.numeric(rownames(rank.top100)[i])
  a = which (id == as.numeric(rownames(degree)))
  write.table(cbind(id, degree[a,1], rank.top100[i,1]), file = 
                "rank.top100.dm.d-0.1.txt", append=TRUE, quote=FALSE, sep ="\t", row.names=FALSE,col.names=FALSE)
}


gS.person.dm <- page_rank(gS, vids = V(gS), directed=FALSE, damping=0.1, personalized=prob.vector.dm)$vector      
person.top100.dm <- as.data.frame(head(sort(gS.person.dm, decreasing=TRUE), 100))

for(i in 1:dim(person.top100.dm)[1])
{  
  id = as.numeric(rownames(person.top100.dm)[i])
  a = which (id == as.numeric(rownames(degree)))
  write.table(cbind(id, degree[a,1], person.top100.dm[i,1]), file = 
                "person.top100.dm.d-0.1.txt", append=TRUE, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
}

# only take top-100 most significant proteins base on tt.dm
sorted.rppa.100.dm <- head(prob.data.dm[order(prob.data.dm[,2], decreasing = TRUE),],100)

#------------------------------------------------------------------------------------------------
rank.data.dm <- read.table("rank.top100.dm.d-0.1.txt", header = FALSE, sep = "\t")
person.data.dm <- read.table("person.top100.dm.d-0.1.txt", header = FALSE, sep = "\t")

common.sig.rank.dm <- intersect (rank.data.dm[,1], prob.data.dm$EntrezID)
length(common.sig.rank.dm)
# 26

common.sig.person.dm <- intersect (person.data.dm[,1], prob.data.dm$EntrezID)
length(common.sig.person.dm)
# 100


common.sig.rank.dm <- intersect (rank.data.dm[,1], sorted.rppa.100.dm$EntrezID)
length(common.sig.rank.dm)
# 17

common.sig.person.dm <- intersect (person.data.dm[,1], sorted.rppa.100.dm$EntrezID)
length(common.sig.person.dm)
# 98
# 98 rppa sig. proteins are also listed in the top-ranked proteins of Personalized PageRank

#------------------------------------------------------------------------------------------------
#for damping factor 0.5
gS.rank.05 <- page_rank(gS, vids = V(gS), directed=FALSE, damping=0.5)$vector      # returns numeric vector with rank scores
rank.top100.05 <- as.data.frame(head(sort(gS.rank.05, decreasing=TRUE), 100))

for(i in 1:dim(rank.top100.05)[1])
{  
  id = as.numeric(rownames(rank.top100.05)[i])
  a = which (id == as.numeric(rownames(degree)))
  write.table(cbind(id, degree[a,1], rank.top100.05[i,1]), file = 
                "rank.top100.dm.d-0.5.txt", append=TRUE, quote=FALSE, sep ="\t", row.names=FALSE,col.names=FALSE)
}


gS.person.dm.05 <- page_rank(gS, vids = V(gS), directed=FALSE, damping=0.5, personalized=prob.vector.dm)$vector      
person.top100.dm.05 <- as.data.frame(head(sort(gS.person.dm.05, decreasing=TRUE), 100))

for(i in 1:dim(person.top100.dm.05)[1])
{  
  id = as.numeric(rownames(person.top100.dm.05)[i])
  a = which (id == as.numeric(rownames(degree)))
  write.table(cbind(id, degree[a,1], person.top100.dm.05[i,1]), file = 
                "person.top100.dm.d-0.5.txt", append=TRUE, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
}

# only take top-100 most significant proteins base on tt.dm
sorted.rppa.100.dm.05 <- head(prob.data.dm[order(prob.data.dm[,2], decreasing = TRUE),],100)


rank.data.dm.05 <- read.table("rank.top100.dm.d-0.5.txt", header = FALSE, sep = "\t")
person.data.dm.05 <- read.table("person.top100.dm.d-0.5.txt", header = FALSE, sep = "\t")

common.sig.rank.dm.05 <- intersect (rank.data.dm.05[,1], prob.data.dm$EntrezID)
length(common.sig.rank.dm.05)
# 30

common.sig.person.dm.05 <- intersect (person.data.dm.05[,1], prob.data.dm$EntrezID)
length(common.sig.person.dm.05)
# 99


common.sig.rank.dm.05 <- intersect (rank.data.dm.05[,1], sorted.rppa.100.dm.05$EntrezID)
length(common.sig.rank.dm.05)
# 19

common.sig.person.dm.05 <- intersect (person.data.dm.05[,1], sorted.rppa.100.dm.05$EntrezID)
length(common.sig.person.dm.05)
# 93

#------------------------------------------------------------------------------------------------
#for damping factor 0.8
gS.rank.08 <- page_rank(gS, vids = V(gS), directed=FALSE, damping=0.8)$vector      # returns numeric vector with rank scores
rank.top100.08 <- as.data.frame(head(sort(gS.rank.08, decreasing=TRUE), 100))

for(i in 1:dim(rank.top100.08)[1])
{  
  id = as.numeric(rownames(rank.top100.08)[i])
  a = which (id == as.numeric(rownames(degree)))
  write.table(cbind(id, degree[a,1], rank.top100.08[i,1]), file = 
                "rank.top100.dm.d-0.8.txt", append=TRUE, quote=FALSE, sep ="\t", row.names=FALSE,col.names=FALSE)
}


gS.person.dm.08 <- page_rank(gS, vids = V(gS), directed=FALSE, damping=0.8, personalized=prob.vector.dm)$vector      
person.top100.dm.08 <- as.data.frame(head(sort(gS.person.dm.08, decreasing=TRUE), 100))

for(i in 1:dim(person.top100.dm.08)[1])
{  
  id = as.numeric(rownames(person.top100.dm.08)[i])
  a = which (id == as.numeric(rownames(degree)))
  write.table(cbind(id, degree[a,1], person.top100.dm.08[i,1]), file = 
                "person.top100.dm.d-0.8.txt", append=TRUE, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
}

# only take top-100 most significant proteins base on tt.dm
sorted.rppa.100.dm.08 <- head(prob.data.dm[order(prob.data.dm[,2], decreasing = TRUE),],100)

rank.data.dm.08 <- read.table("rank.top100.dm.d-0.8.txt", header = FALSE, sep = "\t")
person.data.dm.08 <- read.table("person.top100.dm.d-0.8.txt", header = FALSE, sep = "\t")

common.sig.rank.dm.08 <- intersect (rank.data.dm.08[,1], prob.data.dm$EntrezID)
length(common.sig.rank.dm.08)
# 32

common.sig.person.dm.08 <- intersect (person.data.dm.08[,1], prob.data.dm$EntrezID)
length(common.sig.person.dm.08)
# 91


common.sig.rank.dm.08 <- intersect (rank.data.dm.08[,1], sorted.rppa.100.dm.08$EntrezID)
length(common.sig.rank.dm.08)
# 22

common.sig.person.dm.08 <- intersect (person.data.dm.08[,1], sorted.rppa.100.dm.08$EntrezID)
length(common.sig.person.dm.08)
# 83


## How many proteins are in common on the top-10 sig. proteins?

sorted.rppa.10.dm <- head(prob.data.dm[order(prob.data.dm[,2], decreasing = TRUE),],10)

# for d=0.1
common.sig.rank.dm <- intersect (rank.data.dm[1:10,1], sorted.rppa.10.dm$EntrezID)
length(common.sig.rank.dm )
# 1

# for d=0.5
common.sig.rank.dm.05 <- intersect (rank.data.dm.05[1:10,1], sorted.rppa.10.dm$EntrezID)
length(common.sig.rank.dm.05)
# 1

# for d=0.8
common.sig.rank.dm.08 <- intersect (rank.data.dm.08[1:10,1], sorted.rppa.10.dm$EntrezID)
length(common.sig.rank.dm.08)
# 0

#for d=0.1
common.sig.person.dm <- intersect (person.data.dm[1:10,1], sorted.rppa.10.dm$EntrezID)
length(common.sig.person.dm)
# 9

#for d=0.5
common.sig.person.dm.05 <- intersect (person.data.dm.05[1:10,1], sorted.rppa.10.dm$EntrezID)
length(common.sig.person.dm.05)
# 8

#for d=0.8
common.sig.person.dm.08 <- intersect (person.data.dm.08[1:10,1], sorted.rppa.10.dm$EntrezID)
length(common.sig.person.dm.08)
# 7

## How many proteins are in common on the top-20 sig. proteins?
sorted.rppa.20.dm <- head(prob.data.dm[order(prob.data.dm[,2], decreasing = TRUE),],20)

#for d=0.1
common.sig.rank.dm <- intersect (rank.data.dm[1:20,1], sorted.rppa.20.dm$EntrezID)
length(common.sig.rank.dm)
# 1

#for d=0.5
common.sig.rank.dm.05 <- intersect (rank.data.dm.05[1:20,1], sorted.rppa.20.dm$EntrezID)
length(common.sig.rank.dm.05)
# 1

#for d=0.8
common.sig.rank.dm.08 <- intersect (rank.data.dm.08[1:20,1], sorted.rppa.20.dm$EntrezID)
length(common.sig.rank.dm.08)
# 1

#for d=0.1
common.sig.person.dm <- intersect (person.data.dm[1:20,1], sorted.rppa.20.dm$EntrezID)
length(common.sig.person.dm)
# 19

#for d=0.5
common.sig.person.dm.05 <- intersect (person.data.dm.05[1:20,1], sorted.rppa.20.dm$EntrezID)
length(common.sig.person.dm.05)
# 18

#for d=0.8
common.sig.person.dm.08 <- intersect (person.data.dm.08[1:20,1], sorted.rppa.20.dm$EntrezID)
length(common.sig.person.dm.08)
# 14


## How many proteins are in common between the top-ranked-100 and all sig. proteins?
#for d= 0.1
common.sig.rank.dm <- intersect (rank.data.dm[1:10,1], sorted.rppa.100.dm$EntrezID)
length(common.sig.rank.dm)
#2

#for d= 0.5
common.sig.rank.dm.05 <- intersect (rank.data.dm.05[1:10,1], sorted.rppa.100.dm$EntrezID)
length(common.sig.rank.dm.05)
#3


#for d= 0.8
common.sig.rank.dm.08 <- intersect (rank.data.dm.08[1:10,1], sorted.rppa.100.dm$EntrezID)
length(common.sig.rank.dm.08)
#2

#for d=0.1
common.sig.person.dm <- intersect (person.data.dm[1:10,1], sorted.rppa.100.dm$EntrezID)
length(common.sig.person.dm)
#10


#for d=0.5
common.sig.person.dm.05 <- intersect (person.data.dm.05[1:10,1], sorted.rppa.100.dm$EntrezID)
length(common.sig.person.dm.05)
#10


#for d=0.8
common.sig.person.dm.08 <- intersect (person.data.dm.08[1:10,1], sorted.rppa.100.dm$EntrezID)
length(common.sig.person.dm.08)
#9





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

