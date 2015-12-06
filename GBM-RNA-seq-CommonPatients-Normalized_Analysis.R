# Convert txt file to csv file
setwd("C:/Users/Ece/Desktop/GBM")
tab = read.delim("GBM-RNA-seq-CommonPatients-Normalized.txt")
write.table(tab, file="GBM-RNA-seq-CommonPatients-Normalized.csv",sep=",",col.names=TRUE,row.names=FALSE)


balanced.sample <- function(pData, n) {
  
  all.good = which(pData$z == "C1")
  all.poor = which(pData$z == "C0")
  
  sample.good = sample(all.good, round(length(all.good)*n), replace=FALSE)
  sample.poor = sample(all.poor, round(length(all.poor)*n), replace=FALSE)
  
  training.set = c(pData[sample.good, 1], pData[sample.poor, 1])  # return patient ids -> TCGA.06.2567
  return(training.set)
}

# read GBM RNA Seq
Gbm.RNASeq = read.csv("GBM-RNA-seq-CommonPatients-Normalized.csv",header=T)

# get only TCGA columns
TCGA.start = Gbm.RNASeq[4:38]
Gbm.RNASeq.TCGA.DATA <- TCGA.start

Gbm.RNASeq.patient = colnames(TCGA.start)

# read GBM Survival 
Gbm.Survival = read.csv("GBM-survival-class.csv",header=T)
# replace "-" with "." to match RNASeq TCGA data (TCGA-12-1098 => TCGA.12.1098)
Gbm.Survival$feature = gsub("-",".",as.factor(Gbm.Survival$feature))
Gbm.Survival.RNASeq= Gbm.Survival$feature


# find same TCGA's = patients
Gbm.Patients = intersect(c(Gbm.Survival.RNASeq), c(Gbm.RNASeq.patient))

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

# install genefilter library to use rowttest() funtion
library(genefilter)
# install packages for cross validation
library(caret)
library(pROC)


set.seed(123)

trainPatients <- balanced.sample(Gbm.Survival.TCGA.Data, 0.7) 
testPatients <- setdiff(Gbm.Patients, trainPatients)

trainSet <- Gbm.RNASeq.TCGA.DATA[ ,match (trainPatients, colnames(Gbm.RNASeq.TCGA.DATA))]
testSet <- Gbm.RNASeq.TCGA.DATA[ ,match (testPatients, colnames(Gbm.RNASeq.TCGA.DATA))]


#create label for ttest
trainSet.label  <- factor(Gbm.Survival.TCGA.Data[match (colnames(trainSet), Gbm.Survival.TCGA.Data[,1]), 2])
testSet.label  <- factor(Gbm.Survival.TCGA.Data[match (colnames(testSet), Gbm.Survival.TCGA.Data[,1]), 2])

#apply ttest
trainSet.tt <- rowttests(as.matrix(trainSet), trainSet.label)
trainSet.pval <- cbind(trainSet, trainSet.tt$p.value)

trainSet.tt.stat <- abs(trainSet.tt$statistic)

trainSet.tt.dm <- abs(trainSet.tt$dm)
prob.data.dm <- as.data.frame(cbind(Gbm.RNASeq$EntrezID, trainSet.tt.stat))

prob.data <- as.data.frame(cbind(Gbm.RNASeq$EntrezID, trainSet.tt.stat))
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
  write.table(cbind(id, degree[a,1], rank.top100[i,1]), file = "RNA-seq.rank.top100.d-0.1.txt", append=TRUE, quote=FALSE, sep ="\t", row.names=FALSE,col.names=FALSE)
}


gS.person <- page_rank(gS, vids = V(gS), directed=FALSE, damping=0.1, personalized=prob.vector)$vector      
person.top100 <- as.data.frame(head(sort(gS.person, decreasing=TRUE), 100))

for(i in 1:dim(person.top100)[1])
{	
  id = as.numeric(rownames(person.top100)[i])
  a = which (id == as.numeric(rownames(degree)))
  write.table(cbind(id, degree[a,1], person.top100[i,1]), 
              file = "RNA-seq.person.top100.d-0.1.txt", append=TRUE, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
}

# only take top-100 most significant proteins base on tt.statistic
sorted.rnaseq.100 <- head(prob.data[order(prob.data[,2], decreasing = TRUE),],100)

#------------------------------------------------------------------------------------------------
rank.data <- read.table("RNA-seq.rank.top100.d-0.1.txt", header = FALSE, sep = "\t")
person.data <- read.table("RNA-seq.person.top100.d-0.1.txt", header = FALSE, sep = "\t")

common.sig.rank <- intersect (rank.data[,1], prob.data$EntrezID)
length(common.sig.rank)
# 99

common.sig.person <- intersect (person.data[,1], prob.data$EntrezID)
length(common.sig.person)
# 100


common.sig.rank <- intersect (rank.data[,1], sorted.rnaseq.100$EntrezID)
length(common.sig.rank)
# 2

common.sig.person <- intersect (person.data[,1], sorted.rnaseq.100$EntrezID)
length(common.sig.person)
# 42
# 42 rnaseq sig. proteins are also listed in the top-ranked proteins of Personalized PageRank

#------------------------------------------------------------------------------------------------

#for damping factor 0.5
gS.rank.05 <- page_rank(gS, vids = V(gS), directed=FALSE, damping=0.5)$vector      # returns numeric vector with rank scores
rank.top100.05 <- as.data.frame(head(sort(gS.rank.05, decreasing=TRUE), 100))

for(i in 1:dim(rank.top100.05)[1])
{  
  id = as.numeric(rownames(rank.top100.05)[i])
  a = which (id == as.numeric(rownames(degree)))
  write.table(cbind(id, degree[a,1], rank.top100.05[i,1]), file = 
                "RNA-seq.rank.top100.d-0.5.txt", append=TRUE, quote=FALSE, sep ="\t", row.names=FALSE,col.names=FALSE)
}


gS.person.05 <- page_rank(gS, vids = V(gS), directed=FALSE, damping=0.5, personalized=prob.vector)$vector      
person.top100.05 <- as.data.frame(head(sort(gS.person.05, decreasing=TRUE), 100))

for(i in 1:dim(person.top100.05)[1])
{  
  id = as.numeric(rownames(person.top100.05)[i])
  a = which (id == as.numeric(rownames(degree)))
  write.table(cbind(id, degree[a,1], person.top100.05[i,1]), 
              file = "RNA-seq.person.top100.d-0.5.txt", append=TRUE, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
}

# only take top-100 most significant proteins base on tt.statistic
sorted.rnaseq.100 <- head(prob.data[order(prob.data[,2], decreasing = TRUE),],100)


rank.data.05 <- read.table("RNA-seq.rank.top100.d-0.5.txt", header = FALSE, sep = "\t")
person.data.05 <- read.table("RNA-seq.person.top100.d-0.5.txt", header = FALSE, sep = "\t")

common.sig.rank.05 <- intersect (rank.data.05[,1], prob.data$EntrezID)
length(common.sig.rank.05)
# 99

common.sig.person.05 <- intersect (person.data.05[,1], prob.data$EntrezID)
length(common.sig.person.05)
# 183


common.sig.rank.05 <- intersect (rank.data.05[,1], sorted.rnaseq.100$EntrezID)
length(common.sig.rank.05)
# 2

common.sig.person.05 <- intersect (person.data.05[,1], sorted.rnaseq.100$EntrezID)
length(common.sig.person.05)
# 42

#------------------------------------------------------------------------------------------------

#for damping factor 0.8
gS.rank.08 <- page_rank(gS, vids = V(gS), directed=FALSE, damping=0.8)$vector      # returns numeric vector with rank scores
rank.top100.08 <- as.data.frame(head(sort(gS.rank.08, decreasing=TRUE), 100))

for(i in 1:dim(rank.top100.08)[1])
{  
  id = as.numeric(rownames(rank.top100.08)[i])
  a = which (id == as.numeric(rownames(degree)))
  write.table(cbind(id, degree[a,1], rank.top100.08[i,1]), file = 
                "RNA-seq.rank.top100.d-0.8.txt", append=TRUE, quote=FALSE, sep ="\t", row.names=FALSE,col.names=FALSE)
}


gS.person.08 <- page_rank(gS, vids = V(gS), directed=FALSE, damping=0.8, personalized=prob.vector)$vector      
person.top100.08 <- as.data.frame(head(sort(gS.person.08, decreasing=TRUE), 100))

for(i in 1:dim(person.top100.08)[1])
{  
  id = as.numeric(rownames(person.top100.08)[i])
  a = which (id == as.numeric(rownames(degree)))
  write.table(cbind(id, degree[a,1], person.top100.08[i,1]), 
              file = "RNA-seq.person.top100.d-0.8.txt", append=TRUE, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
}

# only take top-100 most significant proteins base on tt.statistic
sorted.rnaseq.100 <- head(prob.data[order(prob.data[,2], decreasing = TRUE),],100)

rank.data.08 <- read.table("RNA-seq.rank.top100.d-0.8.txt", header = FALSE, sep = "\t")
person.data.08 <- read.table("RNA-seq.person.top100.d-0.8.txt", header = FALSE, sep = "\t")

common.sig.rank.08 <- intersect (rank.data.08[,1], prob.data$EntrezID)
length(common.sig.rank.08)
# 99

common.sig.person.08 <- intersect (person.data.08[,1], prob.data$EntrezID)
length(common.sig.person.08)
# 99

common.sig.rank.08 <- intersect (rank.data.08[,1], sorted.rnaseq.100$EntrezID)
length(common.sig.rank.08)
# 2

common.sig.person.08 <- intersect (person.data.08[,1], sorted.rnaseq.100$EntrezID)
length(common.sig.person.08)
# 2

##------------------------------------------------------------------------

## How many proteins are in common on the top-10 sig. proteins?

sorted.rnaseq.10 <- head(prob.data[order(prob.data[,2], decreasing = TRUE),],10)

# for d=0.1
common.sig.rank <- intersect (rank.data[1:10,1], sorted.rnaseq.10$EntrezID)
length(common.sig.rank)
# 1

# for d=0.5
common.sig.rank.05 <- intersect (rank.data.05[1:10,1], sorted.rnaseq.10$EntrezID)
length(common.sig.rank.05)
# 1

# for d=0.8
common.sig.rank.08 <- intersect (rank.data.08[1:10,1], sorted.rnaseq.10$EntrezID)
length(common.sig.rank.08)
# 1

#for d=0.1
common.sig.person <- intersect (person.data[1:10,1], sorted.rnaseq.10$EntrezID)
length(common.sig.person)
# 1

#for d=0.5
common.sig.person.05 <- intersect (person.data.05[1:10,1], sorted.rnaseq.10$EntrezID)
length(common.sig.person.05)
# 8

#for d=0.8
common.sig.person.08 <- intersect (person.data.08[1:10,1], sorted.rnaseq.10$EntrezID)
length(common.sig.person.08)
# 1

## How many proteins are in common on the top-20 sig. proteins?
sorted.rnaseq.20 <- head(prob.data[order(prob.data[,2], decreasing = TRUE),],20)

#for d=0.1
common.sig.rank <- intersect (rank.data[1:20,1], sorted.rnaseq.20$EntrezID)
length(common.sig.rank)
# 1

#for d=0.5
common.sig.rank.05 <- intersect (rank.data.05[1:20,1], sorted.rnaseq.20$EntrezID)
length(common.sig.rank.05)
# 1

#for d=0.8
common.sig.rank.08 <- intersect (rank.data.08[1:20,1], sorted.rnaseq.20$EntrezID)
length(common.sig.rank.08)
# 1

#for d=0.1
common.sig.person <- intersect (person.data[1:20,1], sorted.rnaseq.20$EntrezID)
length(common.sig.person)
# 11

#for d=0.5
common.sig.person.05 <- intersect (person.data.05[1:20,1], sorted.rnaseq.20$EntrezID)
length(common.sig.person.05)
# 11

#for d=0.8
common.sig.person.08 <- intersect (person.data.08[1:20,1], sorted.rnaseq.20$EntrezID)
length(common.sig.person.08)
# 1


## How many proteins are in common between the top-ranked-100 and all sig. proteins?
#for d= 0.1
common.sig.rank <- intersect (rank.data[1:10,1], sorted.rnaseq.100$EntrezID)
length(common.sig.rank)
# 1

#for d= 0.5
common.sig.rank.05 <- intersect (rank.data.05[1:10,1], sorted.rnaseq.100$EntrezID)
length(common.sig.rank.05)
# 1


#for d= 0.8
common.sig.rank.08 <- intersect (rank.data.08[1:10,1], sorted.rnaseq.100$EntrezID)
length(common.sig.rank.08)
# 1

#for d=0.1
common.sig.person <- intersect (person.data[1:10,1], sorted.rnaseq.100$EntrezID)
length(common.sig.person)
#9


#for d=0.5
common.sig.person.05 <- intersect (person.data.05[1:10,1], sorted.rnaseq.100$EntrezID)
length(common.sig.person.05)
# 1


#for d=0.8
common.sig.person.08 <- intersect (person.data.08[1:10,1], sorted.rnaseq.100$EntrezID)
length(common.sig.person.08)
# 1

