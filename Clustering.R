# Author: Shrikant Pawar, 06/30/19, program for clustering cancer types.

# Following are the datasets used in this study.
Breast
GSE2034 
GSE25066

Ovarian
GSE9899
GSE26712-14764

Colon
GSE39582
GSE14333

Lung
GSE30219
GSE68465

# Load GEO and biobase libraries
library(GEOquery)
library(Biobase)
library(preprocessCore)
library(multiClust)

# Obtain GSE series matrix file from GEO website using getGEO function
gse <- getGEO(GEO="GSE30219")

# Save the gene expression matrix as an object
data.gse <- exprs(gse[[1]])

# Save the patient clinical data to an object
pheno <- pData(phenoData(gse[[1]]))

# Write the gene expression and clinical data to text files
WriteMatrixToFile(tmpMatrix=data.gse, tmpFileName="GSE30219.expression.txt",
blnRowNames=TRUE, blnColNames=TRUE)

WriteMatrixToFile(tmpMatrix=pheno, tmpFileName="GSE30219.clinical.txt",
blnRowNames=TRUE, blnColNames=TRUE)

# Obtain GSE series matrix file from GEO website using getGEo function
gse <- getGEO(GEO="GSE30219")

# Save the gene expression matrix as an object
data.gse <- exprs(gse[[1]])

# Quantile normalization of the dataset
data.norm <- normalize.quantiles(data.gse, copy=FALSE)

# shift data before log scaling to prevent errors from log scaling
# negative numbers
if (min(data.norm)> 0) {} else {mindata.norm=abs(min(data.norm)) + .001;
 data.norm=data.norm + mindata.norm}

# Log2 scaling of the dataset
data.log <- t(apply(data.norm, 1, log2))

# Write the gene expression and clinical data to text files
WriteMatrixToFile(tmpMatrix=data.log,
     tmpFileName="GSE30219.normalized.expression.txt",
     blnRowNames=TRUE, blnColNames=TRUE)

# Obtain gene expression matrix
exp_file <- system.file("extdata", "GSE30219.normalized.expression.txt", package= "multiClust")

# Load the gene expression matrix 
data.exprs <- input_file(input=exp_file)

OR

data.exprs <- read.table("GSE30219.normalized.expression.txt", header = TRUE, sep = "", dec = ".")
	
# View the first few rows and columns of the matrix
data.exprs[1:4,1:4]

--------------CAN MAKE OBJECTS WITH data.exprs FOR ALL CANCERS HERE-------------------------------

Breast <- data.exprs
54675   286
Ovarian <- data.exprs
Ovarian[is.na(Ovarian)] <- 0
54675   295
Colon <- data.exprs
54675   585
Lung <- data.exprs
54675   307

#Merge datasets to form one object with all the chips
OB <- merge(Ovarian, Breast, by=0, all=TRUE)
OBC <- merge(OB, Colon, by=0, all=TRUE)
OBCL <- merge(OBC, Lung, by=0, all=TRUE)
OBCL[is.na(OBCL)] <- 0 
OBCLnew <- OBCL[-c(1,3)]
OBCLnew[is.na(OBCLnew)] <- 0
dataclusterno <- OBCLnew[2:1474,]

dataclusterno2 <- as.data.frame(lapply(dataclusterno, as.numeric))
dataclusterno2[is.na(dataclusterno2)] <- 0

# Call the number_clusters function
## # data.exp is the original expression matrix object ouputted from
## # the input_file function
## # User chooses the gap_statistic option by making gap_statistic equal TRUE
## # The Fixed argument is also set to NULL
# OBTAIN CLUSTER NUMBER FROM BELOW

cluster_num <- number_clusters(data.exp=dataclusterno2, Fixed=NULL,
     gap_statistic=TRUE)

Clustering k = 1,2,..., K.max (= 8): .. done
Bootstrapping, b = 1,2,..., B (= 100)  [one "." per sample]:
.................................................. 50 
.................................................. 100 
[1] "The gap statistic cluster number is: 5"

suppressPackageStartupMessages(library(ctc))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(dendextend))
library(ctc)
library(gplots)
library(dendextend)
library(graphics)
library(grDevices)
library(amap)

# Call the cluster_analysis function
hclust_analysis <- cluster_analysis(sel.exp=dataclusterno2,
    cluster_type="HClust", seed=NULL,
    distance="euclidean", linkage_type="ward.D2", 
    gene_distance="correlation",
    num_clusters=5, data_name="Breast, Colon, Lung, Ovarian Cancer", 
    probe_rank="SD_Rank", probe_num_selection="Fixed_Probe_Num",
    cluster_num_selection="Fixed_Clust_Num")

# Display the first few columns and rows of the object
head(hclust_analysis)

# Produce a sample dendrogram using the selected gene/probe expression object
# In the function, cluster_analysis the dendrogram is not displayed in the console. 

# Make the selected expression object numeric
colname <- colnames(dataclusterno2)
ranked.exprs <- t(apply(dataclusterno2, 1,as.numeric))
colnames(ranked.exprs) <- colname

# Normalization of the selected expression matrix before clustering
norm.sel.exp <- t(apply(ranked.exprs, 1, nor.min.max))

hc <- hclust(dist(x=t(norm.sel.exp), method="euclidean"), method="ward.D2")
dc <- as.dendrogram(hc)

# Portray the samples in each cluster as a different color
dc <- set(dc, "branches_k_color", value=1:5, k=5)

# Make a vector filled with blank spaces to replace sample names with blank space
vec <- NULL
len <- length(labels(dc))
for (i in 1:len) {
  i <- ""
  vec <- c(vec, i)
}

# Remove sample labels from dendrogram
dc <- dendextend::set(dc, "labels", vec)

# Plot the sample dendrogram
plot(dc, main=paste("Four cancer types", "Sample Dendrogram", sep=" "))

# Display orange lines on the plot to clearly separate each cluster
rect.dendrogram(dc, k=3, border='orange')

# Call the cluster_analysis function for KMEANS ANALYSIS
kmeans_analysis <- cluster_analysis(sel.exp=dataclusterno2,
    cluster_type="Kmeans", seed=5,
    distance=NULL, linkage_type=NULL, 
    gene_distance=NULL, num_clusters=5,
    data_name="Breast, Colon, Lung, Ovarian Cancer", probe_rank="SD_Rank",
    probe_num_selection="Fixed_Probe_Num",
    cluster_num_selection="Fixed_Clust_Num")
 
# Display the first few rows and columns of the object
head(kmeans_analysis)

Hierarchial <- read.csv(file="Breast, Colon, Lung, Ovarian Cancer HClust euclidean ward.D2 SD_Rank Fixed_Probe_Num Fixed_Clust_Num Samples.Clusters.csv", header=TRUE, sep=",")
Kmeans <- read.csv(file="Breast, Colon, Lung, Ovarian Cancer Kmeans SD_Rank Fixed_Probe_Num Fixed_Clust_Num Samples.Clusters.csv", header=TRUE, sep=",")

Hierarchial <- Hierarchial[order(Hierarchial[,3]),] 
grps <- factor(Hierarchial[,3]) 

Hierarchial$color[Hierarchial[,3]=="O"] <- "red"
Hierarchial$color[Hierarchial[,3]=="B"] <- "blue"
Hierarchial$color[Hierarchial[,3]=="L"] <- "darkgreen" 
Hierarchial$color[Hierarchial[,3]=="C"] <- "black"

Hierarchial <- as.data.frame(Hierarchial)

# colors in order B C L O as follows
my_cols <- c("#000000", "#0000FF", "#FF0000", "#FFFF00")

dotchart(as.numeric(Hierarchial[,2])+1,labels=Hierarchial[,3],cex=.7,groups= grps,
   main="Cancers Clusters with Hierarchical Clustering",
   xlab="Cluster Number", gcolor=my_cols, color = my_cols[grps]
, xlim=c(1, 6), gpch = 15)

-------------PLOT FOR KMEANS----------------
Kmeans <- Kmeans[order(Kmeans[,3]),] 
grps <- factor(Kmeans[,3]) 

Kmeans$color[Kmeans[,3]=="O"] <- "red"
Kmeans$color[Kmeans[,3]=="B"] <- "blue"
Kmeans$color[Kmeans[,3]=="L"] <- "darkgreen" 
Kmeans$color[Kmeans[,3]=="C"] <- "black"

Kmeans <- as.data.frame(Kmeans)

Kmeans <- Kmeans[-1,]

# colors in order B C L O as follows
my_cols <- c("#000000", "#0000FF", "#FF0000", "#FFFF00")

dotchart(as.numeric(Kmeans[,2]),labels=Kmeans[,3],cex=.7,groups= grps,
   main="Cancers Clusters with K Means Clustering",
   xlab="Cluster Number", gcolor=my_cols, color = my_cols[grps]
, xlim=c(1, 6), gpch = 15)


# From above analysis Breast and Ovarian cancers are similar as the patients
# cluster almost in one cluster from both hierarchial and K means clustering.
# Next is get normal Breast and Ovarian samples, do the differential analysis
# collect the top upregulated and top downregulated genes, may be combine
# them and do RF analysis to identify unique genes selected in test data (3, 5, 10)
# use IMP function, follow the survival analysis for these unique genes to conclude.
# Can validate IMP genes with just fold change on Cancer/Normal

# Read in CANCER data files.
setwd("C:/Users/Bio-user/Desktop/Clustering_AI/Breast Cancer")
affy.data = ReadAffy()
eset.mas5breast = mas5(affy.data)
exprSet.nologsbreast = exprs(eset.mas5breast)
exprSetBreast = log(exprSet.nologsbreast, 2)
breastCancer <- rowMeans(exprSetBreast)
breastCancerMeans <- as.data.frame(breastCancer)

setwd("C:/Users/Bio-user/Desktop/Clustering_AI/Ovarian Cancer")
affy.data = ReadAffy()
eset.mas5ovarian = mas5(affy.data)
exprSet.nologsovarian = exprs(eset.mas5ovarian)
exprSetOvarian = log(exprSet.nologsovarian, 2)
ovarianCancer <- rowMeans(exprSetOvarian)
ovarianCancerMeans <- as.data.frame(ovarianCancer)

# Processing NORMAL samples
library(affy)
library(limma)
setwd("C:/Users/Bio-user/Desktop/Clustering_AI/Breast Normal Samples")
affy.data = ReadAffy()
eset.mas5normbreast = mas5(affy.data)
exprSet.nologsnormbreast = exprs(eset.mas5normbreast)
exprSetBreastnormbreast = log(exprSet.nologsnormbreast, 2)
BreastCancer <- rowMeans(exprSetBreastnormbreast)
BreastCancerMeans <- as.data.frame(BreastCancer)

setwd("C:/Users/Bio-user/Desktop/Clustering_AI/Normal Ovarian Samples")
affy.data = ReadAffy()
eset.mas5normovarian = mas5(affy.data)
exprSet.nologsnormovarian = exprs(eset.mas5normovarian)
exprSetOvariannormovarian = log(exprSet.nologsnormovarian, 2)
OvarianCancer <- rowMeans(exprSetOvariannormovarian)
OvarianCancerMeans <- as.data.frame(OvarianCancer)

# Foldchange to identify significant genes in cancer and normal samples
# 4 objects
# Cancer:
breastCancerMeans
ovarianCancerMeans
# Normal:
BreastCancerMeans
OvarianCancerMeans

BreastFinal <- merge(breastCancerMeans,BreastCancerMeans,by="row.names",all.x=TRUE)
OvarianFinal <- merge(ovarianCancerMeans,OvarianCancerMeans,by="row.names",all.x=TRUE)
BreastFinal[is.na(BreastFinal)] <- 0
OvarianFinal[is.na(OvarianFinal)] <- 0

OvarianFinal$FoldOvarian <- OvarianFinal[,2]-OvarianFinal[,3]
BreastFinal$FoldBreast <- BreastFinal[,2]-BreastFinal[,3]

write.csv(OvarianFinal, file = "OvarianFold.csv",row.names=FALSE)
write.csv(BreastFinal, file = "BreastFold.csv",row.names=FALSE)

-------APPLYING RF-----------------------------
library(randomForest)
set.seed(1234)
print(date())
RFData <- merge(exprSetBreast,exprSetOvarian,by="row.names",all.x=TRUE)
samp2 <- RFData[,-1]
rownames(samp2) <- RFData[,1]
write.csv(samp2, file = "RFData.csv",row.names=TRUE)
# classification, in order
group <- c(rep('B',286),rep('O',133)) 
samp2[is.na(samp2)] <- 0
print(date())
rf <- randomForest(x=t(samp2),y=as.factor(group),ntree=10000)
print(date())  
# Look at variable importance
imp.temp <- abs(rf$importance[,])
t <- order(imp.temp,decreasing=TRUE)
plot(c(1:nrow(samp2)),imp.temp[t],log='x',cex.main=1.5,
 xlab='Gene Rank',ylab='Variable Importance',cex.lab=1.5,
 pch=16,main='Variable Importance')

# Get subset of expression values for 25 most 'important' genes
gn.imp <- names(imp.temp)[t]
gn.25 <- gn.imp[1:25] # vector of top 25 genes, in order
t <- is.element(rownames(samp2),gn.25)
sig.eset <- samp2[t,] 

## Make a heatmap, with group differences obvious on plot
library(RColorBrewer)
hmcol <- colorRampPalette(brewer.pal(11,"PuOr"))(256)
colnames(sig.eset) <- group # This will label the heatmap columns
csc <- rep(hmcol[50],30)
csc[group=='T'] <- hmcol[200]
# column side color will be purple for T and orange for B
hmap <- data.matrix(sig.eset)
heatmap(hmap,scale="row", col=hmcol,ColSideColors=csc)
varImpPlot(rf, n.var=25, main='Variable Importance Results')
imp.temp <- importance(rf)
t <- order(imp.temp,decreasing=TRUE)
gn.imp <- names(imp.temp)[t]

#Surival Analysis for selected high/low expression genes in the patients







########################### CODE END HERE ###########################3#
