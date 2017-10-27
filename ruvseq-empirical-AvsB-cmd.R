#setwd("/Users/stalenygard/ruvseq")
#install.packages("pROC")
library(pROC)

files<-c("GSE47774_SEQC_ILM_AGR.txt",
"GSE47774_SEQC_ILM_BGI.txt",
"GSE47774_SEQC_ILM_CNL.txt",
"GSE47774_SEQC_ILM_COH.txt",
"GSE47774_SEQC_ILM_MAY.txt",
"GSE47774_SEQC_ILM_NVS.txt",
"GSE47774_SEQC_ROC_SQW.txt",
"GSE47774_SEQC_ROC_NYU.txt",
"GSE47774_SEQC_ROC_MGP.txt",
"GSE47774_SEQC_LIF_SQW.txt",
"GSE47774_SEQC_LIF_PSU.txt",
"GSE47774_SEQC_LIF_NWU.txt",
"GSE47774_SEQC_LIF_LIV.txt")

D<-read.table(files[1],sep="\t",header=TRUE)
#for (i in 2:length(files)){
#    print(i)
#    data<-read.table(files[i],sep="\t",header=TRUE)
#    D<-cbind(D,data[,-1])
#}

dim(D)
ix.1<-grep("_A_",names(D))
ix.2<-grep("_B_",names(D))

#source("https://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
library(DESeq2)
ix<-c(ix.1,ix.2)
X<-D[,ix]
X=as.matrix(X)
condition=factor(c(rep("C",length(ix.1)),rep("D",length(ix.2))))
coldata<-data.frame(condition)
row.names(coldata)<-colnames(X)

#source("https://bioconductor.org/biocLite.R")
#library(RUVSeq)
#filter <- apply(X, 1, function(x), length(x[x>5])>=2)
#filtered <- X[filter,]
#set <- newSeqExpressionSet(as.matrix(filtered),phenoData =coldata)


#install.packages("pROC")
library(pROC)
source("https://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
library(DESeq2)

dds <- DESeqDataSetFromMatrix(countData = X,
                              colData = coldata,
                              design= ~ condition)            

library(BiocParallel)
register(MulticoreParam(7))
#BiocParallel::register(BiocParallel::SerialParam())

dds<-DESeq(dds,parallel=TRUE)
res <- results(dds,parallel=TRUE)
#hkg<-read.table("hkg-short-2.txt",sep="\t",header=TRUE)
#hkg.navn<-hkg[,2]
navn<-as.character(D[,1])
data<-read.table("cms_095046.txt",sep="\t",header=TRUE)
mid.NEG<-as.character(data[which(data[,3]=="B"),2])
ix.neg<-match(gsub("-","_",mid.NEG),navn)
mid.POS<-as.character(data[which(data[,3]!="B"),2])
ix.pos<-match(gsub("-","_",mid.POS),navn)
ix.ctrl<-c(ix.neg,ix.pos)
pred1<-abs(res$log2FoldChange[ix.ctrl])
truth<-c(rep(0,length(ix.neg)),rep(1,length(ix.pos)))
roc.res.1<-roc(response=truth,predictor=pred1)
auc(roc.res.1)
plot(roc.res.1)

rm(dds)

#source("https://bioconductor.org/biocLite.R")
#biocLite("RUVSeq")
library(RUVSeq)
filter <- apply(X,1, function(x) length(x[x>5])>=2)
filtered <- X[filter,]

set <- newSeqExpressionSet(as.matrix(filtered),phenoData =coldata)

#establishing empirical
library(edgeR)
design <- model.matrix(~condition, data=pData(set))
y <- DGEList(counts=counts(set), group=coldata$condition)
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)
top <- topTags(lrt, n=nrow(set))$table
empirical <- rownames(set)[which(!(rownames(set) %in% rownames(top)[1:5000]))]

V<-data.frame(condition)
row.names(V)<-colnames(X)
colnames(V)<-"condition"
#rm(X)

L<-10
auc.vec<-rep(0,L) 
#source("https://bioconductor.org/biocLite.R")
#biocLite("RUVSeq")
library(RUVSeq)
for (l in 1:L){

	seqRUVg <- RUVg(X,empirical, k=l)
	#V<-data.frame(condition)
	col2<-cbind(V,seqRUVg$W)
	#col2<-as.matrix(col2)
	dds2 <- DESeqDataSetFromMatrix(countData = as.matrix(X),colData = col2,design=~condition)
	design(dds2)=as.formula(paste(c("~",paste(colnames(seqRUVg$W),rep("+" ,l)),"condition")))
	dds2 <- DESeq(dds2,parallel=TRUE)
	res2 <- results(dds2,parallel=TRUE)
	pred2<-abs(res2$log2FoldChange[ix.ctrl])
	roc.res.2<-roc(response=truth,predictor=pred2)
        auc.vec[l]<-auc(roc.res.2)
	print(auc.vec)
}

print(auc.vec)

pdf("SEQC-RUVSEQ-DESEQ2-emprical-AvsB-ercc-ROC-curves.pdf")
plot(roc.res.1)
plot(roc.res.2,col=2,add=TRUE)
graphics.off()

