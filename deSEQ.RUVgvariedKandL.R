library(DESeq2)
library(edgeR)

###Clinical information with patients/controls as "condition"
colData=base.col
###Count matrix from RNA-Seq
counts=base.count
dim(count)
[1] 39195    47

##filter counts less than 5 in less than 2 samples
filter=apply(counts,1, function(x) length(x[x>5])>=2)
filtered.counts=count[filter,]

###Transformed measurement of the Age variable was obtained by using the scaling factors function before adding it to the DESeq2 model
ageF=cut(colData$Age,breaks=4)
colData$ageF=as.factor(ageF)

#####DESeq2 and RUVSeq with three options to analyse (A,B,C) as in figure 2
###A: DESeq2 Without variance removal
dds<- DESeqDataSetFromMatrix(countData=filtered.counts,
colData=colData,
design=~Condition)
design(dds)=formula(~ageF+gender+Condition)
test<-DESeq(dds)
res<-results(test)

summary(res)
###B: DESeq2 With RUVSeq variance removal
source("http://bioconductor.org/biocLite.R")
biocLite(RUVSeq)
library(RUVSeq)

#RUVg: obtain a set of “in-silico empirical” negative controls
set <- newSeqExpressionSet(as.matrix(filtered.counts),phenoData =colData)

design <- model.matrix(~Condition, data=pData(set))
y <- DGEList(counts=counts(set), group=colData$Condition)
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)
top <- topTags(lrt, n=nrow(set))$table
empirical <- rownames(top)[which(!(rownames(set) %in% rownames(top)[1:5000]))]
#range k=1-10
setk <- RUVg(set, empirical, k=k)
pData(setk)
plotRLE(setk, outline=FALSE, ylim=c(-4, 4), col=colors[colData$Condition])
plotPCA(setk, col=colors[colData$Condition], cex=1.2)


####C:Firstly optimize for the choice of number of unwanted variance factors and empirical genes, then decide the best parameters for RUVSeq normalization and add to DESeq2

L=as.numeric(c(0.8,0.85,0.875))
for (l in L){
    controlRes <- subset(top , FDR > l)
    empirical <- rownames(top)[which(rownames(set) %in% rownames(controlRes))]
    K<-5
    res2<-rep(0,K)
    seqRUVg=rep(0,K)
    for (k in 1:K){
        print(k)
        seqRUVg<- RUVg(as.matrix(filtered.counts),empirical,k=k)
        V=data.frame(colData$Condition,colData$gender,colData$ageF)
        col2<-cbind(V,seqRUVg$W)
        dds2 <- DESeqDataSetFromMatrix(countData = filtered.counts,colData=col2,design=~colData.Condition)
        
        design(dds2)=as.formula(paste(c("~",paste(colnames(seqRUVg$W),rep("+" ,k)),c(paste("colData.ageF","+","colData.gender","+","colData.Condition")))))
        dds2 <- DESeq(dds2)
        res2 <- results(dds2)
        summary(res2)
    }
    print(summary(res2))
}
print(summary(res2))

#Use the best combination of L and k to go back to DEseq2 as above
setk <- RUVg(set, empirical, k=k)
pData(setk)
plotRLE(setk, outline=FALSE, ylim=c(-4, 4), col=colors[colData$Condition])
plotPCA(setk, col=colors[colData$Condition], cex=1.2)


