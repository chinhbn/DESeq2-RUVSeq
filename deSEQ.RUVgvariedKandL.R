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

             
toLatex(sessionInfo())
\begin{itemize}\raggedright
  \item R version 3.4.2 (2017-09-28), \verb|x86_64-apple-darwin15.6.0|
  \item Locale: \verb|C|
  \item Running under: \verb|macOS Sierra 10.12.6|
  \item Matrix products: default
  \item BLAS: \verb|/Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRblas.0.dylib|
  \item LAPACK: \verb|/Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRlapack.dylib|
  \item Base packages: base, datasets, grDevices, graphics, methods,
    stats, utils
  \item Loaded via a namespace (and not attached):
    AnnotationDbi~1.38.2, Biobase~2.36.2, BiocGenerics~0.22.1,
    BiocParallel~1.10.1, Biostrings~2.44.2, DBI~0.7, DESeq2~1.16.1,
    DelayedArray~0.2.7, Formula~1.2-2, GenomeInfoDb~1.12.3,
    GenomeInfoDbData~0.99.0, GenomicAlignments~1.12.2,
    GenomicFeatures~1.28.5, GenomicRanges~1.28.6, Hmisc~4.0-3,
    IRanges~2.10.5, Matrix~1.2-11, RColorBrewer~1.1-2, RCurl~1.95-4.8,
    RSQLite~2.0, Rcpp~0.12.13, Rsamtools~1.28.0, S4Vectors~0.14.7,
    SummarizedExperiment~1.6.5, XML~3.98-1.9, XVector~0.16.0,
    acepack~1.4.1, annotate~1.54.0, backports~1.1.1, base64enc~0.1-3,
    biomaRt~2.32.1, bit~1.1-12, bit64~0.9-7, bitops~1.0-6, blob~1.1.0,
    checkmate~1.8.5, cluster~2.0.6, colorspace~1.3-2, compiler~3.4.2,
    data.table~1.10.4-2, digest~0.6.12, foreign~0.8-69,
    genefilter~1.58.1, geneplotter~1.54.0, ggplot2~2.2.1, grid~3.4.2,
    gridExtra~2.3, gtable~0.2.0, htmlTable~1.9, htmltools~0.3.6,
    htmlwidgets~0.9, knitr~1.17, lattice~0.20-35, latticeExtra~0.6-28,
    lazyeval~0.2.0, locfit~1.5-9.1, magrittr~1.5, matrixStats~0.52.2,
    memoise~1.1.0, munsell~0.4.3, nnet~7.3-12, pROC~1.10.0,
    parallel~3.4.2, plyr~1.8.4, rlang~0.1.2, rpart~4.1-11,
    rtracklayer~1.36.6, scales~0.5.0, splines~3.4.2, stats4~3.4.2,
    stringi~1.1.5, stringr~1.2.0, survival~2.41-3, tibble~1.3.4,
    tools~3.4.2, xtable~1.8-2, zlibbioc~1.22.0
\end{itemize}
