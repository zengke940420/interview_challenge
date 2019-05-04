rm(list=ls())
library("pheatmap")
library("DESeq2")
countdata<- read.table("/Users/2kisa/Dr.Shen_challenge/edited.txt", header = T, row.names = 1)
countdata<- as.matrix(countdata)
condition<-factor(c(rep("Nac",3),rep("PFC",3)))
colData<- data.frame(row.names = colnames(countdata), condition)
dds<- DESeqDataSetFromMatrix(countData = countdata, colData = colData, design = ~condition)
select <- order(rowMeans(counts(dds,normalized=FALSE)),decreasing=TRUE)
rld<- rlogTransformation(dds, fitType="mean")
plotPCA(rld, intgroup = "condition")
dds <- DESeq(dds)
res <- results(dds)
res <- res[order(res$padj), ]
deCount <- table(res$padj<0.05)
write.csv(res, file="difExp_results.csv")
png("DE_pvals.png", 1000, 1000, pointsize=20)
hist(res$pvalue, breaks=50, col="grey")
library("genefilter")
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),20)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(rld))
pheatmap(mat, annotation_col=df)

DE_genes <- as.integer(res$padj <= 0.05)
names(DE_genes) <- rownames(res)

library(goseq)
pwf=nullp(DE_genes[!is.na(DE_genes)],"mm9","ensGene")
#biocLite("org.Mm.eg.db")
GO.wall=goseq(pwf,"mm9","ensGene")

GO.samp=goseq(pwf,"mm9","ensGene",method="Sampling",repcnt=1000)
plot(log10(GO.wall[,2]), log10(GO.samp[match(GO.wall[,1],GO.samp[,1]),2]),
     + xlab="log10(Wallenius p-values)",ylab="log10(Sampling p-values)",
     + xlim=c(-3,0))
abline(0,1,col=3,lty=2)
GO.MF=goseq(pwf,"mm9","ensGene",test.cats=c("GO:MF"))
enriched.GO=GO.wall$category[p.adjust(GO.wall$over_represented_pvalue,method="BH")<.05]
library(GO.db)
for(go in enriched.GO[1:10]){
  + print(GOTERM[[go]])
  + cat("--------------------------------------\n")}
en2eg=as.list(org.Mm.egENSEMBL2EG)
eg2kegg=as.list(org.Mm.egPATH)
grepKEGG=function(id,mapkeys){unique(unlist(mapkeys[id],use.names=FALSE))}
kegg=lapply(en2eg,grepKEGG,eg2kegg)
KEGG=goseq(pwf,gene2cat=kegg)
write.table(KEGG, file="kegg.txt",sep = "\t",row.names = F,col.names = T)
