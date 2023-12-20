rm(list = ls())
setwd("../火山图_热图/")

library(RColorBrewer)
df <- read.csv("../Deseq2/DESeq2_res2.csv",row.names = 1) #导入数据，第一列作为行名
log2FoldChange <- 1 #设置foldchange阈值
cut.log2FoldChange <- 1
pvalue <- 0.05 #设置p阈值
pdf( "df_volcano.pdf") #打开画板
plot(df$log2FoldChange, -log10(df$padj), col="#00000033", pch=19,
     xlab=paste("log2 (fold change)"),
     ylab="-log10 (padj)")
#筛选上下调
up <- subset(df, df$padj< 0.05 & df$log2FoldChange > 1)
down <- subset(df, df$padj< 0.05 & df$log2FoldChange< -1)
#绘制上下调
points(up$log2FoldChange, -log10(up$padj), col=1, bg = brewer.pal(9, "YlOrRd")[6],
       pch=21, cex=1)
points(down$log2FoldChange, -log10(down$padj), col = 1, bg = brewer.pal(11,"RdBu")
       [9], pch = 21,cex=1)
#加上线p、log2FoldChange阈值线
abline(h=-log10(pvalue),v=c(-1*log2FoldChange,log2FoldChange),lty=2,lwd=1)
dev.off()#关闭






####################画热图----------------------------------
rm(list = ls())
rt<-read.table('../Deseq2//Deseq2_rld.txt',sep = '\t',header = T,row.names = 1)

df <- read.csv("../Deseq2//sig DESeq22.csv",row.names = 1)
ex_b_limma3 <- rt[match(rownames(df),rownames(rt)),]#将差异基因矩阵与表达矩阵匹配，得到差异基因表达矩阵
library(pheatmap)
head(ex_b_limma3)

clinical<-read.csv('../Deseq2/clinical.csv')


gse<-colnames(rt)#获得样本名
group_list <- c(clinical$group)


n=t(scale(t(ex_b_limma3)))#用scale函数对数据标准化
n[n>2]=2
n[n< -2]= -2
n[1:4,1:4]

x=df$log2FoldChange
names(x)=rownames(df)
cg=c(names(head(sort(x),30)),names(tail(sort(x),30)))#差异前30和后30个基因
cg=names(tail(sort(abs(x)),10))


group2<-clinical$group[order(clinical$group)]
n<-n[,order(clinical$group)]



#pheatmap(n,show_colnames =F,show_rownames = F)#可以画前几个
annotation_col = data.frame(group=factor(group2))
rownames(annotation_col) = colnames(n) #将ac的行名也就分组信息 给到n的列名，即热图中位于上方的分组信息，这步很重要
pheatmap(n[cg,]  ,show_colnames =F,
         show_rownames = T,
         cluster_cols = F,
         annotation_col=annotation_col,
         fontsize = 8,
         filename = 'df-heatmap.pdf')

